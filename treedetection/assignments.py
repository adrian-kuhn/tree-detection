#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
import csv
import os

from skimage.feature import peak_local_max
from skimage.morphology import label, closing, disk
from skimage.segmentation import watershed
from skimage.measure import regionprops
from skimage.filters import gaussian
from skimage.restoration import denoise_tv_chambolle

from . import settings, tiles

log = logging.getLogger()


class Assignment:
    """
    Representation of an assignment containing all the data (las, dom, dtm) to delineate the individual trees.
    """
    __slots__ = ["tile_number", "tile_buffer", "extent", "params", "resolution", "csv_file", "timberline", "las", "dom",
                 "dtm"]
    """ Using slots for memory performance boost. See: http://book.pythontips.com/en/latest/__slots__magic.html"""

    def __init__(self, tile_number, tile_buffer, extent):
        """
        Initialization of an assignment.

        :param tile_number: The unique tile number
        :type tile_number: Integer
        :param tile_buffer: The buffer in meter to extend this assignment before processing
        :type tile_buffer: Integer
        :param extent: The extent of the tiles belonging to this assignment
        :type extent: :class:`treedetection.tiles.Extent`
        """
        self.tile_number = tile_number
        self.tile_buffer = tile_buffer
        self.extent = extent
        self.params = WatershedParams()
        self.resolution = self.params.resolution
        self.csv_file = settings.get("data.tmp_trees").format(tile_number=self.tile_number)
        self.timberline = settings.get("tree_restrictions.timberline")

        las_path = settings.get("data.las")[tile_number % len(settings.get("data.las"))].format(tile_number=tile_number)
        dom_path = settings.get("data.dom")[tile_number % len(settings.get("data.dom"))].format(tile_number=tile_number)
        dtm_path = settings.get("data.dtm")[tile_number % len(settings.get("data.dtm"))].format(tile_number=tile_number)
        self.las = tiles.LasTile(tile_number, "LAS", las_path, self.resolution, tiles.Neighborhood.CENTER, extent,
                                 self.timberline)
        self.dom = tiles.RasterTile(tile_number, "DOM", dom_path, self.resolution, tiles.Neighborhood.CENTER, extent,
                                    self.timberline)
        self.dtm = tiles.RasterTile(tile_number, "DTM", dtm_path, self.resolution, tiles.Neighborhood.CENTER, extent,
                                    self.timberline)

    def __repr__(self):
        """
        String representation of this object

        :return: The assignment as a string
        """
        return "Assignment {}".format(self.tile_number)

    def __eq__(self, other):
        """
        Check equality of two assignments. Two assignments are equal when processing the same tile number.

        :param other: The assignment to compare this object with
        :type other: :class:`runner.Assignment`
        :return: True, if the assignment represent the same tile number, false otherwise.
        """
        return isinstance(other, self.__class__) and self.tile_number == other.tile_number

    def set_neighbour(self, tile_number, extent, position):
        """
        Append neighboring tiles to this assignment.
        :param tile_number: The unique tile number
        :type tile_number: Integer
        :param extent: The extent of the tiles belonging to this assignment
        :type extent: :class:`treedetection.tiles.Extent`
        :param position: The position of this tile in the cartesian 8-neighborhood
        :type position: :class:`treedetection.tiles.Neighborhood`
        """
        las_path = settings.get("data.las")[tile_number % len(settings.get("data.las"))].format(tile_number=tile_number)
        dom_path = settings.get("data.dom")[tile_number % len(settings.get("data.dom"))].format(tile_number=tile_number)
        dtm_path = settings.get("data.dtm")[tile_number % len(settings.get("data.dtm"))].format(tile_number=tile_number)
        self.las.append(tiles.LasTile(tile_number, "LAS", las_path, self.resolution, position, extent, self.timberline))
        self.dom.append(
            tiles.RasterTile(tile_number, "DOM", dom_path, self.resolution, position, extent, self.timberline))
        self.dtm.append(
            tiles.RasterTile(tile_number, "DTM", dtm_path, self.resolution, position, extent, self.timberline))

    def ready(self):
        """
        Check for existence of LAS, DOM and DTM files to run this assignment
        :return: True, if LAS, DOM and DTM files are existing. False otherwise.
        :rtype: Boolean
        """
        return self.las.exists() and self.dom.exists() and self.dtm.exists()

    def clear_data(self):
        """
        Delete some data after delineation of trees to free some memory.
        """
        del self.las
        del self.dtm
        del self.dom
        del self.params
        log.debug("Removed LAS, DTM and DOM from memory in {}".format(self))

    def run(self):
        """
        Main entry point to run this assignment. Checks for existence of required data before.
        Prepare all the LAS-, DOM- and DTM-files, delineate the trees and export it to a CSV file.
        """
        if not self.ready():
            log.info("Skipping {}. Could not find all required data (LAS, DOM and DTM)".format(self))
        elif os.path.exists(self.csv_file):
            log.info("Skipping {}. Trees already delineated".format(self))
        else:
            log.info("Loading data for {}".format(self.tile_number))
            self.las.prepare(self.tile_buffer)
            self.dom.prepare(self.tile_buffer)
            self.dtm.prepare(self.tile_buffer)
            log.info("Delineating trees in {}".format(self.las.number))
            trees = self._delineate_trees()
            self.clear_data()
            self.to_csv(trees)

    def _delineate_trees(self):
        """
        Private function to get the single trees by watershed segmentation.
        Forest and open field areas are processed separately with distinct parameters.

        :return: A list of trees
        :rtype: Dictionary of :class:`treedetection.assignments.Tree`
        """
        # Closing
        closed = closing(self.las.image, disk(self.params.closing_radius))
        log.debug("Morphologically closed {}".format(self.las.number))

        # Create a mask for regions with trees
        mask = numpy.copy(closed)
        mask[mask != 0] = 1
        del closed

        veg_dom = numpy.ma.array(self.dom.image, mask=(1 - mask).astype(int), fill_value=0).filled()

        # Separating field from forest regions
        regions_field = label(mask)
        regions_forest = numpy.copy(regions_field)
        region_props = regionprops(regions_field, intensity_image=self.dtm.image)
        forest_labels = [r.label for r in region_props if
                         r.filled_area / (
                                 self.params.resolution * self.params.resolution) > self.params.forest_area_threshold or r.mean_intensity > self.params.conifer_height_threshold]
        regions_forest[numpy.isin(regions_forest, forest_labels, invert=True)] = 0
        regions_field[numpy.isin(regions_field, forest_labels)] = 0

        field = numpy.ma.array(veg_dom, mask=regions_forest, fill_value=0).filled()
        forest = numpy.ma.array(veg_dom, mask=regions_field, fill_value=0).filled()
        log.debug("Separated forest and field areas for {}".format(self.las.number))

        del veg_dom

        trees_field = self._watershed(field, self.las.number, "field", self.params.field_denoising_weight,
                                      self.params.field_sigma, self.params.field_truncate,
                                      self.params.field_min_distance, self.params.field_compactness)
        trees_forest = self._watershed(forest, self.las.number, "forest", self.params.forest_denoising_weight,
                                       self.params.forest_sigma, self.params.forest_truncate,
                                       self.params.forest_min_distance, self.params.forest_compactness)

        trees = trees_field + (trees_forest * (numpy.max(trees_field) + 1))

        del field
        del forest
        del trees_field
        del trees_forest

        log.info("Found {} trees in {}".format(len(regionprops(trees)), self.las.number))
        return self._extract_tree_params(trees)

    def _extract_tree_params(self, trees):
        """
        Private method to extract tree params from delineated trees.
        To get the shape the region properties of the trees has to be analyzed with the DOM.
        To get the height an additional analyzing loop with the DTM is required.

        Unrealistic trees will be eliminated according the settings.

        :param trees: Regions of all trees found in the given area.
        :type trees: numpy Array
        :return: Dictionary object containing all trees
        :rtype: Dictionary
        """
        found_trees = {}

        # Tree position (3D)
        for tree in regionprops(trees, intensity_image=self.dom.image):
            # Export with location of weighted_centroid
            centroid = list(tree.weighted_centroid)
            if self.params.tile_buffer * self.params.resolution < centroid[0] < self.dom.image.shape[0] - (
                    self.params.tile_buffer * self.params.resolution) and self.params.tile_buffer * self.params.resolution < \
                    centroid[1] < self.dom.image.shape[1] - (self.params.tile_buffer * self.params.resolution):
                x = centroid[1] / self.params.resolution + self.las.extent.min_x
                y = self.las.extent.max_y - centroid[0] / self.params.resolution
                z = tree.max_intensity

                loc_y, loc_x = numpy.where(tree.intensity_image == z)
                x_height = loc_x[0]
                y_height = loc_y[0]

                chm = tree.intensity_image[numpy.nonzero(tree.intensity_image)]
                minor_axis = max(1 / self.params.resolution, tree.minor_axis_length / self.params.resolution)
                major_axis = max(1 / self.params.resolution, tree.major_axis_length / self.params.resolution)

                found_trees[tree.label] = Tree(tree.label, x, y, z, x_height, y_height,
                                               dom_max=tree.max_intensity,
                                               dom_mean=numpy.mean(chm),
                                               dom_median=numpy.median(chm),
                                               dom_min=numpy.min(chm),
                                               major_axis=major_axis,
                                               minor_axis=minor_axis,
                                               area=tree.filled_area / (
                                                       self.params.resolution * self.params.resolution))

        for tree in regionprops(trees, intensity_image=self.dtm.image):
            if tree.label in found_trees:
                t = found_trees[tree.label]
                ohm = tree.intensity_image[numpy.nonzero(tree.intensity_image)]

                t.dtm_mean = numpy.mean(ohm)
                t.dtm_min = numpy.min(ohm)
                t.dtm_median = numpy.median(ohm)
                t.dtm_max = tree.max_intensity
                t.height = t.dom_max - tree.intensity_image[t.y_height, t.x_height]

                # Eliminate unrealistic trees according tree params
                if t.height < self.params.min_tree_height or \
                        t.major_axis <= self.params.min_major_axis or \
                        t.minor_axis <= self.params.min_minor_axis or \
                        round(t.minor_axis / t.major_axis, 2) <= self.params.min_eccentricity:
                    del found_trees[tree.label]
                    continue

                # Set some unrealistic values to null
                if t.height > self.params.max_tree_height:
                    t.height = None

                if t.major_axis > self.params.max_tree_diameter:
                    t.major_axis = None

        log.info("Keep {} trees in {} after analyzing the parameters".format(len(found_trees), self.las.number))

        del trees
        return found_trees

    @staticmethod
    def _watershed(img, number, region, denoising_weight, sigma, truncate, min_distance, compactness):
        """
        Private function to run the watershed segmentation for a prepared image.

        :param img: The image to run the watershed segmentation
        :type img: Numpy array
        :param number: The unique tile number used for logging
        :type number: Integer
        :param region: The region for with the watershed segmentation is executed. (Used for logging)
        :type region: String
        :param denoising_weight: The weight factor for denoising the image before watershed segmentation. See: https://scikit-image.org/docs/dev/api/skimage.restoration.html#skimage.restoration.denoise_tv_chambolle
        :type denoising_weight: Float
        :param sigma: Sigma value for gaussian blurring the image before watershed segmentation. See: https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian
        :type sigma: Float
        :param truncate: Truncation value for gaussian blurring the image before watershed segmentation. See: https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian
        :type truncate: Float
        :param min_distance: Minimum distance for local maxia representing tree tops.
        :type min_distance: Float
        :param compactness: Compactness of a watershed basin: See https://scikit-image.org/docs/dev/api/skimage.morphology.html#skimage.morphology.watershed
        :type compactness: Float
        :return: The image containing the labeled areas for individual trees.
        :rtype: Numpy array
        """
        if denoising_weight != 0:
            denoised = denoise_tv_chambolle(img, weight=denoising_weight)
            log.debug("Denoised {} area {}".format(region, number))
        else:
            denoised = img

        # Gauss filter
        gauss = gaussian(denoised, sigma, mode='constant', cval=0, preserve_range=True, truncate=truncate)
        log.debug("Gaussian smoothed {} area {}".format(region, number))

        # Create a mask so the segmentation will only occur on the trees
        mask = numpy.copy(img)
        mask[mask != 0] = 1

        # Local maxima
        # Get coordinates of peaks and convert to boolean mask (what indices=False used to do)
        coordinates = peak_local_max(gauss, min_distance=min_distance, exclude_border=False)
        local_max = numpy.zeros_like(gauss, dtype=bool)
        local_max[tuple(coordinates.T)] = True

        markers = label(local_max)
        log.debug("Peaked local maxima in {} area {}".format(region, number))

        # watershed
        labels = watershed(gauss, markers, mask=mask, compactness=compactness)
        log.debug("Applied watershed delineation in {} area {}".format(region, number))
        return labels

    def to_csv(self, trees):
        """
        Export delineated trees to CSV file.

        :param trees: A dictionary objects containing trees.
        :type trees: Dictionary with :class:`treedetection.assignments.Tree`
        """
        log.info("Exporting trees to CSV file {}".format(self.csv_file))
        with open(self.csv_file, 'w', newline='') as csv_file:
            for tree in trees.items():
                writer = csv.writer(csv_file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(tree[1].dump())


class WatershedParams:
    """
    Container for all watershed parameters.
    """

    def __init__(self):
        """
        All parameters are read from settings
        """
        self.tile_buffer = settings.get("algorithm.tile_buffer")
        self.resolution = int(1 / settings.get("algorithm.chm_resolution"))
        self.forest_area_threshold = settings.get("algorithm.forest_area_threshold")
        self.conifer_height_threshold = settings.get("algorithm.conifer_height_threshold")
        self.closing_radius = settings.get("algorithm.closing_radius")

        self.field_denoising_weight = settings.get("algorithm.field_denoising_weight")
        self.field_sigma = settings.get("algorithm.field_sigma")
        self.field_truncate = settings.get("algorithm.field_truncate")
        self.field_min_distance = settings.get("algorithm.field_min_distance")
        self.field_compactness = settings.get("algorithm.field_compactness")

        self.forest_denoising_weight = settings.get("algorithm.forest_denoising_weight")
        self.forest_sigma = settings.get("algorithm.forest_sigma")
        self.forest_truncate = settings.get("algorithm.forest_truncate")
        self.forest_min_distance = settings.get("algorithm.forest_min_distance")
        self.forest_compactness = settings.get("algorithm.forest_compactness")

        self.max_tree_height = settings.get("tree_restrictions.max_tree_height")
        self.max_tree_diameter = settings.get("tree_restrictions.max_tree_diameter")
        self.min_tree_height = settings.get("tree_restrictions.min_tree_height")
        self.min_major_axis = settings.get("tree_restrictions.min_major_axis")
        self.min_minor_axis = settings.get("tree_restrictions.min_minor_axis")
        self.min_eccentricity = settings.get("tree_restrictions.min_eccentricity")


class Tree:
    """
    Representation of a single tree
    """

    def __init__(self, tree_number, x=None, y=None, z=None, x_height=None, y_height=None, dtm_min=None, dtm_mean=None,
                 dtm_median=None, dtm_max=None,
                 dom_min=None, dom_mean=None, dom_median=None, dom_max=None, area=None, minor_axis=None,
                 major_axis=None, height=None):
        """
        Constructor for a single tree.

        :param tree_number: The unique tree number
        :type tree_number: Integer
        :param x: The x-coordinate of this tree (east value according LV95)
        :type x: Float
        :param y: The y-coordinate of this tree (north value according LV95)
        :type y: Float
        :param z: The z-value of this tree (height over sea level)
        :type z: Float
        :param x_height: The x-coordinate of this tree at the position of the maximum height
        :type x_height: Float
        :param y_height: The y-coordinate of this tree at the position of the maximum height
        :type y_height: Float
        :param dtm_min: The minimum value of the digital terrain model in the footprint area of this tree
        :type dtm_min: Float
        :param dtm_mean: The mean value of the digital terrain model in the footprint area of this tree
        :type dtm_mean: Float
        :param dtm_median: The median value of the digital terrain model in the footprint area of this tree
        :type dtm_mean: Float
        :param dtm_max: The maximum value of the digital terrain model in the footprint area of this tree
        :type dtm_max: Float
        :param dom_min: The minimum value of the digital surface model in the footprint area of this tree
        :type dom_min: Float
        :param dom_mean: The mean value of the digital surface model in the footprint area of this tree
        :type dom_mean: Float
        :param dom_median: The median value of the digital surface model in the footprint area of this tree
        :type dom_median: Float
        :param dom_max: The max value of the digital surface model in the footprint area of this tree
        :type dom_max: Float
        :param area: The area in square meters of the footprint belonging to this tree
        :type area: Float
        :param major_axis: The diameter of this tree on the widest position
        :type major_axis: Float
        :param minor_axis: The diameter of this tree on the narrowest position
        :type minor_axis: Float
        :param height: The relative height of this tree calculated by the difference of dom_max and dtm_mean
        :type height: Float
        """
        self.tree_number = tree_number
        self.x = x
        self.y = y
        self.z = z
        self.x_height = x_height
        self.y_height = y_height
        self.dtm_min = dtm_min
        self.dtm_mean = dtm_mean
        self.dtm_median = dtm_median
        self.dtm_max = dtm_max
        self.dom_min = dom_min
        self.dom_mean = dom_mean
        self.dom_median = dom_median
        self.dom_max = dom_max
        self.area = area
        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.height = height

    def __repr__(self):
        """
        String representation of this object

        :return: The tree data as a string
        """
        return "Tree {} ({}, {}, {})".format(self.tree_number, self.x, self.y, self.z)

    def dump(self):
        """
        Dump this tree into a list. Used for exporting to a CSV file.

        :return: All tree params in a list
        :rtype: List<String>
        """
        return [self.tree_number, self.major_axis, self.minor_axis, self.area, self.dtm_min, self.dtm_mean,
                self.dtm_median, self.dtm_max, self.dom_min, self.dom_mean, self.dom_median, self.dom_max,
                self.height, self.x, self.y, self.z]
