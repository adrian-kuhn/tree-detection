#!/usr/bin/env python
# -*- coding: utf-8 -*-

import arcpy
import os
import logging
import numpy
import laspy
import copy
from skimage import io
from enum import Enum

log = logging.getLogger()


class Neighborhood(Enum):
    """
    Enumeration for tile positions in a cartesian two dimensional 8-neighborhood.
    """
    CENTER = 0
    TOP_LEFT = 1
    TOP_CENTER = 2
    TOP_RIGHT = 3
    LEFT = 4
    RIGHT = 5
    BOTTOM_LEFT = 6
    BOTTOM_CENTER = 7
    BOTTOM_RIGHT = 8


class Tile:
    """
    Representation of a tile in a raster area.
    """
    __slots__ = ["number", "name", "path", "resolution", "position", "target_extent", "neighbours", "timberline", "i", "extent", "image"]
    """ Using slots for memory performance boost See: http://book.pythontips.com/en/latest/__slots__magic.html"""

    def __init__(self, number, name, path, resolution, position, target_extent, timberline):
        """
        Constructor for a tile.

        :param number: The unique tile number
        :type number: Integer
        :param name: The name of this tile
        :type name: String
        :param path: The path to the real world file containing the tile data (e.g. a TIFF, LAS file)
        :type path: String
        :param resolution: The resolution factor to 1 meter. E.g. raster data with 0.25m per pixel get the factor 4.
        :type resolution: Number
        :param position: The position of this tile in the cartesian 8-neighborhood
        :type position: :class:`treedetection.tiles.Neighborhood`
        :param target_extent: The target extent for this tile
        :type target_extent: :class:`treedetection.tiles.Extent`
        :param timberline: The timberline to fix false classified LiDAR points
        :type timberline: Integer
        """
        self.number = number
        self.name = name
        self.path = path
        self.resolution = resolution
        self.position = position
        self.target_extent = target_extent
        self.neighbours = []
        self.timberline = timberline

        self.i = 0
        """ Index variable to make the tile and neighbouring tiles iterable """

        self.extent = None
        self.image = None

    def __repr__(self):
        """
        String representation of this object

        :return: The tile as a string
        """
        return "{} tile {}".format(self.name, self.number)

    def __eq__(self, other):
        """
        Equality check for two tiles. Tiles are identical by unique number.

        :param other: The tile to check for identity.
        :type other: :class:`treedetection.runner.Tile`
        :return: True, if the other tile is identical with this tile, false otherwise
        :rtype: Boolean
        """
        return isinstance(other, self.__class__) and self.number == other.number

    def __len__(self):
        """
        Return the number of neighbouring tiles

        :return: Number of tiles as Integer
        """
        return len(self.neighbours)

    def __contains__(self, tile):
        """
        Checks, if the provided tile is a neighbour of this tile

        :param tile: A tile object
        :type tile: :class:`treedetection.runner.Tile`
        :return: True, if the tile is a neighbour, false otherwise
        """
        return tile in self.neighbours

    def __iter__(self):
        """
        Make the tile object iterable by walking through all neighbouring tiles

        :return: The tile itself
        """
        self.i = 0
        return self

    def __next__(self):
        """
        Get the next tile of the neighbour collection

        :return: The next tile collection
        :rtype: :class:`treedetection.runner.Tile`
        """
        if self.i < len(self.neighbours):
            next_tile = self.neighbours[self.i]
            self.i += 1
            return next_tile
        else:
            raise StopIteration()

    def exists(self):
        """
        Check existence of this tile on the file system

        :return: True, if the according file exists, false otherwise
        :rtype: Boolean
        """
        return os.path.exists(self.path)

    def append(self, neighbour):
        """
        Append a tile as neighbour to this tile if not already appended.

        :param neighbour: A new neighbouring tile
        :type neighbour: :class:`treedetection.runner.Tile`
        """
        if not self.__contains__(neighbour):
            if isinstance(neighbour, Tile):
                self.neighbours.append(neighbour)
            else:
                raise AttributeError("Only tile objects can be appended to another tile")

    def is_loaded(self):
        """
        Check if a tile is already loaded into a numpy array.

        :return: True, if tile is loaded. False otherwise.
        :rtype: Boolean
        """
        return self.image is not None

    def is_adapted(self):
        """
        Check if a tile is already adapted to target extent.

        :return: True, if tile is adapted. False otherwise.
        :rtype: Boolean
        """
        return self.extent == self.target_extent

    def is_center(self):
        """
        Check if a tile is a center tile.

        :return: True, if tile is a center tile. False otherwise.
        :rtype: Boolean
        """
        return self.position is Neighborhood.CENTER

    def prepare(self, buffer):
        """
        Preparing workflow with loading, adapting, handling of no data areas, clipping, preparing all neighbour tiles
        and padding.

        :param buffer: The buffer to pad this tile at the end of the preparing workflow.
        """
        self._load()
        self._adapt()
        self._handle_no_data()
        self._clip(buffer)

        for neighbour in self.neighbours:
            neighbour.prepare(buffer)

        self._pad(buffer)

    def _load(self):
        """
        Loading a tile must be implemented in concrete sub class.
        """
        raise NotImplementedError("Abstract method must be implemented in sub class.")

    def _adapt(self):
        """
        Adapt a tile to target extent by applying zero values.
        """
        if self.is_loaded() and not self.is_adapted():
            margin_left = self._get_margin_left()
            margin_top = self._get_margin_top()
            new_img = numpy.zeros((self.target_extent.height * self.resolution, self.target_extent.width * self.resolution))
            new_img[margin_top:margin_top + self.image.shape[0], margin_left:margin_left + self.image.shape[1]] = self.image

            # plt.imsave(arr=new_img, cmap=plt.get_cmap('Greys_r'), fname='field_{}_{}'.format(self.name, self.number))

            self.image = new_img
            self.extent = copy.deepcopy(self.target_extent)
            log.debug("Adapted {} to shape {}".format(self, self.image.shape))

    def _clip(self, buffer):
        """
        Clipping this tile according the position and buffer size.
        :param buffer: Buffer in meter
        :type buffer: Integer
        """
        if self.is_loaded() and not self.is_center():
            buffer = buffer * self.resolution
            clipped_image = self.image
            if self.position is Neighborhood.TOP_CENTER:
                clipped_image = clipped_image[-buffer:]
            elif self.position is Neighborhood.BOTTOM_CENTER:
                clipped_image = clipped_image[:buffer]
            elif self.position is Neighborhood.LEFT:
                clipped_image = clipped_image[:, -buffer:]
            elif self.position is Neighborhood.RIGHT:
                clipped_image = clipped_image[:, :buffer]
            elif self.position is Neighborhood.TOP_LEFT:
                clipped_image = clipped_image[-buffer:, -buffer:]
            elif self.position is Neighborhood.TOP_RIGHT:
                clipped_image = clipped_image[-buffer:, :buffer]
            elif self.position is Neighborhood.BOTTOM_LEFT:
                clipped_image = clipped_image[:buffer, -buffer:]
            elif self.position is Neighborhood.BOTTOM_RIGHT:
                clipped_image = clipped_image[:buffer, :buffer]

            self.image = clipped_image.copy()
            log.debug("Clipped {} to {} position with shape {}".format(self, self.position, self.image.shape))
            # Todo: Fix extent

    def _handle_no_data(self):
        """
        Handle no data areas by setting values to 0.
        """
        if self.is_loaded():
            self.image[self.image < 0] = 0
            log.debug("Set no data areas to 0 in {}".format(self))

    def _pad(self, buffer):
        """
        Padding the tile according the given buffer and neighboring tiles.

        :param buffer: The buffer in meter
        :type buffer: Integer
        """
        if self.is_loaded() and self.is_center():
            pixel_buffer = buffer * self.resolution
            expanded_img = numpy.pad(self.image, pixel_buffer, 'constant', constant_values=0)
            for neighbour in self.neighbours:
                if neighbour.exists():
                    margin_left = 0
                    margin_top = 0
                    if neighbour.position is Neighborhood.TOP_CENTER:
                        margin_left = pixel_buffer
                        margin_top = 0
                    elif neighbour.position is Neighborhood.BOTTOM_CENTER:
                        margin_left = pixel_buffer
                        margin_top = expanded_img.shape[0] - pixel_buffer
                    elif neighbour.position is Neighborhood.LEFT:
                        margin_left = 0
                        margin_top = pixel_buffer
                    elif neighbour.position is Neighborhood.RIGHT:
                        margin_left = expanded_img.shape[1] - pixel_buffer
                        margin_top = pixel_buffer
                    elif neighbour.position is Neighborhood.TOP_LEFT:
                        margin_left = 0
                        margin_top = 0
                    elif neighbour.position is Neighborhood.TOP_RIGHT:
                        margin_left = expanded_img.shape[1] - pixel_buffer
                        margin_top = 0
                    elif neighbour.position is Neighborhood.BOTTOM_LEFT:
                        margin_left = 0
                        margin_top = expanded_img.shape[0] - pixel_buffer
                    elif neighbour.position is Neighborhood.BOTTOM_RIGHT:
                        margin_left = expanded_img.shape[1] - pixel_buffer
                        margin_top = expanded_img.shape[0] - pixel_buffer

                    height = neighbour.image.shape[0]
                    width = neighbour.image.shape[1]
                    expanded_img[margin_top:margin_top + height, margin_left:margin_left + width] = neighbour.image

            self.image = expanded_img
            self.extent.pad(buffer)
            log.debug("Padded {} with neighbouring tiles to shape {}".format(self, self.image.shape))

            self._clear_neighbours()

    def _clear_neighbours(self):
        """
        Delete neighboring tiles. Used for memory performance boost after padding this tile.
        """
        del self.neighbours
        log.debug("Removed neighbouring tiles from memory in {}".format(self))

    def _get_margin_top(self):
        """
        Calculation of margins must be implemented in concrete sub class.
        """
        raise NotImplementedError("Abstract method must be implemented in sub class.")

    def _get_margin_left(self):
        """
        Calculation of margins must be implemented in concrete sub class.
        """
        raise NotImplementedError("Abstract method must be implemented in sub class.")


class RasterTile(Tile):
    def _load(self):
        """
        Load a raster tile from file system into a numpy array. Calculate the extent afterwards.
        """
        if self.exists():
            self.image = io.imread(self.path)
            self.extent = Extent(extent=arcpy.Describe(self.path).extent)
            log.debug("Loaded {}".format(self))

    def _get_margin_top(self):
        """
        Get the margin top for adapting a raster tile.
        :return: The margin in number of pixel
        :rtype: Integer
        """
        return int(round((self.target_extent.max_y - self.extent.max_y) * self.resolution))

    def _get_margin_left(self):
        """
        Get the margin left for adapting a raster tile.
        :return: The marging in number of pixel
        :rtype: Integer
        """
        return int(round((self.extent.min_x - self.target_extent.min_x) * self.resolution))


class LasTile(Tile):
    """
    Represents a tile containing a LiDAR point cloud (*.las)
    """

    def _load(self):
        """
        Load a las tile with `laspy` and filter to high vegetation
        Calculate a 2D histogram by projecting the point cloud to the ground.
        To be compatible with the DOM and DTM, the numpy array has to be rotated by -90 degree.
        """
        if self.exists():
            las = laspy.file.File(self.path, mode="r")
            veg = numpy.where(numpy.logical_and(las.Classification == 5, las.z < self.timberline))
            veg_points = las.points[veg]
            veg_x = veg_points["point"]["X"] * las.header.scale[0] + las.header.offset[0]
            veg_y = veg_points["point"]["Y"] * las.header.scale[1] + las.header.offset[1]
            las.close()

            min_x = numpy.min(veg_x) if len(veg_x) > 0 else 0
            max_x = numpy.max(veg_x) if len(veg_x) > 0 else 0
            min_y = numpy.min(veg_y) if len(veg_y) > 0 else 0
            max_y = numpy.max(veg_y) if len(veg_y) > 0 else 0
            self.extent = Extent(min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)

            bin_x = self.extent.width
            bin_y = self.extent.height

            if bin_x > 0 and bin_y > 0:
                his, bx, by = numpy.histogram2d(veg_x, veg_y, bins=[(bin_x * self.resolution), (bin_y * self.resolution)])
                self.image = numpy.rot90(his)
            else:
                self.image = numpy.zeros((0, 0))
            log.debug("Loaded {}, filtered to high vegetation and calculated a 2D histogram".format(self))
            del las, veg, veg_points, veg_x, veg_y

    def _get_margin_top(self):
        """
        Get the margin top for adapting a las tile.
        :return: The margin in number of pixel
        :rtype: Integer
        """
        return int(round((self.target_extent.max_y - self.extent.max_y))) * self.resolution

    def _get_margin_left(self):
        """
        Get the margin left for adapting a las tile.
        :return: The margin in number of pixel
        :rtype: Integer
        """
        return int(round((self.extent.min_x - self.target_extent.min_x))) * self.resolution


class Extent:
    """
    Represents an extent. Used as wrapper class for ESRI extents.
    """
    def __init__(self, extent=None, min_x=None, max_x=None, min_y=None, max_y=None):
        """
        Init an extent with an ESRI extent object or explicit edge values.
        :param extent: An ESRI extent object
        :type extent: ESRI extent
        :param min_x: X-Coordinate of the lower left corner
        :type min_x: Integer
        :param max_x: X-Coordinate of the upper right corner
        :type max_x: Integer
        :param min_y: Y-Coordinate of the lower left corner
        :type min_y: Integer
        :param max_y: Y-Coordinate of the upper right corner
        :type max_y: Integer
        """
        if extent:
            min_x = extent.XMin
            max_x = extent.XMax
            min_y = extent.YMin
            max_y = extent.YMax

        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y
        self.width = int(round(max_x - min_x))
        self.height = int(round(max_y - min_y))

    def __eq__(self, other):
        """
        Check equality of two extents. Two extents are equal when width and size are equal

        :param other: The extent to compare this object with
        :type other: :class:`runner.Extent`
        :return: True, if this extent represent the same extent, false otherwise.
        """
        return isinstance(other, self.__class__) and self.width == other.width and self.height == other.height

    def _calculate_extent(self):
        """
        Recalculate the width and height of the extent.
        """
        self.width = int(round(self.max_x - self.min_x))
        self.height = int(round(self.max_y - self.min_y))

    def pad(self, buffer):
        """
        Pad the extent with a given buffer. Recalculate the extent width and height afterwards.
        :param buffer: The buffer in meter
        :type buffer: Float
        """
        self.min_x -= buffer
        self.max_x += buffer
        self.min_y -= buffer
        self.max_y += buffer
        self._calculate_extent()
