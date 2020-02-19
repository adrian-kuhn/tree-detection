#!/usr/bin/env python
# -*- coding: utf-8 -*-

import arcpy
import os
import logging
import logging.config
import logging.handlers
from multiprocessing import Pool
import datetime
import itertools

from . import settings, tiles, exporter, assignments

log = logging.getLogger()


class Controller:
    def __init__(self):
        """
        Init the controller by initializing an empty set of assignments
        """

        self.assignments = []
        """ List of all assignments to process """

        self.exporter = exporter.TreeExporter()
        """ Instance of a tree exporter used to export the trees to a file geo database after """

    def start_detection(self):
        """
        Main entry point to start the tree detection runner.
        Prepares all assignments and tiles, starts the process (serially or in parallel) end exports the detected trees.
        """
        if self.ready():
            start_time = datetime.datetime.now()
            self._load_assignments()
            self.assignments = [ass for ass in self.assignments if ass.tile_number == 1150412]
            self._run_serial() if settings.get("processing.processes") == 1 else self._run_parallel()
            self.exporter.to_gdb()
            self.exporter.clip_boundary()
            end_time = datetime.datetime.now()
            log.info("Extracted trees in {}".format(end_time - start_time))

    def ready(self):
        ready = True
        for p in settings.get("data.las"):
            if not os.path.isdir(os.path.dirname(p)):
                log.error("Processor not ready. {} not existing".format(os.path.dirname(p)))
                ready = False
        for p in settings.get("data.dom"):
            if not os.path.isdir(os.path.dirname(p)):
                log.error("Processor not ready. {} not existing".format(os.path.dirname(p)))
                ready = False
        for p in settings.get("data.dtm"):
            if not os.path.isdir(os.path.dirname(p)):
                log.error("Processor not ready. {} not existing".format(os.path.dirname(p)))
                ready = False
        if not arcpy.Exists(os.path.join(settings.get("data.tiling.gdb"), settings.get("data.tiling.feature_class"))):
            log.error("Processor not ready. {} not existing".format(os.path.join(settings.get("data.tiling.gdb"), settings.get("data.tiling.feature_class"))))
            ready = False
        return self.exporter.ready() and ready

    def get_assignment(self, number):
        """
        Get a single assignment by its tile number.

        :param number: The tile number
        :type number: Integer
        :return: The assignment for the asked tile number.
        :rtype: :class:`treedetection.assignments.Assignment`. None in case of nonexistent tile.
        """
        return next((a for a in self.assignments if a.tile_number == number), None)

    def _load_assignments(self):
        """
        Connect to the feature class with all tiles and init an assignment per tile.
        Calculate neighborhood for all tiles and assign the neighboring tiles to the assignment.
        The neighbourhood table association will be stored in memory and deleted after the extraction.
        The position of the neighbors are calculated with fixed numeric calculations according the tile number.
        """
        self.tiling_gdb = settings.get("data.tiling.gdb")
        self.tiling_fc = os.path.join(self.tiling_gdb, settings.get("data.tiling.feature_class"))
        self.tiling_field = settings.get("data.tiling.field_name")

        polygons = arcpy.da.SearchCursor(self.tiling_fc, [self.tiling_field, "SHAPE@"])
        for polygon in polygons:
            self.assignments.append(assignments.Assignment(polygon[0], settings.get("algorithm.tile_buffer"), tiles.Extent(extent=polygon[1].extent)))
        del polygons

        log.info("Found {} tiles in {}".format(len(self.assignments), self.tiling_fc))

        nbr_table_name = "in_memory\\neighbours"

        arcpy.PolygonNeighbors_analysis(self.tiling_fc, nbr_table_name, self.tiling_field, both_sides="BOTH_SIDES")
        log.info("Calculated polygon neighbourhood in {} and wrote result to {}".format(self.tiling_fc, nbr_table_name))

        top_left_diffs = settings.get("data.tiling.top_left_diffs")
        top_right_diffs = settings.get("data.tiling.top_right_diffs")
        bottom_left_diffs = settings.get("data.tiling.bottom_left_diffs")
        bottom_right_diffs = settings.get("data.tiling.bottom_right_diffs")

        edge_length_north = settings.get("data.tiling.edge_length_north")
        edge_length_east = settings.get("data.tiling.edge_length_east")

        cursor = arcpy.da.SearchCursor(nbr_table_name, ["src_{}".format(self.tiling_field), "nbr_{}".format(self.tiling_field), "LENGTH"])
        for row in cursor:
            src = row[0]
            nbr = row[1]
            length = row[2]
            assignment = self.get_assignment(src)
            nbr_assignment = self.get_assignment(nbr)

            if length == edge_length_east and nbr > src:
                assignment.set_neighbour(nbr, nbr_assignment.extent, tiles.Neighborhood.RIGHT)
            elif length == edge_length_east and nbr < src:
                assignment.set_neighbour(nbr, nbr_assignment.extent, tiles.Neighborhood.LEFT)
            elif length == edge_length_north and nbr > src:
                assignment.set_neighbour(nbr, nbr_assignment.extent, tiles.Neighborhood.BOTTOM_CENTER)
            elif length == edge_length_north and nbr < src:
                assignment.set_neighbour(nbr, nbr_assignment.extent, tiles.Neighborhood.TOP_CENTER)
            elif src - nbr in top_left_diffs:
                assignment.set_neighbour(nbr, nbr_assignment.extent, tiles.Neighborhood.TOP_LEFT)
            elif src - nbr in top_right_diffs:
                assignment.set_neighbour(nbr, nbr_assignment.extent, tiles.Neighborhood.TOP_RIGHT)
            elif src - nbr in bottom_left_diffs:
                assignment.set_neighbour(nbr, nbr_assignment.extent, tiles.Neighborhood.BOTTOM_LEFT)
            elif src - nbr in bottom_right_diffs:
                assignment.set_neighbour(nbr, nbr_assignment.extent, tiles.Neighborhood.BOTTOM_RIGHT)
            else:
                raise AttributeError("Cannot calculate neighborhood of tile with number {}".format(src))

        del cursor
        arcpy.Delete_management(nbr_table_name)
        log.debug("Deleted in memory feature class {}".format(nbr_table_name))

    def _run_parallel(self):
        """
        Private method to run the tree detection with multiple processes.
        The number of parallel processes is defined in the settings file.
        """
        log.info("Run parallel with {} processes".format(settings.get("processing.processes")))
        pool = Pool(settings.get("processing.processes"))
        self.assignments = pool.starmap(self._execute, zip(self.assignments, itertools.repeat(settings.get("logging"))))

    def _run_serial(self):
        """
        Private method to run the tree detection serially (with one single process)
        """
        log.info("Run serially")
        counter = 1
        for assignment in self.assignments:
            log.info("Processing {} of {} assignments ({})".format(counter, len(self.assignments), assignment))
            assignment.run()
            counter += 1

    @staticmethod
    def _execute(assignment, log_config):
        """
        Private static method to run an assignment in a separate process.
        Used in case of multiple processes as pool.map method parameter.
        Has to init an own logger, because the logger wouldn't be shared across multiple processes.
        Changes the name of the logfile by append its unique process id to prevent
        race conditions while writing to the same log file.

        :param assignment: One single assignment to process
        :type :class:`treedetection.assignments.Assignment`
        :param log_config: A dictionary with a logger configuration.
        :type log_config: Dictionary
        :return: The processed assignment
        :rtype: :class:`treedetection.assignments.Assignment`
        """
        if str(os.getpid()) not in log_config["handlers"]["file"]["filename"]:
            log_config["handlers"]["file"]["filename"] = log_config["handlers"]["file"]["filename"].replace(".log", "_{}.log".format(os.getpid()))
            logging.config.dictConfig(log_config)
        assignment.run()
        return assignment
