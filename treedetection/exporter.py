import os
import arcpy
import csv
import logging
import datetime
import glob

from . import settings
log = logging.getLogger()


class TreeExporter:
    """
    Generic class to export the result to a ESRI file geo database.
    """
    def __init__(self):
        self.gdb_name = settings.get("exporter.gdb.name")
        self.gdb_folder = settings.get("exporter.folder")
        self.gdb_path = os.path.join(self.gdb_folder, self.gdb_name)
        self.fc_name = settings.get("exporter.gdb.feature_class")
        self.fc_path = os.path.join(self.gdb_path, self.fc_name)

        self.csv_files = settings.get("data.tmp_trees")
        self.gdb_fields = settings.get("exporter.gdb.fields")

        self.canton_boundary = settings.get("data.canton_boundary")

    def ready(self):
        """
        Check for locks in exporting feature classes and check, if every path is accessible.
        :return: True, if the exporter is ready. False otherwise
        :rtype: Boolean
        """
        ready = True
        if not os.path.isdir(self.gdb_folder):
            log.error("Exporter not ready. {} is not existing.".format(self.gdb_folder))
            ready = False
        if not os.path.isdir(os.path.dirname(self.csv_files)):
            log.error("Exporter not ready. {} is not existing.".format(os.path.dirname(self.csv_files)))
            ready = False
        if not arcpy.Exists(self.canton_boundary):
            log.error("Exporter not ready. {} is not existing.".format(os.path.dirname(self.canton_boundary)))
            ready = False
        if self._is_locked():
            ready = False
            log.error("Exporter not ready. {} is locked".format(self.gdb_path))

        return ready

    def _is_locked(self):
        """
        Check for locks in the resulting file geo database.

        :return: True, if the file geo database is locked. False otherwise.
        :rtype: Boolean
        """
        lock_file = self._get_lock_file()
        if not lock_file:
            return False

        try:
            os.remove(lock_file)
            return False
        except OSError:
            return True

    def _get_lock_file(self):
        """
        Private method to search for lock files within a ESRI file geo database.

        :return: The name of the lock file in case of existence. An empty string otherwise.
        :rtype: String
        """
        full_file_paths = glob.glob(self.gdb_path + "\\*.lock")
        for f in full_file_paths:
            if f.endswith(".lock"):
                return f
        return ""

    def to_gdb(self):
        """
        Parse the folder with CSV files containing single trees per row.
        Write every tree as a point feature class to the file geo database and set a unique id for every tree.
        """
        log.debug("Start exporting trees to feature class {}".format(self.fc_path))
        self._init_gdb()
        self._init_fc()

        fields = [f["name"] for f in self.gdb_fields]
        fields.append("SHAPE@XYZ")
        cursor = arcpy.da.InsertCursor(self.fc_path, fields)

        counter = 0
        for root, dirs, files in os.walk(os.path.dirname(self.csv_files)):
            for name in files:
                if name.endswith(".csv"):
                    with open(os.path.join(root, name), 'r', newline='') as csv_file:
                        reader = csv.reader(csv_file, delimiter=',', quotechar='|')
                        for row in reader:
                            counter += 1
                            cursor.insertRow([counter, row[1] if row[1] else None, row[2], row[3], row[4],
                                              row[5], row[6], row[7], row[8], row[9], row[10], row[11],
                                              row[12] if row[12] else None,
                                              (float(row[13]), float(row[14]), float(row[15]))])

        log.info("Exported {} trees to {}".format(counter, self.fc_path))
        del cursor

    def _init_gdb(self):
        """
        Init the file geo database containing the resulting trees. Creates the configured folder, if not existing.
        """
        if not os.path.isdir(self.gdb_folder):
            os.mkdir(os.path.dirname(self.gdb_folder))

        if not arcpy.Exists(self.gdb_path):
            arcpy.CreateFileGDB_management(self.gdb_folder, self.gdb_name)
            log.info("Created fGDB {0}".format(self.gdb_name))

    def _init_fc(self):
        """
        Init the feature class for storing the trees. Trees will be stored in a point feature class with z value.

        Delete existing feature class if overwrite setting is set to true.
        Make the feature class unique by timestamp otherwise.
        """
        if arcpy.Exists(self.fc_path):
            if settings.get("exporter.overwrite"):
                arcpy.Delete_management(self.fc_path)
                log.info("Deleted feature class {0}".format(self.fc_name))
            else:
                self.fc_name = "{}_{}".format(self.fc_name, datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
                self.fc_path = os.path.join(self.gdb_path, self.fc_name)

        arcpy.CreateFeatureclass_management(self.gdb_path, self.fc_name, "POINT", has_z="Yes", spatial_reference=settings.get("exporter.gdb.spatial_reference"))
        log.info("Created point feature class {0}".format(self.fc_name))

        for field in self.gdb_fields:
            arcpy.AddField_management(self.fc_path, field["name"], field["type"])
            log.info("Added field {}({}) in {}".format(field["name"], field["type"], self.fc_name))

    def clip_boundary(self):
        """
        Clip resulting feature class to the boundary of the canton.
        Delete original feature class after clipping.
        """
        fc_path_clipped = "{}_clipped".format(self.fc_path)

        rows_before = int(arcpy.GetCount_management(self.fc_path).getOutput(0))
        arcpy.Clip_analysis(self.fc_path, self.canton_boundary, fc_path_clipped)
        rows_after = int(arcpy.GetCount_management(fc_path_clipped).getOutput(0))
        log.info("Clipped to canton boundary. Deleted {} trees.".format(rows_before - rows_after))

        arcpy.Delete_management(self.fc_path)
        arcpy.Rename_management(fc_path_clipped, self.fc_path)
