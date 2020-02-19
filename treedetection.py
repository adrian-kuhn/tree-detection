# !/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging.config

from treedetection import settings, processing


def main():
    """
    Main function checks for python 3, initialize the settings and logging and starts the tree detection.
    """
    if sys.version_info.major != 3:
        print('Please activate the conda package and run on python 3')
        return

    settings.init()
    logging.config.dictConfig(settings.get("logging"))
    r = processing.Controller()
    r.start_detection()


if __name__ == '__main__':
    main()
