#!/usr/bin/env python
# -*- coding: utf-8 -*-

import codecs
import os
import sys

from setuptools import find_packages, setup

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = "\n" + f.read()

about = {}

with open(os.path.join(here, "tree-detection", "__version__.py")) as f:
    exec(f.read(), about)

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist bdist_wheel upload")
    sys.exit()

required = [
    "arcgis > 1.6",
    "numpy >= 1.16.2",
    "scikit-image >= 0.16",
    "laspy >= 1.2"
]

setup(
    name="tree-detection",
    version=about["__version__"],
    description="CLI application to delineate single trees from LiDAR data, DOM and DTM",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Adrian Kuhn",
    author_email="adrian.kuhn@lu.ch",
    url="https://app.geo.lu.ch/redmine/projects/tree_detection",
    packages=find_packages(exclude=["tests"]),
    entry_points={},
    package_data={},
    python_requires=">=3.6",
    setup_requires=[],
    install_requires=required,
    extras_require={},
    include_package_data=True,
    license="Other",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: Other/Proprietary License",
        "Natural Language :: English",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Utilities"
    ],
    cmdclass={},
)
