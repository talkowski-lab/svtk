#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

from setuptools import setup

setup(
    name='svtools',
    version='0.1',
    description='Tools for manipulating structural variant files',
    author='Matthew Stone',
    author_email='mstone5@mgh.harvard.edu',
    packages=['svtools'],
    package_data={'svtools': ['data/standard_template.vcf',
                              'data/vcfcluster_template.vcf']},
    scripts=[
        'scripts/standardize_vcf',
        'scripts/bedcluster',
        'scripts/vcfcluster'],
    install_requires=[
        'numpy',
        'scipy',
        'pysam',
    ]
)
