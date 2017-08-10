#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
import pysam
import svtools.utils as svu


def vcf2bed(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('bed')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = pysam.VariantFile(args.vcf)
    bt = svu.vcf2bedtool(vcf, split_bnd=False, include_samples=True,
                         include_strands=False)

    if args.bed in 'stdout -'.split():
        sys.stdout.write(str(bt))
    else:
        bt.saveas(args.bed)
