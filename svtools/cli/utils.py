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
    parser.add_argument('--split-bnd', action='store_true', default=False)
    parser.add_argument('--no-samples', dest='include_samples',
                        action='store_false', default=True)
    parser.add_argument('--include-strands', action='store_true', default=False)

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = pysam.VariantFile(args.vcf)
    bt = svu.vcf2bedtool(vcf,
                         split_bnd=args.split_bnd,
                         include_samples=args.include_samples,
                         include_strands=args.include_strands)

    bt.saveas(args.bed)
