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
import pybedtools as pbt
import svtools.annotation as anno


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('vcf', help='Structural variants.')
    parser.add_argument('gencode_annotation', help='Gencode annotation bed.')
    parser.add_argument('annotated_vcf', help='Annotated variants.')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = pysam.VariantFile(args.vcf)
    gencode = pbt.BedTool(args.gencode_annotation)

    anno.annotate_vcf(vcf, gencode, args.annotated_vcf)


if __name__ == '__main__':
    main(sys.argv[1:])
