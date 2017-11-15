#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Standardize a VCF of SV calls.

Each record corresponds to a single SV breakpoint and will have the following
INFO fields, with specified constraints:
  SVTYPE:  SV type [DEL,DUP,INV,BND]
  CHR2:    Secondary chromosome [Must be lexicographically greater than CHROM]
  END:     SV end position (or position on CHR2 in translocations)
  STRANDS: Breakpoint strandedness [++,+-,-+,--]
  SVLEN:   SV length (-1 if translocation)
  SOURCE:  Source algorithm
"""

import argparse
import sys
import pkg_resources
from svtk.standardize import VCFStandardizer
from pysam import VariantFile


def any_called(record):
    null_GTs = [(0, 0), (None, None), (0, ), (None, )]

    def _is_called(sample):
        return record.samples[sample]['GT'] not in null_GTs

    return any([_is_called(sample) for sample in record.samples])


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk standardize',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Raw VCF.')
    parser.add_argument('fout', help='Standardized VCF.')
    parser.add_argument('source', help='Source algorithm. '
                        '[delly,lumpy,manta,wham,melt]')
    parser.add_argument('-p', '--prefix', help='If provided, variant names '
                        'will be overwritten with this prefix.')
    parser.add_argument('--include-reference-sites', action='store_true',
                        default=False, help='Include records where all '
                        'samples are called 0/0 or ./.')
    parser.add_argument('--standardizer', help='Path to python file with '
                        'custom standardizer definition. (Not yet supported.)')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    template = pkg_resources.resource_filename('svtk',
                                               'data/standard_template.vcf')
    template = VariantFile(template)
    vcf = VariantFile(args.vcf)

    # Template header includes all necessary FILTER, INFO, and FORMAT fields
    # Just need to add samples from VCF being standardized
    header = template.header
    for sample in vcf.header.samples:
        header.add_sample(sample)

    # Tag source in header
    meta = '##FORMAT=<ID={0},Number=1,Type=Integer,Description="Called by {1}"'
    meta = meta.format(args.source, args.source.capitalize())
    header.add_line(meta)
    header.add_line('##source={0}'.format(args.source))

    fout = VariantFile(args.fout, mode='w', header=header)

    standardizer = VCFStandardizer.create(args.source, vcf, fout)
    idx = 1
    for record in standardizer.standardize_vcf():
        if any_called(record) or args.include_reference_sites:
            if args.prefix is not None:
                record.id = '{0}_{1}'.format(args.prefix, idx)
                idx += 1

            fout.write(record)

    #  for std_rec in standardize_vcf(vcf, fout):
        #  fout.write(std_rec)

    fout.close()
    vcf.close()


if __name__ == '__main__':
    main(sys.argv[1:])
