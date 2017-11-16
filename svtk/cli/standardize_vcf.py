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
    parser.add_argument('--contigs', type=argparse.FileType('r'),
                        help='Reference fasta index (.fai). If provided, '
                        'contigs in index will be used in VCF header. '
                        'Otherwise all GRCh37 contigs will be used in header. '
                        'Variants on contigs not in provided list will be '
                        'removed.')
    parser.add_argument('--min-size', type=int, default=50,
                        help='Minimum SV size to report [50].')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Add contigs to header if provided
    if args.contigs:
        template = pkg_resources.resource_filename(
                'svtk', 'data/no_contigs_template.vcf')
        template = VariantFile(template)
        header = template.header
        contig_line = '##contig=<ID={contig},length={length}>'
        for line in args.contigs:
            contig, length = line.split()[:2]
            header.add_line(contig_line.format(**locals()))
    # Use GRCh37 by default
    else:
        template = pkg_resources.resource_filename(
                'svtk', 'data/GRCh37_template.vcf')
        template = VariantFile(template)
        header = template.header

    vcf = VariantFile(args.vcf)

    # Template header includes all necessary FILTER, INFO, and FORMAT fields
    # Just need to add samples from VCF being standardized
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
        # Skip variants on non-whitelisted contigs
        if record.chrom not in fout.header.contigs:
            continue
        if record.info['CHR2'] not in fout.header.contigs:
            continue

        # Apply size filter (but keep breakends (SVLEN=-1))
        if 0 < record.info['SVLEN'] < args.min_size:
            continue

        # Only report sites with called samples unless requested otherwise
        if any_called(record) or args.include_reference_sites:
            if args.prefix is not None:
                record.id = '{0}_{1}'.format(args.prefix, idx)
                idx += 1

            fout.write(record)

    fout.close()
    vcf.close()


if __name__ == '__main__':
    main(sys.argv[1:])
