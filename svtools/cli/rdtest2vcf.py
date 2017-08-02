#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Convert an RdTest-formatted bed to the standard VCF format.

Each record corresponds to a single CNV interval and will have the following
INFO fields, with specified constraints:
  SVTYPE:  SV type [DEL,DUP]
  CHR2:    Secondary chromosome; set to same as CHROM
  END:     SV end position
  STRANDS: Breakpoint strandedness [DEL:+-,DUP:-+]
  SVLEN:   SV length
  SOURCE:  Tagged with "depth"
"""

import argparse
import sys
import pkg_resources
from collections import namedtuple
from pysam import VariantFile


def RdtestParser(bed):
    CNV = namedtuple('CNV', 'chrom start end name samples svtype'.split())

    for line in bed:
        if line.startswith('#'):
            continue
        data = line.strip().split()

        chrom = data[0]
        start = int(data[1])
        end = int(data[2])
        name = data[3]
        samples = data[4].split(',')
        svtype = data[5].upper()

        yield CNV(chrom, start, end, name, samples, svtype)


def rdtest2vcf(bed, vcf):
    for cnv in RdtestParser(bed):
        record = vcf.new_record()
        record.chrom = cnv.chrom

        if cnv.start < 1:
            record.pos = 1
        else:
            record.pos = cnv.start
        record.id = cnv.name

        record.ref = 'N'
        record.alts = ('<{0}>'.format(cnv.svtype), )

        record.filter.add('PASS')

        # Add required INFO fields
        record.info['SVTYPE'] = cnv.svtype
        record.info['CHR2'] = cnv.chrom
        record.stop = cnv.end
        record.info['SVLEN'] = cnv.end - cnv.start
        record.info['SOURCES'] = ['depth']
        if cnv.svtype == 'DEL':
            record.info['STRANDS'] = '+-'
        elif cnv.svtype == 'DUP':
            record.info['STRANDS'] = '-+'

        # Seed with ref genotypes
        for sample in vcf.header.samples:
            record.samples[sample]['GT'] = (0, 0)
            record.samples[sample]['depth'] = 0

        # Call any samples with variant as heterozygous
        called = 0
        for sample in cnv.samples:
            if sample in vcf.header.samples:
                called += 1
                record.samples[sample]['GT'] = (0, 1)
                record.samples[sample]['depth'] = 1

        if called > 0:
            vcf.write(record)


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtools standardize',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bed', type=argparse.FileType('r'),
                        help='RdTest-formatted bed file. '
                        '(chrom, start, end, name, samples, svtype)')
    parser.add_argument('samples', help='List of all samples present in '
                        'variant callset.')
    parser.add_argument('fout', help='Standardized VCF.')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Get template header
    template = pkg_resources.resource_filename('svtools',
                                               'data/standard_template.vcf')
    template = VariantFile(template)
    header = template.header

    # Get list of samples
    with open(args.samples) as slist:
        samples = sorted([s.strip() for s in slist.readlines()])

    # Template header includes all necessary FILTER, INFO, and FORMAT fields
    # Just need to add list of samples
    for sample in samples:
        header.add_sample(sample)

    # Tag source in header
    meta = ('##FORMAT=<ID=depth,Number=1,Type=Integer,'
            'Description="Called by read-depth algorithms">')
    header.add_line(meta)
    header.add_line('##source=depth')

    fout = VariantFile(args.fout, mode='w', header=header)

    rdtest2vcf(args.bed, fout)

if __name__ == '__main__':
    main(sys.argv[1:])
