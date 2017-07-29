#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Annotate a VCF of structural variants with a list of genomic elements.
"""

import argparse
import pybedtools as pbt
import svtools.utils as svu


def intersection_type(variant, element):
    """
    list of str
        chrom, start, end, ID
    """
    variant_start, variant_end = [int(x) for x in variant[1:3]]
    element_start, element_end = [int(x) for x in element[1:3]]

    if variant_start > element_start and variant_end < element_end:
        return 'BOTH-INSIDE'
    elif variant_start > element_start and variant_start < element_end:
        return 'ONE-INSIDE'
    elif variant_end > element_start and variant_end < element_end:
        return 'ONE-INSIDE'
    else:
        return 'SPAN'


def disruption_type(hit_type, svtype):
    disruptions = {
        'DEL': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'DISRUPTING'},
        'DUP': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'COPY_GAIN'},
        'INV': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'SPAN'},
        'BND': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'DISRUPTING'},
        }

    return disruptions[svtype][hit_type]


def annotate_gencode_elements(sv, gencode):
    """
    Gene disrupting if:
       1) Breakpoint hits exon
       2) Inverted segment lies inside gene AND hits exon
       3) One end of inverted segment lies inside gene AND other does not

    Parameters
    ----------
    sv : pbt.BedTool
        SV breakpoints and CNV intervals
    gencode : pbt.BedTool
        Gencode annotations
    """

    # Number of fields in SV bedtool
    N_BED_FIELDS = 6

    # Check intersection with gene boundaries
    sect = sv.intersect(gencode, wa=True, wb=True)

    for hit in sect.intervals:
        variant = hit.fields[:N_BED_FIELDS]
        variant_ID = variant[3]
        svtype = variant[4]

        # Gencode data
        element = hit.fields[N_BED_FIELDS:]
        gene_ID = element[3]
        element_type = element[7]

        hit_type = intersection_type(variant, element)
        disrupt_type = disruption_type(hit_type, svtype)

        yield (variant_ID, svtype, gene_ID, element_type, hit_type,
               disrupt_type)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Structural variants.')
    parser.add_argument('gencode_annotation', help='Gencode annotation bed.')
    parser.add_argument('fout', type=argparse.FileType('w'))
    args = parser.parse_args()

    sv = svu.vcf2bedtool(args.vcf)
    gencode = pbt.BedTool(args.gencode_annotation)

    annotations = annotate_gencode_elements(sv, gencode)

    header = 'name svtype gene_id element_type hit_type disrupt_type'.split()
    args.fout.write('\t'.join(header) + '\n')

    for hit in annotations:
        entry = '\t'.join(hit) + '\n'
        args.fout.write(entry)


if __name__ == '__main__':
    main()
