#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Classification of reciprocal translocations.
"""


def classify_complex_translocation(tloc1, tloc2, mh_buffer=100):
    """
    Classify the complex class of an inversion and associated CNVs.

    Parameters
    ----------
    FF : pysam.VariantRecord
        FF inversion breakpoint.
    RR : pysam.VariantRecord
        RR inversion breakpoint.
    cnvs : list of pysam.VariantRecord
        List of overlapping CNVs.

    Returns
    -------
    svtype : str
        Complex SV class.
    """

    # Derivative chromosomes are labeled 1 and 2
    # Reference chromosomes are labeled A and B
    # A1 = the breakend of ref chrom A on derivative chrom 1

    if tloc1.chrom != tloc2.chrom or tloc1.info['CHR2'] != tloc2.info['CHR2']:
        return 'TLOC_MISMATCH_CHROM'

    # get positions
    A1 = tloc1.pos
    A2 = tloc2.pos
    B1 = tloc1.info['END']
    B2 = tloc2.info['END']
    strand1 = tloc1.info['STRANDS']
    strand2 = tloc2.info['STRANDS']

    if strand1 == '++' and strand2 == '--':
        if (A2 > A1 - mh_buffer) and (B2 > B1 - mh_buffer):
            return 'CTX_PQ/QP'
    elif strand1 == '+-' and strand2 == '-+':
        if (A2 > A1 - mh_buffer) and (B1 > B2 - mh_buffer):
            return 'CTX_PP/QQ'
    elif strand1 == '-+' and strand2 == '+-':
        if (A1 > A2 - mh_buffer) and (B2 > B1 - mh_buffer):
            return 'CTX_PP/QQ'
    elif strand1 == '--' and strand2 == '++':
        if (A1 > A2 - mh_buffer) and (B1 > B2 - mh_buffer):
            return 'CTX_PQ/QP'

    return 'CTX_INS'
