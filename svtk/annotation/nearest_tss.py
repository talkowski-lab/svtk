#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import warnings
from .annotate_intersection import split_gencode_fields


def annotate_nearest_tss(sv, gencode):
    """
    Annotate each variant record with its nearest TSS.

    Parameters
    ----------
    sv : pbt.BedTool
        columns = (chrom, start, end, name, svtype, strands)
    gencode : pbt.BedTool
        Gencode gene annotations (GTF)

    Returns
    -------
    nearest_tss : pd.DataFrame
        Columns = (name, svtype, gene_name, effect)
    """

    def _make_tss(feature):
        strand = feature.fields[6]
        if strand == '-':
            feature.start = feature.end

        feature.end = feature.start + 1
        return feature

    transcripts = gencode.filter(lambda r: r.fields[2] == 'transcript')
    tss = transcripts.each(_make_tss).saveas()

    nearest_tss = sv.sort().closest(tss.sort(), D='b', id=True).saveas()

    # Suppress warning about column names
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        nearest_tss = nearest_tss.to_dataframe()

    def _get_gene_name(field_str):
        fields = split_gencode_fields(field_str)
        return fields['gene_name']

    nearest_tss['gene_name'] = nearest_tss[14].apply(_get_gene_name)
    nearest_tss = nearest_tss[[3, 4, 'gene_name']].copy()
    nearest_tss.columns = 'name svtype gene_name'.split()
    nearest_tss['effect'] = 'NEAREST_TSS'

    return nearest_tss
