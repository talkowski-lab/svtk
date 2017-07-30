# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import pandas as pd
import pysam
from .gencode_elements import annotate_gencode_elements
from .classify_effect import classify_effect
from .nearest_tss import annotate_nearest_tss
import svtools.utils as svu


def annotate_gencode(sv, gencode):
    """
    Gene disrupting if:
       1) Breakpoint hits exon
       2) Inverted segment lies inside gene AND hits exon
       3) One end of inverted segment lies inside gene AND other does not

    Notes:
        * Intragenic variants are annotated with the nearest TSS.
        * Inversions which entirely span a gene are considered "no_effect" and
        excluded from annotation
        * Duplications which overlap the start OR end of a gene boundary
        without fully spanning the boundary are considered "no_effect" and
        excluded from annotation

    Parameters
    ----------
    sv : pbt.BedTool
        SV breakpoints and CNV intervals
    gencode : pbt.BedTool
        Gencode annotations
    """

    # Annotate Gencode hits and predicted effects
    hits = annotate_gencode_elements(sv, gencode)
    effects = classify_effect(hits)
    effects = effects.loc[effects.effect != 'GENE_OTHER'].copy()

    # Annotate nearest TSS
    tss = annotate_nearest_tss(sv, gencode)

    # Only include TSS if no genic hit observed
    tss = tss.loc[~tss.name.isin(effects.name)].copy()

    # Merge annotations
    effects = pd.concat([effects, tss])

    # Replace ENSEMBL gene IDs with gene names
    gene_key = {}
    with open(gencode.fn) as gencode_f:
        for line in gencode_f:
            data = line.strip().split('\t')
            if data[7] != 'gene':
                continue
            fields = data[9].strip(';').split('; ')
            gene_id = fields[0].split()[1].strip('"')
            gene_name = fields[4].split()[1].strip('"')
            gene_key[gene_id] = gene_name
    effects['gene_name'] = effects['gene_id'].replace(gene_key)

    # Aggregate genic effects by variant ID
    effects = effects.pivot_table(index='name',
                                  columns='effect',
                                  values='gene_name',
                                  aggfunc=lambda s: ','.join(sorted(set(s))),
                                  fill_value='NA')

    return effects


GENCODE_INFO = [
    '##INFO=<ID=LOF,Number=.,Type=String,Description="Gene(s) on which the SV is predicted to have a loss-of-function effect.">',
    '##INFO=<ID=COPY_GAIN,Number=.,Type=String,Description="Gene(s) on which the SV is predicted to have a copy-gain effect.">',
    '##INFO=<ID=INTRONIC,Number=.,Type=String,Description="Gene(s) where the SV was found to lie entirely within an intron.">',
    '##INFO=<ID=DUP_PARTIAL,Number=.,Type=String,Description="Gene(s) which are partially overlapped by an SV\'s duplication, such that an unaltered copy is preserved.">',
    '##INFO=<ID=INV_SPAN,Number=.,Type=String,Description="Gene(s) which are entirely spanned by an SV\'s inversion.">',
    '##INFO=<ID=NEAREST_TSS,Number=.,Type=String,Description="Nearest transcription start site to intragenic variants.">',
    '##INFO=<ID=INTRAGENIC,Number=0,Type=Flag,Description="SV does not overlap coding sequence.">'
]


def annotate_vcf(vcf, gencode, annotated_vcf):
    """
    Parameters
    ----------
    vcf : pysam.VariantFile
    gencode : pbt.BedTool
        Gencode annotation bed
    annotated_vcf : str
        Path to output VCF
    """

    # Add metadata lines for annotations
    header = vcf.header
    for line in GENCODE_INFO:
        header.add_line(line)

    # Open output file
    fout = pysam.VariantFile(annotated_vcf, 'w', header=header)

    # Annotate genic hits
    sv = svu.vcf2bedtool(vcf.filename)

    effects = annotate_gencode(sv, gencode)
    effects = effects.to_dict(orient='index')

    # Add results to variant records and save
    for record in vcf:
        anno = effects[record.id]
        for info, genelist in anno.items():
            if genelist != 'NA':
                record.info[info] = genelist

        if 'NEAREST_TSS' in record.info:
            record.info['INTRAGENIC'] = True

        fout.write(record)

    fout.close()
