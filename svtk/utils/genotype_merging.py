#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""


def choose_best_genotype(sample, records):
    """
    Return record where sample has best non-reference genotype, if available

    Parameters
    ----------
    sample : str
    records : list of pysam.VariantRecord

    Returns
    -------
    best_record : pysam.VariantRecord
    """

    best_GT = (0, 0)
    best_GQ = 0 
    best_record = None

    # Pick best non-reference genotype
    for record in records:
        if record.samples[sample]['GQ'] > best_GQ:
            # if record is non-reference , use it
            # or if it's a higher GQ for a reference call, use it
            if record.samples[sample]['GT'] != (0, 0) or best_GT == (0, 0):
                best_GT = record.samples[sample]['GT']
                best_GQ = record.samples[sample]['GQ']
                best_record = record

    return best_record


def update_best_genotypes(new_record, records, is_multiallelic=False):
    """
    For each sample in record, update GT and other formats with best genotype

    Parameters
    ----------
    new_record : pysam.VariantRecord
    records : list of SVRecord
    is_multiallelic : bool
    """
    records = [r.record for r in records]

    for sample in new_record.samples:
        best_record = choose_best_genotype(sample, records)
        for key in best_record.format.keys():
            if key in new_record.header.formats.keys():
                # If any record is multiallelic, replace non-multiallelic
                # genotypes with multiallelic equivalents
                if key == 'GT' and is_multiallelic:
                    GT = best_record.samples[sample][key]
                    if GT == (0, 1):
                        new_record.samples[sample][key] = (0, 2)
                    elif GT == (1, 1):
                        new_record.samples[sample][key] = (2, 2)
                    else:
                        new_record.samples[sample][key] = GT
                else:
                    new_record.samples[sample][key] = best_record.samples[sample][key]
