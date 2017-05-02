# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.
"""
std_wham.py

Standardize WHAM record.
"""


from .standardize import VCFStandardizer


@VCFStandardizer.register('wham')
class WhamStandardizer(VCFStandardizer):
    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Wham record.
        """

        # colons in the ID break PyVCF parsing
        std_rec.id = '_'.join(std_rec.id.split(':'))

        # commas conflict with VCF format field with num='.'
        std_rec.id = ';'.join(std_rec.id.split(','))

        # SV types already standard
        std_rec.info['SVTYPE'] = raw_rec.info['SVTYPE']

        # No CTX events
        std_rec.info['CHR2'] = std_rec.chrom
        std_rec.info['END'] = raw_rec.info['END']

        # Strand not provided for inv/tloc
        if std_rec.info['SVTYPE'] == 'DEL':
            strands = '+-'
        elif std_rec.info['SVTYPE'] == 'DUP':
            strands = '-+'
        else:
            strands = '.'
        std_rec.info['STRANDS'] = strands

        # SVLEN is a list and can be negative
        std_rec.info['SVLEN'] = abs(raw_rec.info['SVLEN'][0])

        std_rec.info['SOURCES'] = ['wham']

        return std_rec
