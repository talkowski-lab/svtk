# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.
"""
std_lumpy.py

Standardize Lumpy records.
"""


from .standardize import VCFStandardizer, parse_bnd_pos


@VCFStandardizer.register('lumpy')
class LumpyStandardizer(VCFStandardizer):
    def standardize_vcf(self):
        for record in self.raw_vcf:
            # Split inversion events into their constituent breakpoints
            # For each strandedness in record, make a corresponding std record
            strands = record.info['STRANDS']
            for i, strand in enumerate(strands):
                record.info['STRANDS'] = (strand, )
                std_rec = self.std_vcf.new_record()
                std_rec = self.standardize_record(std_rec, record)

                # Some variants have stranded pairs that don't match their
                # SV type
                if std_rec.info['SVTYPE'] == 'DEL':
                    if std_rec.info['STRANDS'] != '+-':
                        continue
                if std_rec.info['SVTYPE'] == 'DUP':
                    if std_rec.info['STRANDS'] != '-+':
                        continue
                if std_rec.info['SVTYPE'] == 'INV':
                    if std_rec.info['STRANDS'] not in '++ --'.split():
                        continue

                # Tag split record
                std_rec.id += 'abcd'[i]

                yield std_rec

    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Lumpy record.

        1) Add CHR2, END
        """

        std_rec.info['SVTYPE'] = raw_rec.info['SVTYPE']

        # Parse CHR2 and END
        if std_rec.info['SVTYPE'] == 'BND':
            chr2, end = parse_bnd_pos(std_rec.alts[0])
        else:
            chr2, end = raw_rec.chrom, raw_rec.stop

        std_rec.info['CHR2'] = chr2
        std_rec.stop = end

        # Add SVLEN
        if std_rec.chrom == std_rec.info['CHR2']:
            std_rec.info['SVLEN'] = end - std_rec.pos
        else:
            std_rec.info['SVLEN'] = -1

        # Strip per-strand counts
        std_rec.info['STRANDS'] = raw_rec.info['STRANDS'][0].split(':')[0]

        std_rec.info['SOURCES'] = ['lumpy']

        return std_rec

    def standardize_format(self, std_rec, raw_rec):
        """
        Parse called samples from TAGS field
        """

        source = std_rec.info['SOURCES'][0]

        # Any sample in TAGS field is considered to be called
        for sample in raw_rec.samples:
            if raw_rec.samples[sample]['SU'] >= 4:
                std_rec.samples[sample]['GT'] = (0, 1)
                std_rec.samples[sample][source] = 1
            else:
                std_rec.samples[sample]['GT'] = (0, 0)
                std_rec.samples[sample][source] = 0

        return std_rec
