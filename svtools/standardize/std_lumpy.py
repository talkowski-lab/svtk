"""
std_lumpy.py

Standardize Lumpy records.

# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.
"""


from .standardize import VCFStandardizer, parse_bnd_pos


@VCFStandardizer.register('lumpy')
class LumpyStandardizer(VCFStandardizer):
    def standardize_vcf(self):
        for record in self.raw_vcf:
            # Split inversion events into their constituent breakpoints
            # For each strandedness in record, make a corresponding std record
            strands = record.info['STRANDS']
            for strand in strands:
                record.info['STRANDS'] = (strand, )
                std_rec = self.std_vcf.new_record()
                yield self.standardize_record(std_rec, record)

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
            chr2, end = raw_rec.chrom, raw_rec.info['END']

        std_rec.info['CHR2'] = chr2
        std_rec.info['END'] = end

        # Add SVLEN
        if std_rec.chrom == std_rec.info['CHR2']:
            std_rec.info['SVLEN'] = end - std_rec.pos
        else:
            std_rec.info['SVLEN'] = -1

        # Strip per-strand counts
        std_rec.info['STRANDS'] = raw_rec.info['STRANDS'][0].split(':')[0]

        std_rec.info['SOURCE'] = 'lumpy'

        return std_rec
