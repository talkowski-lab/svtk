"""
std_manta.py

Standardize a Manta record.

Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
Distributed under terms of the MIT license.
"""

# Note about parsing strands from BND ALT:
# ALT Meaning
# t[p[ piece extending to the right of p is joined after t (+-)
# t]p] reverse comp piece extending left of p is joined after t (++)
# ]p]t piece extending to the left of p is joined before t (-+)
# [p[t reverse comp piece extending right of p is joined before t (--)


def standardize_manta(raw_rec, std_rec):
    """
    Standardize Manta record.

    1) Replace colons in ID with underscores (otherwise can break VCF parsing)
    2) Define CHR2 and END
    3) Add strandedness
    4) Add SVLEN
    """

    # Colons in the ID can break parsing
    std_rec.id = '_'.join(std_rec.id.split(':'))

    svtype = raw_rec.info['SVTYPE']
    std_rec.info['SVTYPE'] = svtype

    # Define CHR2 and END
    if svtype == 'BND':
        alt = raw_rec.alts[0].strip('ATCGN')
        # Strip brackets separately, otherwise GL contigs will be altered
        alt = alt.strip('[]')
        chr2, end = alt.split(':')
        end = int(end)
    else:
        chr2 = raw_rec.chrom
        end = raw_rec.info['END']
    std_rec.info['CHR2'] = chr2
    std_rec.info['END'] = end

    # Strand parsing
    if svtype == 'INV':
        if 'INV3' in raw_rec.info.keys():
            strands = '++'
        else:
            strands = '--'
    elif svtype == 'BND':
        alt = raw_rec.alts[0]
        if alt.endswith('['):
            strands = '+-'
        elif alt.endswith(']'):
            strands = '++'
        elif alt.startswith(']'):
            strands = '-+'
        elif alt.startswith('['):
            strands = '--'
    elif svtype == 'DEL':
        strands = '+-'
    elif svtype == 'DUP':
        strands = '-+'
    elif svtype == 'INS':
        strands = '.'
    std_rec.info['STRANDS'] = strands

    if svtype == 'BND' and std_rec.chrom != std_rec.info['CHR2']:
        std_rec.info['SVLEN'] = -1
    else:
        std_rec.info['SVLEN'] = std_rec.info['END'] = std_rec.pos

    std_rec.info['SOURCE'] = 'manta'

    return std_rec
