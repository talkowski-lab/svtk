#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Resolve complex SV from inversion/translocation breakpoints and CNV intervals.
"""

import argparse
import sys
from svtools.cxsv import link_cpx, resolve_cpx


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtools link-cpx',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Breakpoint VCFs.')
    parser.add_argument('bed', type=argparse.FileType('w'),
                        help='Resolved complex variants.')
    parser.add_argument('unresolved', type=argparse.FileType('w'),
                        help='Unresolved complex breakpoints.')
    parser.add_argument('-p', '--prefix', default='CPX_',
                        help='Variant prefix [CPX_]')

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    resolved_classes = ('INV delINV INVdel delINVdel dupINV INVdup dupINVdup '
                        'delINVdup dupINVdel DUP5/INS3 DUP3/INS5 COMPLEX_INS')
    resolved_classes = resolved_classes.split()

    clusters = link_cpx(args.vcf)
    for i, cluster in enumerate(clusters):
        cpx_type, entry = resolve_cpx(cluster)
        entry = entry.format(name=args.prefix + str(i + 1))

        if cpx_type in resolved_classes:
            args.bed.write(entry)
        else:
            args.unresolved.write(entry)


if __name__ == '__main__':
    main(sys.argv[1:])
