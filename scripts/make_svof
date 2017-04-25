#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Convert a clustered bed file with metrics to an SVOF
"""

import argparse
from collections import namedtuple, defaultdict, Counter
import numpy as np

bed_fields = 'chrom start end name sample source svtype cname cfreq'.split()
bed_fields = bed_fields + ['INFO']
BedRecord = namedtuple('BedRecord', bed_fields)


class SVCall(object):
    """Stand in for old API"""
    def __init__(self, chrA, posA, chrB, posB, name, sample, source, svtype):
        self.chrA = chrA
        self.posA = int(posA)
        self.chrB = chrB
        self.posB = int(posB)
        self.name = name
        self.sample = sample
        self.source = source
        self.svtype = svtype
        self.cname = None
        self.INFO = None


def parse_clusters(bed, metriclist, sep='\t'):
    prev = None
    calls = []

    for line in bed:
        if line.startswith('#'):
            continue

        data = line.strip().split(sep)
        metrics = {k: v for k, v in zip(metriclist, data[9:])}
        rec = BedRecord(*data[:9], metrics)

        call = SVCall(rec.chrom, rec.start, rec.chrom, rec.end,
                      rec.name, rec.sample, rec.source, rec.svtype)
        call.cname = rec.cname
        call.INFO = rec.INFO

        if prev is None or prev.cname == rec.cname:
            calls.append(call)
            prev = rec
            continue
        else:
            yield calls
            calls = [call]
            prev = rec

    yield calls


merged_fields = 'chrA posA posB chrB name svtype sample sources'.split()
merged_fields = merged_fields + ['INFO']
MergedCall = namedtuple('MergedCall', merged_fields)


def merge(calls, svtype=None):
    chrA = calls[0].chrA
    chrB = calls[0].chrB
    try:
        posA = int(np.median([call.posA for call in calls]))
    except:
        import pdb
        pdb.set_trace()
    posB = int(np.median([call.posB for call in calls]))
    name = calls[0].cname

    sample = calls[0].sample
    sources = sorted(call.source for call in calls)

    INFO = defaultdict(list)
    for call in calls:
        for k, v in call.INFO.items():
            if v == '':
                continue
            try:
                INFO[k].append(float(v))
            except ValueError:
                INFO[k].append(v)

    med_INFO = defaultdict(str)
    for k, vlist in INFO.items():
        if len(vlist) == 0:
            med_INFO[k] = ''
        else:
            try:
                # TODO: adapt SVRecord merge
                #  med_INFO[k] = str(np.median(vlist))
                med_INFO[k] = str(np.sum(vlist))
            # String types
            # FT
            except TypeError:
                try:
                    med_INFO[k] = ','.join(vlist)
                    msg = 'String values found in numeric INFO: {name}, {fmt}'
                    msg = msg.format(name=name, fmt=k)
                    msg = msg + '\n' + med_INFO[k]
                    #  print(msg)
                # GSDEPTHRATIO
                except TypeError:
                    med_INFO[k] = ','.join([str(v) for v in vlist])
                    msg = 'String values found in numeric INFO: {name}, {fmt}'
                    msg = msg + '. Float values also found'
                    msg = msg.format(name=name, fmt=k)
                    msg = msg + '\n' + med_INFO[k]
                    #  print(msg)

    if svtype is None:
        svtype = calls[0].svtype

    m = MergedCall(chrA, posA, posB, chrB, name, svtype, sample, sources,
                   med_INFO)
    return m


def make_svof(calls, sources, metriclist, split_sv=False):
    calldict = defaultdict(list)
    entries = []

    for call in calls:
        calldict[call.sample].append(call)
    for sample in calldict:
        if split_sv:
            svdict = defaultdict(list)
            for call in calldict[sample]:
                svdict[call.svtype].append(call)
            merged = [merge(svdict[svtype], svtype) for svtype in svdict]
        else:
            merged = [merge(calldict[sample], 'cnv')]

        for call in merged:
            # check for overlapping calls from same caller
            base = [str(s) for s in call[:7]]

            obs = Counter(call.sources)
            obs = [str(obs[s]) for s in sources]

            metrics = [call.INFO[k] for k in metriclist]
            #  metrics = ['%.03f' % m for m in metrics]

            entries.append('\t'.join(base + obs + metrics))

    return '\n'.join(entries)


def main():
    parser = argparse.ArgumentParser(
        description="")
    parser.add_argument('clusterbed', type=argparse.FileType('r'))
    parser.add_argument('sources', help='comma delimited')
    parser.add_argument('svof', type=argparse.FileType('w'))
    parser.add_argument('-s', '--split', action='store_true', default=False)
    args = parser.parse_args()

    sources = args.sources.split(',')

    header = 'chrA posA posB chrB name svtype sample'.split()
    header = header + sources

    h = next(args.clusterbed)
    h = h.strip().split()
    metrics = h[9:]
    header = header + metrics

    args.svof.write('\t'.join(header) + '\n')

    for calls in parse_clusters(args.clusterbed, metrics):
        args.svof.write(make_svof(calls, sources, metrics, args.split) + '\n')


if __name__ == '__main__':
    main()
