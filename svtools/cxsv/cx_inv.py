#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
from collections import namedtuple, defaultdict
import numpy as np
from natsort import natsorted
import svtools.utils as svu


fields = 'chr start end name strands'.split()
Inv = namedtuple('Inv', fields)

fields = 'chr start end name cnv_type'.split()
CNV = namedtuple('CNV', fields)


def breakpoints_match(FF, RR, svtype, mh_buffer=50):
    del_order = ((RR.start > FF.start - mh_buffer) and
                 (FF.end > RR.start) and
                 (RR.end > FF.end - mh_buffer))

    dup5_order = ((RR.start < FF.start) and
                  (FF.start < FF.end) and
                  (FF.end < RR.end + mh_buffer))

    dup3_order = ((FF.start < RR.start + mh_buffer) and
                  (RR.start < RR.end) and
                  (RR.end < FF.end))

    dupINVdup_order = (RR.start < FF.start < RR.end < FF.end)

    if svtype in 'delINV INVdel delINVdel'.split():
        return del_order
    elif svtype in 'dupINV dupINVdel'.split():
        return dup5_order
    elif svtype in 'INVdup delINVdup'.split():
        return dup3_order
    else:
        return dupINVdup_order


def classify_2_cnv(FF, RR, cnv5, cnv3, min_frac=0.5):
    cnv_type5 = cnv5.cnv_type.lower()
    cnv_type3 = cnv3.cnv_type.lower()

    if cnv_type5 == 'del':
        interval5 = (FF.start, RR.start)
    else:
        interval5 = (RR.start, FF.start)
    frac5 = recip_overlap(cnv5.start, cnv5.end, *interval5)

    if cnv_type3 == 'del':
        interval3 = (FF.end, RR.end)
    else:
        interval3 = (RR.end, FF.end)
    frac3 = recip_overlap(cnv3.start, cnv3.end, *interval3)

    if frac5 >= min_frac and frac3 >= min_frac:
        svtype = cnv_type5 + 'INV' + cnv_type3
    elif frac5 >= min_frac and frac3 < min_frac:
        svtype = classify_1_cnv(FF, RR, cnv5)
    elif frac5 < min_frac and frac3 >= min_frac:
        svtype = classify_1_cnv(FF, RR, cnv3)
    else:
        return 'CNV_2_FAIL'

    return svtype


def classify_1_cnv(FF, RR, cnv, min_frac=0.5,
                   min_bkpt_cnv_size=500, max_bkpt_cnv_size=4000):

    cnv_type = cnv.cnv_type.lower()
    if cnv_type == 'del':
        interval5 = (FF.start, RR.start)
        interval3 = (FF.end, RR.end)
    else:
        interval5 = (RR.start, FF.start)
        interval3 = (RR.end, FF.end)

    # First check if paired CNV were likely merged
    start = min(FF.start, RR.start)
    end = max(FF.end, RR.end)
    total_frac = recip_overlap(cnv.start, cnv.end, start, end)
    frac5 = overlap_frac(*interval5, cnv.start, cnv.end)
    frac3 = overlap_frac(*interval3, cnv.start, cnv.end)

    if total_frac > 0.9 and frac5 > 0.95 and frac3 > 0.95:
        svtype = cnv_type + 'INV' + cnv_type  # + '_merged'
        return svtype

    frac5 = recip_overlap(cnv.start, cnv.end, *interval5)
    frac3 = recip_overlap(cnv.start, cnv.end, *interval3)

    if frac5 >= min_frac and frac3 < min_frac:
        svtype = cnv_type + 'INV'

        dist3 = RR.end - FF.end
        if min_bkpt_cnv_size <= dist3 < max_bkpt_cnv_size:
            svtype = svtype + 'del'
        elif min_bkpt_cnv_size <= -dist3 < max_bkpt_cnv_size:
            svtype = svtype + 'dup'

    elif frac5 < min_frac and frac3 >= min_frac:
        svtype = 'INV' + cnv_type

        dist5 = RR.start - FF.start
        if min_bkpt_cnv_size <= dist5 < max_bkpt_cnv_size:
            svtype = 'del' + svtype
        elif min_bkpt_cnv_size <= -dist5 < max_bkpt_cnv_size:
            svtype = 'dup' + svtype

    else:
        return 'CNV_1_unclassified'

    return svtype


def filter_multiple_cnvs(FF, RR, cnvs, min_frac=0.5):
    """
    For cases with 3 or more overlapping CNV, try to remove spurious hits
    by forcing 50% reciprocal with any of the possible CNV intervals. If
    multiple CNVs are present for a candidate interval (e.g. 5' deletion),
    their coordinates are merged by taking the median.

    Parameters
    ----------
    FF : pysam.VariantRecord
        FF inversion breakpoint
    RR : pysam.VariantRecord
        RR inversion breakpoint
    cnvs : list of pysam.VariantRecord
        List of CNVs overlapping breakpoints

    Returns
    -------
    cnvs : list of pysam.VariantRecord
        Filtered and merged CNVs
    """

    # Identify eligible intervals for flanking CNV, defined by inv breakpoints
    del5 = (FF.pos, RR.pos)
    del3 = (FF.info['END'], RR.info['END'])
    dup5 = (RR.pos, FF.pos)
    dup3 = (RR.info['END'], FF.info['END'])

    # Determine if CNV supports 5' CNV, 3' CNV, spans event, or fails overlap
    def _test_overlap(cnv):
        svtype = cnv.info['SVTYPE']
        if svtype == 'DEL':
            frac5 = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], *del5)
            frac3 = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], *del3)
        else:
            frac5 = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], *dup5)
            frac3 = svu.reciprocal_overlap(cnv.pos, cnv.info['END'], *dup3)

        if frac5 >= min_frac and frac3 >= min_frac:
            return svtype + '_53'
        elif frac5 >= min_frac:
            return svtype + '_5'
        elif frac3 >= min_frac:
            return svtype + '_3'
        else:
            return 'no_hit'

    # Collect CNV of same overlap type (e.g., 5' deletion) for merging
    cnvlists = defaultdict(list)
    for cnv in cnvs:
        cnvtype = _test_overlap(cnv)
        if cnvtype == 'no_hit':
            continue
        cnvlists[cnvtype].append(cnv)

    # Keep original CNV if only one present,
    # else merge by taking median start/end
    cnvs = []
    for overlap in cnvlists.keys():
        if len(cnvlists[overlap]) == 1:
            cnvs.append(cnvlists[overlap][0])
        else:
            cnvlist = cnvlists[overlap]
            # Overwrite values in first VariantRecord
            # (can't add list of IDs yet)
            merged_cnv = cnvlist[0]

            # get coordinates
            start = int(np.median([c.pos for c in cnvlist]))
            end = int(np.median([c.info['END'] for c in cnvlist]))
            name = '__'.join([c.name for c in cnvlist])

            merged_cnv.pos = start
            merged_cnv.info['END'] = end
            merged_cnv.id = name

            cnvs.append(merged_cnv)

    return sorted(cnvs, key=lambda record: record.pos)


# TODO:
# Define complex inversion class to store classification, breakpoints, and CNVs

def classify_complex_inversion(FF, RR, cnvs):
    """
    Parameters
    ----------
    FF : pysam.VariantRecord
        FF inversion breakpoint
    RR : pysam.VariantRecord
        RR inversion breakpoint
    cnvs : list of pysam.VariantRecord
        List of overlapping CNVs
    """

    if len(cnvs) > 2:
        cnvs = filter_multiple_cnvs(FF, RR, cnvs)

    if len(cnvs) == 0:
        svtype = 'SIMPLE'
    elif len(cnvs) == 1:
        svtype = classify_1_cnv(FF, RR, cnvs[0])
    elif len(cnvs) == 2:
        svtype = classify_2_cnv(FF, RR, cnvs)
    else:
        svtype = 'MULT_CNVS'

    if breakpoints_match(self.FF, self.RR, svtype, mh_buffer=50):
        return svtype
    else:
        return 'COMPLEX_INS'

#  class DoubleEndLink(Link):
    #  def __init__(self, link):
        #  self.samples = link.samples
        #  self.cnvs = sorted(link.cnvs)
        #  self.invs = sorted(link.invs)
        #  self._svtype = None

        #  if self.invs[0].strands == 'FF':
            #  self.FF = self.invs[0]
            #  self.RR = self.invs[1]
        #  else:
            #  self.FF = self.invs[1]
            #  self.RR = self.invs[0]

        #  self.start = min(self.FF.start, self.RR.start)
        #  self.chrom = self.FF.chr

    #  def get_svtype(self, min_frac=0.5, min_span_frac=0.9,
                   #  min_bkpt_cnv_size=500, max_bkpt_cnv_size=4000):
        #  """
        #  In cases where there is one supporting depth CNV, check if
        #  breakpoints at other end support small CNV

        #  Parameters
        #  ----------
        #  min_bkpt_cnv_size : int
            #  Minimum distance between breakpoints to report CNV w/o depth
        #  max_bkpt_cnv_size : int
            #  Maximum distance between breakpoints to report CNV w/o depth
        #  """

        #  #  check_name = 'DLM.519_20_9087'
        #  #  if self.FF.name == check_name:
            #  #  import ipdb
            #  #  ipdb.set_trace()



#  def parse_complex_sv(inv_poly_overlap):
    #  for link in parse_flanking_cnv(inv_poly_overlap):
        #  #  check_name = 'DLM.519_11_11928'
        #  #  if sorted(link.invs)[0].name == check_name:
            #  #  import ipdb
            #  #  ipdb.set_trace()
        #  if link.svtype == 'CANDIDATE':
            #  link = DoubleEndLink(link)
            #  svtype = link.svtype
            #  exclude = 'CNV_1_unclassified CNV_2_FAIL COMPLEX_MULTI_CNV'.split()
            #  if svtype not in exclude:
                #  yield link

