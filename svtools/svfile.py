"""
svfile.py

Wrap the pysam API to permit clustering of standardized SV VCF records.

Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
Distributed under terms of the MIT license.
"""

import os
from collections import OrderedDict, defaultdict, deque
from functools import partial
import numpy as np
import vcf
from vcf.model import _Record, _Call, _Breakend, _SV, make_calldata_tuple

from pysam import VariantFile

from .genomeslink import GSNode
from .constants import NULL_GT, STD_FORMATS, MEDIANS


class UnsupportedFiletypeError(Exception):
    """Unsupported or unexpected filetype"""


class UnboundedFetchError(Exception):
    """Only start position of tabix fetch region specified, no end"""


class TabixNotFoundError(Exception):
    """No Tabix index found"""


class PyVCFParsingError(Exception):
    """Unexpected datatype from PyVCF"""


class MissingInfoError(Exception):
    """Missing required INFO field"""


class SVFile(object):
    def __init__(self, filename):
        """
        Wrapper for standardized VCF files.
        """
        self.filename = filename

        self.reader = VariantFile(filename)

        # Confirm all standard INFO fields are present
        required_info = 'SVTYPE CHR2 END STRANDS SVLEN SOURCE'.split()
        for info in required_info:
            if info not in self.reader.header.info.keys():
                msg = "Required INFO field {0} not found in file {1}"
                msg = msg.format(info, filename)
                raise MissingInfoError(msg)

    def fetch(self, chrom, start=None, end=None):
        """
        Fetch only calls from specified region of SVFile.

        Requires Tabix index. Updates SVFile in place.

        Parameters
        ----------
        chrom : str
        start : int, optional
        end : int, optional
            Required if start specified

        Raises
        ------
        UnboundedFetchError
        TabixNotFoundError
        """
        if start is not None and end is None:
            msg = 'Start {}:{} specified but no end coordinate provided'
            msg = msg.format(chrom, start)
            raise UnboundedFetchError(msg)

        try:
            self.reader = self.reader.fetch(chrom, start, end)
        except ValueError:
            raise TabixNotFoundError(self.filename)

    @property
    def samples(self):
        """
        Returns
        ------
        samples : list of str
        """
        if hasattr(self.reader.header, 'samples'):
            return list(self.reader.header.samples)
        else:
            return []

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        record = next(self.reader)
        return SVRecord(record)


def recip(startA, endA, startB, endB, frac):
    if frac == 0:
        return True

    start = max(startA, startB)
    end = min(endA, endB)
    olen = end - start
    lenA = endA - startA
    lenB = endB - startB

    try:
        lapA = olen / float(lenA)
        lapB = olen / float(lenB)
    except ZeroDivisionError:
        return False

    return (olen > 0) and (lapA >= frac) and (lapB >= frac)


class SVRecord(_Record, GSNode):
    """
    Extend VCF record to add clusterability
    """

    def __init__(self, record):
        """
        record : vcf._Record
            Must specify 'CHR2' and 'END' in INFO
        """

        _Record.__init__(self, record.CHROM, record.POS, record.ID, record.REF,
                         record.ALT, record.QUAL, record.FILTER, record.INFO,
                         record.FORMAT, record._sample_indexes, record.samples)

        GSNode.__init__(self, self.CHROM, self.POS, self.INFO['CHR2'],
                        self.INFO['END'], self.ID)

        self.source = self.INFO['SOURCES'][0]

    def clusters_with(self, other, dist, frac=0.0, match_strands=False):
        if match_strands:
            strand = self.INFO['STRANDS'] == other.INFO['STRANDS']
        else:
            strand = True

        if self.svtype == 'tloc' and other.svtype == 'tloc':
            overlap = True
        else:
            overlap = recip(self.posA, self.posB, other.posA, other.posB, frac)

        return(super().clusters_with(other, dist) and
               self.svtype == other.svtype and overlap and strand)

    @property
    def svtype(self):
        """
        Returns
        -------
        svtype : str
            One of {del, dup, inv tloc}
        """
        SIMPLE_SVTYPES = {
            'DEL': 'del',
            'DUP': 'dup',
            'DUP:TANDEM': 'dup',
            'INV': 'inv',
            'TRA': 'tloc',
            'INS': 'ins',
            'RPL': 'del',
            'BND': 'inv'}

        if self.CHROM != self.INFO['CHR2']:
            return 'tloc'
        else:
            return SIMPLE_SVTYPES[self.INFO['SVTYPE']]

    @staticmethod
    def merge_pos(records):
        """
        Compute aggregate POS/END of clustered SVRecords.

        Defaults to computing median of POS and END over constituent
        records. CIPOS/CIEND are calculated as (MIN - MEDIAN, MAX - MEDIAN)
        over POS and END.

        Parameters
        ----------
        records : list of SVRecord

        Returns
        -------
        POS : int
        END : int
        CIPOS : list of int
        CIEND : list of int
        """

        # Position bounds
        MIN_POS = min(rec.POS for rec in records)
        MAX_POS = max(rec.POS for rec in records)
        MIN_END = min(rec.INFO['END'] for rec in records)
        MAX_END = max(rec.INFO['END'] for rec in records)

        POS = int(np.median([rec.POS for rec in records]))
        END = int(np.median([rec.INFO['END'] for rec in records]))
        CIPOS = [MIN_POS - POS, MAX_POS - POS]
        CIEND = [MIN_END - END, MAX_END - END]

        return POS, END, CIPOS, CIEND

    @classmethod
    def merge(cls, records, samples, sources, ID='.',
              MEDIANS=MEDIANS, summary_fn=np.median):
        """
        Aggregates VCF INFO fields across clustered records.

        * Original record IDs are preserved in a new INFO field
        * Fields common to multiple sources are renamed to
        ${FIELDNAME}_${SOURCENAME}, and master field is assigned to median/sum

        Parameters
        ----------
        others : list of SVRecord
        ID : str
            ID for merged record
        samples : list of str
            List of sample IDs to include in list of calls
        sources : list of str
            List of source algorithms to parse INFO/FORMAT from
        MEDIANS : dict, {str : (list of str)}
            Numeric fields are summed by default. Fields listed here will be
            combined using an alternative summary function. Keys are source
            program names, values are lists of field names
        summary_fn : function
            of summing (defaults to np.median)
            Only merge positions and IDs, not metrics
        """

        # Secondary records have duplicate POS/END/INFO
        records = [rec for rec in records if not ('SECONDARY' in rec.INFO)]

        if len(records) == 0:
            return None

        record = records[0]

        # Constant INFO fields
        CHROM = record.CHROM
        REF = 'N'
        CHR2 = record.INFO['CHR2']
        SVTYPE = record.INFO['SVTYPE']

        POS, END, CIPOS, CIEND = cls.merge_pos(records)

        # ALT
        if SVTYPE == 'BND':
            # Worry about strands later
            orient = remote = True
            ALT = [_Breakend(CHR2, END, orient, remote, 'N', True)]
        else:
            ALT = [_SV(SVTYPE)]

        # SVLEN for intra-chromosomal is -1
        if CHROM == CHR2:
            SVLEN = END - POS
        else:
            SVLEN = -1

        # QUAL, FILTER currently unused
        QUAL = None
        FILTER = None

        # INFO
        INFO = OrderedDict([
            ('SVTYPE', SVTYPE),
            ('CHR2', CHR2),
            ('END', END),
            ('SVLEN', SVLEN),
            ('CIPOS', CIPOS),
            ('CIEND', CIEND)])

        if 'STRANDS' in record.INFO:
            INFO['STRANDS'] = record.INFO['STRANDS']

        # Count number of samples with support from each caller
        sources = sorted(sources)
        source_support = defaultdict(set)
        for record in records:
            # Check genotype, allow './.' calls
            for call in record.samples:
                if call.data.GT not in NULL_GT[record.source]:
                    source_support[record.source].add(call.sample)
        for source in sources:
            INFO[source] = len(source_support[source])

        # Make FORMAT field from FORMAT fields of sources
        # Tag with source
        ORIG_IDS = [source + '_IDs' for source in sources]
        FORMAT = [(fmt, src) for src in sources for fmt in STD_FORMATS[src]]
        FORMAT = ['_'.join(fmt_src) for fmt_src in FORMAT]
        FORMAT = ':'.join(sources + ORIG_IDS + FORMAT)
        # PyVCF forces a 'GT' field as first in FORMAT, but we
        # don't call a merged genotype yet
        FORMAT = 'GT:' + FORMAT

        MergedCall = make_calldata_tuple(FORMAT.split(':'))

        # Aggregate any call information for each sample/source pairing
        source_calldata = defaultdict(list)
        original_IDs = defaultdict(list)
        for record in records:
            for call in record.samples:
                key = (call.sample, record.source)
                source_calldata[key].append(call.data)
                # Track IDs of variants with 0/0 call for evidence lookup
                original_IDs[key].append(record.ID)

        # Collapse multiple calls from same source
        sample_calldata = deque()
        for sample in samples:
            calldata = {}
            for source in sources:
                # List of original call IDs
                calldata[source + '_IDs'] = original_IDs[(sample, source)]

                # Merge data across multiple calls by same program
                data = source_calldata[sample, source]

                try:
                    data = cls.merge_calldata(data, source, MEDIANS[source])
                except:
                    import pdb
                    pdb.set_trace()
                calldata.update(data)

            # Call het if called by any source
            if any(calldata[source] for source in sources):
                calldata['GT'] = '0/1'
            else:
                calldata['GT'] = '0/0'

            FIELD = [calldata[field] for field in FORMAT.split(':')]
            sample_calldata.append(MergedCall(*FIELD))

        # sample_indexes
        sample_indexes = {sample: i for i, sample in enumerate(samples)}

        record = _Record(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,
                         sample_indexes, [])

        # Convert CallData tuples back to PyVCF _Call format
        sc = zip(samples, sample_calldata)
        calls = [_Call(record, sample, calldata) for sample, calldata in sc]
        record.samples = calls

        source_support = sorted(set([record.source for record in records]))
        record.INFO['SOURCES'] = source_support
        #  source_support = ','.join(source_support)

        return SVRecord(record)

    @staticmethod
    def merge_field_vals(vals, dtype, summary_fn):
        """
        Merge list of genotype field values

        Parameters
        ----------
        vals : list
        dtype : type
            type of value (or type of subvalue if values are lists/tuples)
        summary_fn : function
            Function to apply to numeric fields

        Returns
        -------
        val : type of input values
        """

        # Eliminate None values. If all None, return default null value
        vals = [val for val in vals if val is not None]
        if len(vals) == 0:
            if dtype == str:
                return ''
            else:
                return dtype(0)

        # Map merge_fields to each sub-item in tuple or list field
        if isinstance(vals[0], list) or isinstance(vals[0], tuple):
            ltype = type(vals[0])
            map_merge = partial(SVRecord.merge_field_vals, dtype=dtype,
                                summary_fn=summary_fn)
            val = [map_merge(sub_vals) for sub_vals in zip(*vals)]
            return ltype(val)

        # Join string fields together for now
        elif isinstance(vals[0], str):
            return ','.join(vals)

        # Take median or sum numerics
        elif isinstance(vals[0], float) or isinstance(vals[0], int):
            # check this
            vals = [v for v in vals if v is not None]
            if len(vals) == 0:
                return dtype(0)
            else:
                return dtype(summary_fn(vals))

        else:
            msg = 'Fields must be list, tuple, str, int, or float'
            raise PyVCFParsingError(msg)

    @staticmethod
    def null_genotypes(source):
        """
        Generate null genotype field values
        str = '', float = 0.0, int = 0

        Returns
        -------
        genotypes : dict
            Null values
        """
        genotypes = {}
        for fmt in STD_FORMATS[source].values():
            if fmt.type == 'String':
                genotypes[fmt.id] = ''
            elif fmt.type == 'Integer':
                if fmt.num > 1:
                    genotypes[fmt.id] = [0] * fmt.num
                else:
                    genotypes[fmt.id] = 0
            elif fmt.type == 'Float':
                if fmt.num > 1:
                    genotypes[fmt.id] = [0.0] * fmt.num
                else:
                    genotypes[fmt.id] = 0.0

        return genotypes

    @classmethod
    def merge_calldata(cls, calldata, source, MEDIANS):
        """
        Merge a list of CallData from vcf Records into consensus field values

        Strings are joined by a comma, numerics are summed or medianed

        Parameters
        ----------
        calldata : list of vcf.CallData

        Returns
        -------
        merged_data : dict
            {field: merged_value} for field in calldata attributes
        """

        # Make null call
        merged_data = {}
        called = [c for c in calldata if c.GT not in NULL_GT[source]]
        merged_data[source] = len(called)

        # Zero out genotype fields
        genotypes = cls.null_genotypes(source)

        # Return null call if no calls made in this
        if len(calldata) == 0:
            pass

        # If only one call, no need to merge
        elif len(calldata) == 1:
            call = calldata[0]
            for field in genotypes:
                genotypes[field] = getattr(call, field)

        # Merge data from overlapping calls
        else:
            for field in genotypes:
                if field in MEDIANS:
                    summary_fn = np.median
                else:
                    summary_fn = np.sum
                vals = [getattr(call, field) for call in calldata]

                dtype = STD_FORMATS[source][field].type
                if dtype == 'Float':
                    dtype = float
                elif dtype == 'Integer':
                    dtype = int
                elif dtype == 'String':
                    dtype = str

                val = cls.merge_field_vals(vals, dtype, summary_fn)
                genotypes[field] = val

        # Rename fields with a suffix indicating source caller
        for field in [k for k in genotypes.keys()]:
            genotypes[field + '_' + source] = genotypes.pop(field)

        merged_data.update(genotypes)
        return merged_data

    def __hash__(self):
        return id(self)
