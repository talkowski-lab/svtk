from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libcalignmentfile cimport AlignmentFile

cpdef bint is_excluded(AlignedSegment read):
    cdef bint exclude = (read.is_unmapped or
                         read.mate_is_unmapped or
                         read.is_secondary or
                         read.is_duplicate or
                         read.is_supplementary)
    return exclude

cpdef bint is_soft_clipped(AlignedSegment read):
    return (read.cigartuples[0][0] == 4) ^ (read.cigartuples[-1][0] == 4)

def collect_splits(bam):
    cdef AlignedSegment read

    for read in bam:
        if is_excluded(read):
            continue
        if is_soft_clipped(read):
            yield read

cpdef int get_clip_direction(AlignedSegment read):
    """
    Calculate split direction based on CIGAR ops - (LEFT, RIGHT, MIDDLE)

    Parameters
    ----------
    read : pysam.AlignedSegment

    Returns
    -------
    direction : str [RIGHT,LEFT,MIDDLE, NO_CLIP]
        Direction of soft clip
    """

    cdef int first_op = read.cigartuples[0][0]
    cdef int last_op = read.cigartuples[-1][0]

    if first_op == 4 and last_op == 4:
        return 0  # 'MIDDLE'
    elif first_op == 4:
        return 1  # 'LEFT'
    elif last_op == 4:
        return 2  # 'RIGHT'
    else:
        return 3  # 'NO_CLIP'

cpdef int get_clip_position(AlignedSegment read, int direction):
    cdef int pos = read.pos

    cdef int op
    cdef int length

    cdef int LEFT = 1

    if direction == LEFT:
        for op, length in read.cigartuples:
            if op == 0:
                pos += length

    return pos
