# SVTools

Utilities for manipulating structural variation calls.

## Installation

```
$ git clone https://github.com/talkowski-lab/svtools.git
$ cd svtools
$ pip install -e .
```

## Available commands

```
SVTools: Utilities for manipulating structural variation

usage: svtools [-h] <subcommand> [options]

[ Preprocessing ]
    standardize   Convert SV calls to a standardized format.

[ Algorithm integration ]
    vcfcluster     Cluster SV calls from a list of VCFs. (Generally PE/SR.)
    bedcluster     Cluster SV calls from a BED. (Generally depth.)

[ Statistics ]
    count-svtypes  Count instances of each svtype in each sample in a VCF

[ Read-depth analysis ]
    bincov         Calculate normalized genome-wide depth of coverage.
    rdtest*        Calculate comparative coverage statistics at CNV sites.

[ Split-read analysis ]
    sr-count*      Pile up clipped read counts genome-wide
    sr-test        Calculate enrichment of clipped reads at SV breakpoints

* Not yet implemented

optional arguments:
  -h, --help  show this help message and exit
```
