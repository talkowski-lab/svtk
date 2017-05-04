# SVTools

Utilities for manipulating structural variation calls.

```
$ svtools -h
SVTools: Utilities for manipulating structural variation

usage: svtools [-h] <subcommand> [options]

[ Preprocessing ]
    standardize   Convert SV calls to a standardized format.

[ Algorithm integration ]
    vcfcluster    Cluster SV calls from a list of VCFs. (Generally PE/SR.)
    bedcluster    Cluster SV calls from a BED. (Generally depth.)

[ Uncategorized or not yet implemented ]
    bincov        Calculate normalized genome-wide depth of coverage.
    rdtest        Calculate comparative coverage statistics at CNV sites.

optional arguments:
  -h, --help  show this help message and exit
```
