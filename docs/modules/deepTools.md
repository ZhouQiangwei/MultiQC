---
Name: deepTools
URL: http://deeptools.readthedocs.io/en/latest/index.html
Description: >
deepTools is a suite of python tools particularly developed
for the efficient analysis of high-throughput sequencing data,
such as ChIP-seq, RNA-seq or MNase-seq.
---

Supported commands:

* `plotCorrelation`
* `plotPCA`
* `plotFingerprint`
* `computeGCBias`
* `plotCoverage`

since the log file names are set by the user and the files don't contain the sample name, MultiQC will find the files through:

```yaml
deepTools:
    Corr:
        fn: '*Corr*tab'
    PCA:
        fn: '*PCA*tab'
    fgpr:
        fn: '*fgpr*tab'
    gc_bias:
        fn: '*gcFreq.txt'
    coverage:
        fn: '*Cov*tab'
```
