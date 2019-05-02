# Genomic File Edit
----
## Problem Statement
We would like to have tools to efficiently make changes to genomic files that we collected.

## Detailed Statement:
The parameters are:

 - By "Changes" we mean, at the moment anonymisation.  
 - By "Genomic files" we mostly mean FastQ, BAM, VCF. These files might be from different stages of bioinformatics pipeline, so they might require different ways of handling.
 - By "Efficiently" we mean time and memory efficient. Most of the algorithm should be O(n) in time, but we would like to O(1) in memory, i.e. we would like to process these files by streaming from input to output. Furthermore, for files that might be compressed (gzipped or bgzipped), we would like to avoid a decompression step. It would be even better if we can stream into a compressed version of the file already.
 - By "Tools" : probably an executable Python script or a pip-installable commandline tool. We might be using these in Docker.



---
# Example run
---
1. Use `$ python utility.py --help` to see documentation and input options.
2. Sample input:
```
$ python utility.py  \
    --old "OLD STRING TO BE REPLACED WITH NEW STRING"
    --new "NEW STRING"
    --infile "/path/to/input/file" # can be gzip or bgzip files
    --outfile "/path/to/output/file" # will overwrite. use suitable extensions.
    --outfile_compress gzip
```
