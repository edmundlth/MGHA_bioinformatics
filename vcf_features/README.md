# VCF Features Analysis
----
This package provides tools to extract features from VCF files and produce analysis report for those features. 


---
## Features provided
---
At the moment, the only feature provided is the mutation load of a specified list of genes. See pseudocode below. 


---
## Analysis Provided
---
We generate a HTML report containing basic statistics of the provided VCFs and generate (2 and 3 dimensional) PCA plots for the mutation load. 

---
# Code
---

Given 
  * a list of genes (provided in a bed file of their exon position for example).
  * a list of variant information (provided in a vcf file). 
compute the mutation load for each gene. 

Ultimately we want to output a dataframe with the following rows and columns: 

|         | Gene1 | Gene2 | ... | GeneM |
|---------|-------|-------|-----|-------|
| Sample1 | load  |       |     |       |
| Sample2 |       |       |     |       |
| ...     |       |       |     |       |
| SampleN |       |       |     |       |


## Pseudocode
---
Input: 
 1. a bed file specifying genomic intervals labelled by their gene names. (overlapping intervals allowed)
 2. a single sample vcf file. 

Output: 
 1. a list indexed by gene names containing the mutation load of each gene. 

Steps: 
 1. build an interval tree from bedfile (for each chromosome?) 
 2. for each variant in the vcf file
   1. query the interval tree as to what intervals the variant intersects.
   2. for each intervals that intersects the variant 
     1. update the number of mutation of the gene associated to that interval by add `1` for each allele that is different from the reference. 

---
# Example run
---
### On commandline
----
1. Use `$ vcf_features --help` to see documentation and input options. 
2. Sample input: 
```
$ vcf_features \
    --bed ./example_gene_exons.bed \
    --outdir analysis_output \
    --vcf \
    /path/to/data/sample1.vcf.gz \
    /path/to/data/sample2.vcf.gz \
    /path/to/data/sample3.vcf.gz \
    /path/to/data/sample4.vcf.gz \
    /path/to/data/sample5.vcf \
    /path/to/data/sample6.vcf \
    /path/to/data/sample7.vcf.gz \
    --label \
    Class1 \
    Class1 \
    Class2 \
    Class2 \
    Class2 \
    Class3 \
    Class3 \


! open analysis_output/report.html
```
