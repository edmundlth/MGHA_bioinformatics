#!/usr/bin/env python

from distutils.core import setup

setup(name="MGHA_bioinformatics",
      version="0.2.0",
      author="Edmund Lau",
      author_email="edmundlth95@gmail.com",
      packages=["anonymise_fastq", "vcf_features"],
      package_dir={
        "anonymise_fastq" : "anonymise_fastq", 
        "vcf_features" : "vcf_features"},
      package_data={"vcf_features" : ["template_report.html"]},
      entry_points={
             "console_scripts":[
               "anonymise_fastq = anonymise_fastq.main:main", 
               "vcf_features = vcf_features.main:main"
             ]
            },
      url="https://github.com/edmundlth/MGHA_bioinformatics",
      install_requires=[
        "Biopython >= 1.73", 
        "matplotlib>=2.0",
        "numpy>=1.14.5",
        "pandas>=0.23.1",
        "scikit-learn>=0.19.1",
        "intervaltree>=3.0.2"
      ],
     )
