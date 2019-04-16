#!/usr/bin/env python

from distutils.core import setup

setup(name="anonymise_fastq",
      version="1.0.0",
      author="Edmund Lau",
      author_email="edmundlth95@gmail.com",
      packages=["anonymise_fastq"],
      package_data={},
      entry_points={
             "console_scripts":["anonymise_fastq = anonymise_fastq.main:main"]
            },
      url="https://github.com/edmundlth/MGHA_bioinformatics",
      install_requires=["Biopython >= 1.73"],
     )
