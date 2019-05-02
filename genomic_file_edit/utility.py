#!/usr/bin/env python
"""
#!! Docs
"""
import os
import gzip
from typing import TextIO, Optional
from Bio.bgzf import BgzfReader, BgzfWriter
import warnings
from argparse import ArgumentParser

def vcf_edit_header_by_line(old: str,
                            new: str,
                            infilename: str,
                            outfilename: str,
                            outfile_compress: Optional[str] = None) -> None:
    """
    This mimics the vcf_edit function writen in vcf_edit.py in
    github/genovic/gos-ensign.

    The differences are
      - we allow gzipped input file. The function will detect gzip-file automatically.
      - we allow streaming into compressed output file. Compression mode can be either gzip or bgzip.
    """
    if infilename == outfilename:
        raise Exception("Disallow in-place changes. Use different outfilename")

    with infile_handler(infilename) as infilehandle, \
            outfile_handler(outfilename,
                            compression=outfile_compress) as outfilehandle:
        for line in infilehandle:
            if line.startswith('#'):
                new_line = line.replace(old, new)
            else:
                new_line = line
            outfilehandle.write(new_line)
    return


def is_gzip(filepath: str) -> bool:
    """
    Check if a the file specified by filepath
    is a gzipped file. We do this by checking
    whether the first two bytes is the magic
    numbers included in the file header.
    We check them bytewise to avoid issue with
    endianness.

    Note that this wouldn't detect if the gzip
    file is a concatenation of multiple "members".
    See gzip specification at:
     https://www.ietf.org/rfc/rfc1952.txt
    """
    if not os.path.isfile(filepath):
        warnings.warn("The file %s does not exist" % filepath)
        return False

    with open(filepath, 'rb') as filehandle:
        byte1 = filehandle.read(1)
        byte2 = filehandle.read(1)
    if byte1 == b'\x1f' and byte2 == b'\x8b':
        return True
    else:
        return False


def infile_handler(filepath: str) -> TextIO:
    """
    Detect if the file specified by `filepath` is gzip-compressed
    and open the the file in read mode using appriate open handler.
    """
    if is_gzip(filepath):
        return gzip.open(filepath, "rt")
    else:
        return open(filepath, "rt")


def outfile_handler(filepath: str,
                    compression: Optional[str] = None) -> TextIO:
    """
    Return a file handle in write mode using the appropriate
    handle depending on the compression mode.
    Valid compression mode:
        compress = None | "None" | "gzip" | "gz" | "bgzip" | "bgz"
    If compress = None or other input, open the file normally.
    """
    if os.path.isfile(filepath):
        warnings.warn("Overwriting the existing file: %s" % filepath)

    if compression is None:
        return open(filepath, mode="wt")
    elif type(compression) == str:
        if compression.lower() in ["gzip", "gz"]:
            return gzip.open(filepath, mode="wt")
        elif compression.lower() in ["bgzip", "bgz"]:
            return BgzfWriter(filepath)
        elif compression.lower() == "none":
            return open(filepath, mode="wt")
    else:
        raise Exception("`compression = %s` invalid." % str(compression))


#############################################################################
# Commandline argument parser and main function
#############################################################################


def parse_args():
    """Replace old text with new text in the header of a VCF file"""
    parser = ArgumentParser(description="Replace old text with new"
                            " text in the header of a VCF file")
    parser.add_argument("--old",
                        required=True,
                        type=str,
                        help="old string (to be replaced)")
    parser.add_argument("--new",
                        required=True,
                        type=str,
                        help="new string (to replace old)")
    parser.add_argument("--outfile",
                        required=True,
                        type=str,
                        help="output VCF file path")
    parser.add_argument("--infile",
                        required=True,
                        type=str,
                        help="input VCF file path")
    parser.add_argument("--outfile_compress",
                        type=str,
                        default=None,
                        help="Specify whether the output file should be"
                             "compressed. Valid arguments")
    return parser.parse_args()


def main():
    user_inputs = parse_args()
    if 'gz' in user_inputs.outfile.split('.') \
    and not user_inputs.outfile_compress:
        warnings.warn("Outfilename %s has '.gz' extension"
                      " but user specifed compression method is %s"
                      % (user_inputs.outfile, str(outfile_compress)))

    vcf_edit_header_by_line(user_inputs.old,
             user_inputs.new,
             user_inputs.infile,
             user_inputs.outfile,
             user_inputs.outfile_compress)
    return


if __name__ == '__main__':
    main()
