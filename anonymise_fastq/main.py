"""
#!!! Need documentation. See README.md for now.
"""




import logging
import sys
import argparse

from anonymise import anonymise_process


# Uncomment for memory_profiler to work on this script
#@profile
def main():
	"""
	Main
	"""
	parser = argparser()
	user_inputs = parser.parse_args()
	anonymise_process(user_inputs.infilepath, 
		              user_inputs.outfilepath, 
    	              fmt=user_inputs.idformat, 
    	              retained_fields=user_inputs.retained_fields,
    	              out_sep=':',
    	              outfilemode='wt')
	return 0



def argparser():
    """
    Parse commandline arguments into Namespace() object.
    """
    parser = argparse.ArgumentParser(
            description="Anonymise FastQ ID string.")
    parser.add_argument("--infilepath", metavar="FILENAME", 
                        type=str, 
                        help="File path to the input `.fastq` or `.fastq.gz` file.")
    parser.add_argument("--outfilepath", metavar="FILENAME", 
    	                type=str, 
    	                default="anonymise_fastq_output.fastq",
    	                help="Name of the output file. "
    	                     "Use `.fastq.gz` for zipped output. "
    	                     "Default: ")
    parser.add_argument("--idformat", metavar="FORMAT", 
    	                type=str, 
    	                default="<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>/<sense>", 
    	                help="A string specifying the format "
    	                     "of the ID string for each read. "
    	                     "Each field is given a unique name enclosed by angle bracket. "
                             "Default: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>/<sense>")
    parser.add_argument("--retained_fields", metavar="FIELDS", 
    	                type=str,
    	                nargs='*',
    	                default=[],
    	                help="A list of field names as specified in the ID format string "
    	                     "for which the anonymising process would retain. "
    	                     "Default: Empty")
    return parser

if __name__ == "__main__":
    main()
