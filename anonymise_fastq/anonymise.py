# Biopython imports
import Bio.SeqIO as seq_io
import Bio.Seq as seq

# Standard library import
import gzip
import sys
import os
import random


def count_fastq(file_path):
    """
    Read in FastQ entries from a single FastQ file and 
    output the number of reads and the total sequence length. 
    """
    count = 0
    total_len = 0
    with open_handler(file_path) as file:
        for fastQid, seq, qual in seq_io.QualityIO.FastqGeneralIterator(file):
            count += 1
            total_len += len(seq)
    return (count, total_len)


def get_all_fastQid(file_path, zipped=False):
    fastQids = []
    with open_handler(file_path) as file:
        for fastQid, seq, qual in seq_io.QualityIO.FastqGeneralIterator(file):
            fastQids.append(fastQid)
    return fastQids


def is_gzipped(file_path):
    extensions = os.path.basename(file_path).split('.')[1:]
    for ex in ['gzip', 'gz']:
        if ex in extensions:
            return True
    return False

def open_handler(file_path, filemode='rt'):
    if is_gzipped(file_path):
        return gzip.open(file_path, filemode)
    else:
        return open(file_path, filemode)


#########################
# Anonymising single input id
#########################

def parse_id(input_id, fmt):
    """
    Parse `input_id` into a dictionary {field name : value}
    base on the provided format `fmt`
    
    Parameters:
    -----------
    input_id :: String
      - a single ID of a FastQ file entry 
        WITHOUT the initial '@' character. 
      - e.g. @HWI-D00119:50:H7AP8ADXX:1:1101:1318:44446
    fmt :: String
      - a string of field names enclosed in angle brackets. 
      - e.g. (from illumina documentation 
      https://support.illumina.com/content/dam/illumina-support/
      help/BaseSpaceHelp_v2/Content/Vault/Informatics/
      Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm): 
      
      <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> 
      <read>:<is filtered>:<control number>:<sample number>
      
      - DO NOT include initial '@' character.
      
    """
    # provide an ending character to fmt
    input_id.rstrip()
    fmt.rstrip()
    fmt += '\n'
    input_id += '\n'
    record = {}
    i = 0
    id_index = 0
    while i < len(fmt) or '<' in fmt[i:]:
        if fmt[i] == '<':
            next_close_brac_index = fmt.index('>', i)
            field_name = fmt[i +1 : next_close_brac_index] 
            i += len(field_name) + 2
            
            if '<' in fmt[i:]: 
                next_open_brac_index = fmt.index('<', i)
                next_separator = fmt[i : next_open_brac_index]
                i += len(next_separator)
            else:
                next_separator = '\n'  
                i = len(fmt)
            
            
            if next_separator not in input_id[id_index:]:
                print(next_separator, input_id[id_index:], input_id, fmt)
                
            id_separator_index = input_id.index(next_separator, id_index)
            field_value = input_id[id_index : id_separator_index]
            id_index += len(field_value) + len(next_separator)
            
            if field_name not in record:
                record[field_name] = field_value
            else: 
                raise RuntimeError("Repeated field name \n     %s \n in format \n     %s" % (field_name, fmt))
        else: 
            raise RuntimeError("input_id: %s \n fmt: %s" % (input_id, fmt))
    return record


def anonymise_id(input_id, fmt, retained_fields=[], out_sep=':'):
    """
    #!! documentation see sample run below.
    """
    field_dict = parse_id(input_id, fmt)
    
    anonymised = ""
    for field_name in retained_fields:
        if field_name not in field_dict:
            raise RuntimeError("Field name '%s' not found in format \n   %s" % (field_name, fmt))
        anonymised += out_sep + str(field_dict[field_name]) 
        
    # extra processing here
    anonymised = str(abs(hash(input_id)))+ ("_%i" % random.randint(1, 10000000)) + anonymised
    return anonymised






################################
# Anonymising entire FastQ file
################################
def anonymise_fastq(file_path, fmt, retained_fields=[], out_sep=':'):
    """
    #!! see sample run below. 
    """
    with open_handler(file_path) as file:
        reads = seq_io.parse(file, "fastq")
        anonymised_reads = []
        for read in reads:
            anon = anonymise_id(read.id, fmt, retained_fields, out_sep)
            read.id = anon
            read.name = anon
            anonymised_reads.append(read)
    return anonymised_reads


## Streaming both input and output

def anonymise_process(infilepath, outfilepath, fmt, retained_fields=[], out_sep=':', outfilemode='wt'):
    """
    Stream from input fastq to output fastq. 
    #!!! Change to using zipfile output writer bgzf!!
    """
    with open_handler(infilepath) as infile:
        with open_handler(outfilepath, filemode=outfilemode) as outfile:
            count = 0
            reads = seq_io.parse(infile, "fastq")
            for read in reads:
                anon = anonymise_id(read.id, fmt, retained_fields, out_sep)
                read.id = anon
                read.name = anon
                seq_io.write(read, outfile, "fastq")
                count += 1
    return count




