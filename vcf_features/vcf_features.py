
import pandas as pd

import intervaltree

import csv 
import gzip
import os




def gene_variants(df_vcf, chrom_interval_tree):
    """
    #!!! documentation. 
    """
    sample_name = df_vcf.columns[-1] #!!! Need to make this more general / robust
    gene_variant_record = {}
    no_chrom = []
    no_intersections = []
    
    for index, variant in df_vcf.iterrows():
        chrom = variant[0]        
        if chrom not in chrom_interval_tree:
            no_chrom.append(index)
            #print("Chromosome %s not in supplied gene list. \nCurrent variant: \n%s" % (chrom, variant))
            continue
        
        intersections = chrom_interval_tree[chrom][variant["POS"]]
        if not intersections: # no intersections
            no_intersections.append(index)
            #print("Variant not contained in genes provided.\nCurrent variant: \n%s" % variant)
            continue
            
        for interval in intersections:
            gene_name = interval.data
            if gene_name not in gene_variant_record:
                gene_variant_record[gene_name] = []
            #genotype = get_genotype(variant[sample_name], fmt=variant["FORMAT"])
            gene_variant_record[gene_name].append(variant)
    return gene_variant_record, df_vcf.iloc[no_chrom, :], df_vcf.iloc[no_intersections, :] #!!! need deep copy?
        





## Utilities

###########
def get_genotype(sample_string, fmt=''):
    """
    #!!! Docs!
    """
    if not fmt:
        gt_index = 0
    else:
        fields = fmt.split(':')
        if "GT" in fields:
            gt_index = fields.index("GT")
        else:
            raise RuntimeError("No genotype tag detected.")
    return sample_string.split(':')[gt_index]

############
def num_header_row(vcfpath):
    """
    #!!! Add documentation
    #!!! Check vcf docs to verity this is robust
         e.g. first row might not have '#' even if it has header. 
    """
    # find out number of header rows
    with open_handler(vcfpath, 'rt') as infile: 
        num = 0
        for row in infile:
            if row[0] != '#':
                num += 1
            else:
                break
    return num

def vcf_to_df(vcfpath):
    """
    #!!! Docs!
    """
    # num = num_header_row(vcfpath)
    with open_handler(vcfpath, 'rt') as infile:
        header = []
        n_header = 0
        for row in infile:
            row = row.rstrip()
            if row[0] == '#':
                n_header += 1
                header.append(row)
            else:
                break
        
        if '\t' in header[-1]: 
            names = header[-1].split('\t')
        
    df_vcf = pd.read_csv(open_handler(vcfpath, "rt"), 
                         skiprows=n_header, 
                         sep='\t', 
                         names=names)
    return header, df_vcf

def bed_to_df(bedfilepath, header=None):
    """
    Given a path name to a BED file, 
    parse the first 4 columns into 
    a pandas.DataFrame object. 
    """
    df_bed = pd.read_csv(open_handler(bedfilepath, 'rt'), 
                     delimiter='\t', header=header)
    return df_bed


def build_chrom_interval_tree(bedfilepath):
    """
    Build an interval tree for each chromosome
    in the BED file. 
    """
    with open(bedfilepath) as bedfile:
        tree_dict = {}
        for row in csv.reader(bedfile, delimiter='\t'):
            chrom, start, end, name = row[:4]
            start = int(start)
            end = int(end)
            if chrom not in tree_dict:
                tree_dict[chrom] = intervaltree.IntervalTree()
            tree_dict[chrom][start : end] = name
    return tree_dict

##############

def is_gzipped(file_path):
    """
    Detect if the file is gzipped. 
    A path name is gzipped if and only if 
      - it is dot separated.
      - and one of dot separated field is 'gz' or 'gzip'. 
    """
    extensions = os.path.basename(file_path).split('.')[1:]
    for ex in ['gzip', 'gz']:
        if ex in extensions:
            return True
    return False

def open_handler(file_path, filemode='rt'):
    """
    Return a filehandle base on whether 
    it is a regular file or a gzipped file. 
    """
    if is_gzipped(file_path):
        return gzip.open(file_path, filemode)
    else:
        return open(file_path, filemode)

    