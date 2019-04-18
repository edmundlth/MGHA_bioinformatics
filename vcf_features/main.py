
import logging
import sys
import os
import argparse

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA


from .vcf_features import gene_variants, vcf_to_df, build_chrom_interval_tree, get_genotype





def main():
    """
    Main
    """
    # parse input
    parser = argparser()
    user_inputs = parser.parse_args()
    if not user_inputs.label:
        user_inputs.label = ["None"] * len(user_inputs.vcf)


    
    # Initialise recording
    output_record = {
    "BED file"  : user_inputs.bed,
    "VCF files" : user_inputs.vcf,
    "labels"    : user_inputs.label
    }
    

    # parse input files
    chrom_interval_tree = build_chrom_interval_tree(user_inputs.bed)
    vcf_info_dict = {}
    record = {}
    label_dict = {}
    for i in range(len(user_inputs.vcf)):
        vcffilepath = user_inputs.vcf[i]
        vcf_label = user_inputs.label[i]
        vcf_name = os.path.basename(vcffilepath)
        label_dict[vcf_name] = vcf_label

        vcf_header, df_vcf = vcf_to_df(vcffilepath)
        gene_variant_record, df_no_chrom, df_no_intersection = gene_variants(df_vcf, chrom_interval_tree)
        vcf_info_dict[vcf_name] = (gene_variant_record, df_no_chrom, df_no_intersection, vcf_label, df_vcf.shape[0])

        record[vcf_name] = {}
        for gene in gene_variant_record:
            genotypes = [get_genotype(var[-1], var[-2]) for var in gene_variant_record[gene]]
            record[vcf_name][gene] = sum(
                sum(allele.isdigit() and int(allele) > 0 
                    for allele in gt.split('/')
                    )
                for gt in genotypes
                )
    df_load = pd.DataFrame.from_dict(record, orient='index', dtype=int).sort_index(axis=1).fillna(value=0)

    output_record["df_load"] = df_load
    output_record["vcf_info_dict"] = vcf_info_dict
    output_record["chrom_interval_tree"] = chrom_interval_tree


    # Do PCA on mutation load data (contained in df_load).
    X_reduced = PCA(n_components=3).fit_transform(df_load)
    y = pd.Series([label_dict[name] for name in df_load.index])


    # Create the output directory if it hasn't already existed. 
    if not os.path.isdir(user_inputs.outdir):
        os.mkdir(user_inputs.outdir)


    # Plot PCA. 
    fig2d, ax2d = plot_pca(X_reduced, y, n_dim=2, figsize=(5, 5))
    fig2dname = os.path.join(user_inputs.outdir, "load_pca2d.png")
    fig2d.savefig(fig2dname, format='png', bbox='tight')
    output_record["load_pca2d_filename"] = "load_pca2d.png"

    fig3d, ax3d = plot_pca(X_reduced, y, n_dim=3, figsize=(5, 5))
    fig3dname = os.path.join(user_inputs.outdir, "load_pca3d.png")
    fig3d.savefig(fig3dname, format='png', bbox='tight')
    output_record["load_pca3d_filename"] = "load_pca3d.png"


    # Generate a HTML report. 
    report_html = write_report(output_record)
    with open(os.path.join(user_inputs.outdir, "report.html"), 'w') as outfile:
        outfile.write(report_html)
    return 0



def argparser():
    """
    Parse commandline arguments into Namespace() object. 
    This will provide help message for users when they 
    input `$ <command name> --help`. 
    """
    parser = argparse.ArgumentParser(
            description="Extract analysable features from provide VCF files.")
    parser.add_argument("--vcf", metavar="FILENAME", 
                        type=str, 
                        nargs='*', 
                        required=True,
                        help="Input normalised VCF files. "
                             "Files can be gzipped. "
                             "DO NOT input gVCF.")
    parser.add_argument("--bed", metavar="FILENAME", 
                        type=str, 
                        required=True,
                        help="A BED file specifying genomic regions of interest."
                             "The file can be gzipped. "
                             "Each reagion is accompanied by a name "
                             "(typically the name of associated gene). ")
    parser.add_argument("--label", metavar="LABEL", 
                        type=str,
                        nargs='*',
                        default=None,
                        help="Class label (e.g. by disease) for each vcf files. "
                             "Must be provided in the same order as the vcf files. "
                             "Default: None")
    parser.add_argument("--outdir", metavar="DIR", 
                        type=str, 
                        default="./analysis_output/", 
                        help="Output directory. Default to './analysis_output/'.")
    return parser


def plot_pca(X_reduced, y, n_dim=2, figsize=(10, 10)):
    """
    Given a the PCA reduced feature vectors `X_reduced`
    together with their class labels `y`, 
    generate a 2 or 3 dimensional scatter plot. 
    """
    fig = plt.figure(1, figsize=figsize)
    if n_dim == 3:
        ax = fig.add_subplot(111, projection='3d')
    elif n_dim == 2:
        ax = fig.add_subplot(111)
    else:
        raise RuntimeError("Invalid number of dimension.")

    n_classes = len(set(y))
    for index, label, color in zip(range(n_classes), 
                               set(y), 
                               plt.cm.Set1(np.linspace(0, 1, num=n_classes))):   
        row_index = y.index[y == label]
        if n_dim == 3:
            ax.scatter(X_reduced[row_index, 0], 
                       X_reduced[row_index, 1], 
                       X_reduced[row_index, 2], 
                       color=color,
                       label=label,
                       s=40)
            ax.set_zlabel("PC3")
        elif n_dim == 2:
            ax.scatter(X_reduced[row_index, 0], 
                       X_reduced[row_index, 1], 
                       color=color,
                       label=label,
                       s=40)
    ax.set_title("First %i Principal Components" % n_dim)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(loc='best')
    return fig, ax



############################
def write_report(output_record):
    """
    Given the output of various intermediate data
    of `main()` recorded in the dictionary `output_record`, 
    generate a HTML string as report. 
    """
    formating_data = {
    "bedfilepath"           : output_record["BED file"], 
    "label_file_html_table" : label_file_html_table(output_record["VCF files"], output_record["labels"]),
    "vcf_info_html_table"   : vcf_info_html_table(output_record["vcf_info_dict"]), 
    "fig_pca2d"             : """<img src="%s" alt="2D PCA plot">""" % output_record["load_pca2d_filename"],
    "fig_pca3d"             : """<img src="%s" alt="3D PCA plot">""" % output_record["load_pca3d_filename"],
    }

    this_dir, this_filename = os.path.split(__file__)
    template_path = os.path.join(this_dir, "template_report.html")
    with open(template_path, 'r') as infile:
        template = infile.read()
        output = template.format(**formating_data)
    return output

def label_file_html_table(vcf_files, labels):
    """
    Given a list of vcf file paths in `vcf_files`, 
    and their corresponding class label in `labels`, 
    generate a label:vcf_file HTML table. 
    """
    outstring = ""
    for vfile, label in zip(vcf_files, labels):
        outstring += """<tr>
        <td> {} </td>
        <td> {} </td>
        </tr>\n
        """.format(str(label), vfile)
    return outstring

def vcf_info_html_table(vcf_info_dict):
    """
    Generate a HTML table summarising 
    vcf information collected in `vcf_info_dict`. 
    """
    outstring = ""
    for vcf_name in vcf_info_dict: 
        gene_variant_record, df_no_chrom, df_no_intersection, vcf_label, n_variant = vcf_info_dict[vcf_name]
        n_no_chrom = df_no_chrom.shape[0]
        n_no_intersection = df_no_intersection.shape[0]

        outstring += """
        <tr>
        <td> {} </td>
        <td> {} </td>
        <td> {} </td>
        <td> {} </td>
        <td> {} </td>
        </tr>\n
        """.format(vcf_name, 
                   vcf_label, 
                   n_variant,
                   n_no_chrom,
                   n_no_intersection)
    return outstring



if __name__ == "__main__":
    main()
