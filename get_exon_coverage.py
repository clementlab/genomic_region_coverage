import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import subprocess

from collections import defaultdict
from matplotlib import pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages

def get_genes_from_file(gene_exon_list_file):
    with open(gene_exon_list_file, 'r') as fh_geneList:
        gene_list = [line.strip() for line in fh_geneList]
    print('Read ' + str(len(gene_list)) + ' genes from the gene list file.')
    return gene_list

def get_exons_from_refgene(gene_list, ref_gene_file, chr_index, strand_index, tx_start_index, tx_end_index, cds_start_index, cds_end_index, exon_index, intron_index, gene_id_index, exon_status_index):
    gene_exons = defaultdict(list) # gene_id > list(exons) where exons is [starts: dict of starts, ends: dict of ends]
    with open(ref_gene_file, 'r') as fh_refGene:
        for line in fh_refGene:
            line_arr = line.strip().split("\t")
            line_chr = line_arr[chr_index]
            line_start = int(line_arr[tx_start_index])
            line_end = int(line_arr[tx_end_index])
            line_strand = line_arr[strand_index]
            line_exons = []
            if line_arr[exon_index] != "":
                line_exons = [int(x) for x in line_arr[exon_index].rstrip(',').split(",")]
            line_introns = []
            if line_arr[intron_index] != "":
                line_introns = [int(x) for x in line_arr[intron_index].rstrip(',').split(",")]
            line_cds_start = int(line_arr[cds_start_index])
            line_cds_end = int(line_arr[cds_end_index])
            line_exon_count = len(line_exons)
            line_gene_id = line_arr[gene_id_index]
            
            if line_gene_id not in gene_list:
                continue

            this_gene_exons = []

            # if on the watson strand
            if line_strand == "+":
                # process introns/exons
                for i in range(line_exon_count):
                    this_gene_exons.append((line_chr,line_exons[i],line_introns[i]))

            # if on the crick strand
            else:
                # process introns/exons
                for i in range(line_exon_count-1, -1, -1):
                    this_gene_exons.append((line_chr,line_exons[i],line_introns[i]))

            gene_id_key = (line_gene_id, line_chr + ":" + str(line_start) + "-" + str(line_end))
            gene_exons[gene_id_key] = this_gene_exons

    print('Finished reading refGene.txt. Got exon information for ' + str(len(gene_exons)) + ' genes.')
    
    return gene_exons

def read_command_output(command):
    """
    Runs a shell command and returns an iter to read the output

    Args:
        command: shell command to run

    Returns:
        iter to read the output
    """

    p = subprocess.Popen(command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,shell=True,
#            encoding='utf-8',universal_newlines=True)
            universal_newlines=True,
            bufsize=-1) #bufsize system default
    return iter(p.stdout.readline, b'')

def get_exon_coverage(bam_input, exon_info):
    """
    Gets coverage information for a given exon in a BAM file
    
    Args:
        bam_input: path to the BAM file
        exon_info: tuple of (chr, start, end) for the exon
        
    Returns:
        dict of coverage information for the exon
    """

    nan_value = np.nan
    per_nucleotide_depth_command = f'samtools depth -a -r {exon_info[0]}:{str(exon_info[1])}-{str(exon_info[2])} {bam_input}'
    per_nucleotide_depth_lines = subprocess.check_output(per_nucleotide_depth_command, shell=True).decode('utf-8').split("\n")
    per_nucleotide_depth_values = np.array([int(x.split("\t")[2]) for x in per_nucleotide_depth_lines if x != ""])

    max_per_nucleotide_depth = np.max(per_nucleotide_depth_values) if len(per_nucleotide_depth_values) > 0 else nan_value
    mean_per_nucleotide_depth = np.mean(per_nucleotide_depth_values) if len(per_nucleotide_depth_values) > 0 else nan_value
    median_per_nucleotide_depth = np.median(per_nucleotide_depth_values) if len(per_nucleotide_depth_values) > 0 else nan_value
    min_per_nucleotide_depth = np.min(per_nucleotide_depth_values) if len(per_nucleotide_depth_values) > 0 else nan_value

    #bq30_per_nucleotide_depth_command = f'samtools depth --min-BQ 30 -a -r {exon_info[0]}:{str(exon_info[1])}-{str(exon_info[2])} {bam_input}'
    bq30_per_nucleotide_depth_command = f'samtools depth -q 30 -a -r {exon_info[0]}:{str(exon_info[1])}-{str(exon_info[2])} {bam_input}'
    bq30_per_nucleotide_depth_lines = subprocess.check_output(bq30_per_nucleotide_depth_command, shell=True).decode('utf-8').split("\n")
    bq30_per_nucleotide_depth_values = np.array([int(x.split("\t")[2]) for x in bq30_per_nucleotide_depth_lines if x != ""])

    bq30_max_per_nucleotide_depth = np.max(bq30_per_nucleotide_depth_values) if len(bq30_per_nucleotide_depth_values) > 0 else nan_value
    bq30_mean_per_nucleotide_depth = np.mean(bq30_per_nucleotide_depth_values) if len(bq30_per_nucleotide_depth_values) > 0 else nan_value
    bq30_median_per_nucleotide_depth = np.median(bq30_per_nucleotide_depth_values) if len(bq30_per_nucleotide_depth_values) > 0 else nan_value
    bq30_min_per_nucleotide_depth = np.min(bq30_per_nucleotide_depth_values) if len(bq30_per_nucleotide_depth_values) > 0 else nan_value

    #mq20_per_nucleotide_depth_command = f'samtools depth --min-MQ 20 -a -r {exon_info[0]}:{str(exon_info[1])}-{str(exon_info[2])} {bam_input}'
    mq20_per_nucleotide_depth_command = f'samtools depth -Q 20 -a -r {exon_info[0]}:{str(exon_info[1])}-{str(exon_info[2])} {bam_input}'
    mq20_per_nucleotide_depth_lines = subprocess.check_output(mq20_per_nucleotide_depth_command, shell=True).decode('utf-8').split("\n")
    mq20_per_nucleotide_depth_values = np.array([int(x.split("\t")[2]) for x in mq20_per_nucleotide_depth_lines if x != ""])

    mq20_max_per_nucleotide_depth = np.max(mq20_per_nucleotide_depth_values) if len(mq20_per_nucleotide_depth_values) > 0 else nan_value
    mq20_mean_per_nucleotide_depth = np.mean(mq20_per_nucleotide_depth_values) if len(mq20_per_nucleotide_depth_values) > 0 else nan_value
    mq20_median_per_nucleotide_depth = np.median(mq20_per_nucleotide_depth_values) if len(mq20_per_nucleotide_depth_values) > 0 else nan_value
    mq20_min_per_nucleotide_depth = np.min(mq20_per_nucleotide_depth_values) if len(mq20_per_nucleotide_depth_values) > 0 else nan_value

    total_reads_processed = 0
    count_reads_unaligned = 0
    count_reads_not_primary_alignment = 0
    count_reads_mate_unaligned = 0
    count_reads_mate_on_other_chr = 0
    count_reads_mate_gt_10kb_away = 0
    count_reads_soft_clipped = 0
    count_reads_hard_clipped = 0
    count_reads_fw_strand = 0

    for line in read_command_output(f'samtools view {bam_input} {exon_info[0]}:{str(exon_info[1])}-{str(exon_info[2])}'):
        if line.strip() == "": break
        total_reads_processed += 1

        line_els = line.strip().split("\t")

        read_is_not_primary_alignment = int(line_els[1]) & 0x100
        if read_is_not_primary_alignment:
            count_reads_not_primary_alignment += 1

        line_unmapped = int(line_els[1]) & 0x4
        if line_unmapped:
            count_reads_unaligned += 1

        mate_unmapped = int(line_els[1]) & 0x8
        if mate_unmapped:
            count_reads_mate_unaligned += 1

        is_rv = int(line_els[1]) & 0x10
        if not is_rv:
            count_reads_fw_strand += 1

        cigar = line_els[5]
        if 'S' in cigar:
            count_reads_soft_clipped += 1
        if 'H' in cigar:
            count_reads_hard_clipped += 1

        mate_chr = line_els[6]
        mate_pos = int(line_els[7])

        if mate_chr != "=":
            count_reads_mate_on_other_chr += 1

        if mate_chr != "=" and mate_chr != line_els[2]:
            count_reads_mate_on_other_chr += 1

        if mate_chr == line_els[2] and abs(mate_pos - int(line_els[3])) > 10000:
            count_reads_mate_gt_10kb_away += 1 

    pct_reads_unaligned = count_reads_unaligned / total_reads_processed if total_reads_processed > 0 else nan_value
    pct_reads_not_primary_alignment = count_reads_not_primary_alignment / total_reads_processed if total_reads_processed > 0 else nan_value
    pct_reads_mate_unaligned = count_reads_mate_unaligned / total_reads_processed if total_reads_processed > 0 else nan_value
    pct_reads_mate_on_other_chr = count_reads_mate_on_other_chr / total_reads_processed if total_reads_processed > 0 else nan_value
    pct_reads_mate_gt_10kb_away = count_reads_mate_gt_10kb_away / total_reads_processed if total_reads_processed > 0 else nan_value
    pct_reads_soft_clipped = count_reads_soft_clipped / total_reads_processed if total_reads_processed > 0 else nan_value
    pct_reads_hard_clipped = count_reads_hard_clipped / total_reads_processed if total_reads_processed > 0 else nan_value
    pct_reads_fw_strand = count_reads_fw_strand / total_reads_processed if total_reads_processed > 0 else nan_value

    exon_size = exon_info[2] - exon_info[1]
    per_base_reads_processed = total_reads_processed / exon_size
    per_base_unaligned = count_reads_unaligned / exon_size
    per_base_not_primary_alignment = count_reads_not_primary_alignment / exon_size
    per_base_mate_unaligned = count_reads_mate_unaligned / exon_size
    per_base_mate_on_other_chr = count_reads_mate_on_other_chr / exon_size
    per_base_mate_gt_10kb_away = count_reads_mate_gt_10kb_away / exon_size
    per_base_soft_clipped = count_reads_soft_clipped / exon_size
    per_base_hard_clipped = count_reads_hard_clipped / exon_size

    return {
        'min_per_nucleotide_depth': min_per_nucleotide_depth,
        'mean_per_nucleotide_depth': mean_per_nucleotide_depth,
        'median_per_nucleotide_depth': median_per_nucleotide_depth,
        'max_per_nucleotide_depth': max_per_nucleotide_depth,

        'bq30_min_per_nucleotide_depth': bq30_min_per_nucleotide_depth,
        'bq30_mean_per_nucleotide_depth': bq30_mean_per_nucleotide_depth,
        'bq30_median_per_nucleotide_depth': bq30_median_per_nucleotide_depth,
        'bq30_max_per_nucleotide_depth': bq30_max_per_nucleotide_depth,

        'mq20_min_per_nucleotide_depth': mq20_min_per_nucleotide_depth,
        'mq20_mean_per_nucleotide_depth': mq20_mean_per_nucleotide_depth,
        'mq20_median_per_nucleotide_depth': mq20_median_per_nucleotide_depth,
        'mq20_max_per_nucleotide_depth': mq20_max_per_nucleotide_depth,

        'total_reads_processed': total_reads_processed,

        'count_reads_unaligned': count_reads_unaligned,
        'count_reads_not_primary_alignment': count_reads_not_primary_alignment,
        'count_reads_mate_unaligned': count_reads_mate_unaligned,
        'count_reads_mate_on_other_chr': count_reads_mate_on_other_chr,
        'count_reads_mate_gt_10kb_away': count_reads_mate_gt_10kb_away,
        'count_reads_soft_clipped': count_reads_soft_clipped,
        'count_reads_hard_clipped': count_reads_hard_clipped,
        'count_reads_fw_strand': count_reads_fw_strand,

        'pct_reads_unaligned': pct_reads_unaligned,
        'pct_reads_not_primary_alignment': pct_reads_not_primary_alignment,
        'pct_reads_mate_unaligned': pct_reads_mate_unaligned,
        'pct_reads_mate_on_other_chr': pct_reads_mate_on_other_chr,
        'pct_reads_mate_gt_10kb_away': pct_reads_mate_gt_10kb_away,
        'pct_reads_soft_clipped': pct_reads_soft_clipped,
        'pct_reads_hard_clipped': pct_reads_hard_clipped,
        'pct_reads_fw_strand': pct_reads_fw_strand,

        'per_base_unaligned': per_base_unaligned,
        'per_base_not_primary_alignment': per_base_not_primary_alignment,
        'per_base_mate_unaligned': per_base_mate_unaligned,
        'per_base_mate_on_other_chr': per_base_mate_on_other_chr,
        'per_base_mate_gt_10kb_away': per_base_mate_gt_10kb_away,
        'per_base_soft_clipped': per_base_soft_clipped,
        'per_base_hard_clipped': per_base_hard_clipped
    }

def write_gene_exon_coverage_info_to_file(genes_with_exons, gene_exon_coverage_info, output_root):
    """
    Writes gene exon coverage info to a file at output_root.gene_exon_coverage_info.txt
    
    Args:
        genes_with_exons: dict of gene > list of exons
        gene_exon_coverage_info: dict of gene > list of exon coverage info
        output_root: root for report files and plots
        
    Returns:
        None
    """
    with open(output_root + ".gene_exon_coverage_info.txt", 'w') as fh_gene_exon_coverage_info:
        fh_gene_exon_coverage_info.write("\t".join([
            "gene",
            "gene_coordinates",
            "exon_number",
            "chr",
            "start",
            "end",
            "min_per_nucleotide_depth",
            "mean_per_nucleotide_depth",
            "median_per_nucleotide_depth",
            "max_per_nucleotide_depth",
            "bq30_min_per_nucleotide_depth",
            "bq30_mean_per_nucleotide_depth",
            "bq30_median_per_nucleotide_depth",
            "bq30_max_per_nucleotide_depth",
            "mq20_min_per_nucleotide_depth",
            "mq20_mean_per_nucleotide_depth",
            "mq20_median_per_nucleotide_depth",
            "mq20_max_per_nucleotide_depth",
            "total_reads_processed",
            "count_reads_unaligned",
            "count_reads_not_primary_alignment",
            "count_reads_mate_unaligned",
            "count_reads_mate_on_other_chr",
            "count_reads_mate_gt_10kb_away",
            "count_reads_soft_clipped",
            "count_reads_hard_clipped",
            "count_reads_fw_strand",
            "pct_reads_unaligned",
            "pct_reads_not_primary_alignment",
            "pct_reads_mate_unaligned",
            "pct_reads_mate_on_other_chr",
            "pct_reads_mate_gt_10kb_away",
            "pct_reads_soft_clipped",
            "pct_reads_hard_clipped",
            "pct_reads_fw_strand",
            "per_base_unaligned",
            "per_base_not_primary_alignment",
            "per_base_mate_unaligned",
            "per_base_mate_on_other_chr",
            "per_base_mate_gt_10kb_away",
            "per_base_soft_clipped",
            "per_base_hard_clipped"
        ]) + "\n")

        for gene in genes_with_exons:
            for exon_idx, exon in enumerate(genes_with_exons[gene]):
                exon_info = gene_exon_coverage_info[gene][exon_idx]
                fh_gene_exon_coverage_info.write("\t".join([
                    gene[0],
                    gene[1],
                    str(exon_idx),
                    exon[0],
                    str(exon[1]),
                    str(exon[2]),
                    str(exon_info["min_per_nucleotide_depth"]),
                    str(exon_info["mean_per_nucleotide_depth"]),
                    str(exon_info["median_per_nucleotide_depth"]),
                    str(exon_info["max_per_nucleotide_depth"]),
                    str(exon_info["bq30_min_per_nucleotide_depth"]),
                    str(exon_info["bq30_mean_per_nucleotide_depth"]),
                    str(exon_info["bq30_median_per_nucleotide_depth"]),
                    str(exon_info["bq30_max_per_nucleotide_depth"]),
                    str(exon_info["mq20_min_per_nucleotide_depth"]),
                    str(exon_info["mq20_mean_per_nucleotide_depth"]),
                    str(exon_info["mq20_median_per_nucleotide_depth"]),
                    str(exon_info["mq20_max_per_nucleotide_depth"]),
                    str(exon_info["total_reads_processed"]),
                    str(exon_info["count_reads_unaligned"]),
                    str(exon_info["count_reads_not_primary_alignment"]),
                    str(exon_info["count_reads_mate_unaligned"]),
                    str(exon_info["count_reads_mate_on_other_chr"]),
                    str(exon_info["count_reads_mate_gt_10kb_away"]),
                    str(exon_info["count_reads_soft_clipped"]),
                    str(exon_info["count_reads_hard_clipped"]),
                    str(exon_info["count_reads_fw_strand"]),
                    str(exon_info["pct_reads_unaligned"]),
                    str(exon_info["pct_reads_not_primary_alignment"]),
                    str(exon_info["pct_reads_mate_unaligned"]),
                    str(exon_info["pct_reads_mate_on_other_chr"]),
                    str(exon_info["pct_reads_mate_gt_10kb_away"]),
                    str(exon_info["pct_reads_soft_clipped"]),
                    str(exon_info["pct_reads_hard_clipped"]),
                    str(exon_info["pct_reads_fw_strand"]),
                    str(exon_info["per_base_unaligned"]),
                    str(exon_info["per_base_not_primary_alignment"]),
                    str(exon_info["per_base_mate_unaligned"]),
                    str(exon_info["per_base_mate_on_other_chr"]),
                    str(exon_info["per_base_mate_gt_10kb_away"]),
                    str(exon_info["per_base_soft_clipped"]),
                    str(exon_info["per_base_hard_clipped"])
                ]) + "\n")


def plot_gene_exon_coverage_info(genes_with_exons, gene_info, output_root):
    """
    Plots gene exon coverage info based on calculated values. Produces a pdf file at output_root.gene_exon_coverage_info.pdf
    
    Args:
        genes_with_exons: dict of gene > list of exons
        gene_info: dict of gene > list of exon coverage info
        output_root: root for report files and plots
        
    Returns:
        None
    """
    print('Plotting gene exon coverage info')

    max_gene_exon = 0
    for gene in genes_with_exons:
        if len(genes_with_exons[gene]) > max_gene_exon:
            max_gene_exon = len(genes_with_exons[gene])

    values_to_plot_with_descriptions = [('min_per_nucleotide_depth','Minimum per-nucleotide read depth per exon'), 
                      ('mean_per_nucleotide_depth','Mean per-nucleotide read depth per exon'),
                      ('median_per_nucleotide_depth','Median per-nucleotide read depth per exon'),
                      ('max_per_nucleotide_depth','Maximum per-nucleotide read depth per exon'), 
                      ('bq30_min_per_nucleotide_depth','Minimum per-nucleotide read depth for bases with quality >30 per exon'),
                      ('bq30_mean_per_nucleotide_depth','Mean per-nucleotide read depth for bases with quality >30 per exon'),
                      ('bq30_median_per_nucleotide_depth','Median per-nucleotide read depth for bases with quality >30 per exon'),
                      ('bq30_max_per_nucleotide_depth','Maximum per-nucleotide read depth for bases with quality >30 per exon'), 
                      ('mq20_min_per_nucleotide_depth','Minimum per-nucleotide read depth for reads with mapping quality >20 per exon'),
                      ('mq20_mean_per_nucleotide_depth','Mean per-nucleotide read depth for reads with mapping quality >20 per exon'),
                      ('mq20_median_per_nucleotide_depth','Median per-nucleotide read depth for reads with mapping quality >20 per exon'),
                      ('mq20_max_per_nucleotide_depth','Maximum per-nucleotide read depth for reads with mapping quality >20 per exon'), 
                      ('total_reads_processed','Total reads per exon'),
                      ('count_reads_unaligned','Count of reads unaligned per exon'),
                      ('count_reads_not_primary_alignment','Count of reads not primary alignment per exon'),
                      ('count_reads_mate_unaligned','Count of reads with mate unaligned per exon'),
                      ('count_reads_mate_on_other_chr','Count of reads with mate on other chromosome per exon'), 
                      ('count_reads_mate_gt_10kb_away','Count of reads with mate >10kb away per exon'),
                      ('count_reads_soft_clipped','Count of soft clipped reads per exon'),
                      ('count_reads_hard_clipped','Count of hard clipped reads per exon'),
                      ('count_reads_fw_strand','Count of reads on forward strand per exon'), 
                      ('pct_reads_unaligned','Percent of unaligned reads per exon'),
                      ('pct_reads_not_primary_alignment','Percent of not primary alignment reads per exon'),
                      ('pct_reads_mate_unaligned','Percent of reads with mate unaligned per exon'),
                      ('pct_reads_mate_on_other_chr','Percent of reads with mate on another chromosome per exon'), 
                      ('pct_reads_mate_gt_10kb_away','Percent of reads with mate >10kb away per exon'),
                      ('pct_reads_soft_clipped','Percent of soft clipped reads per exon'),
                      ('pct_reads_hard_clipped','Percent of hard clipped reads per exon'),
                      ('pct_reads_fw_strand','Percent of reads on forward strand per exon'), 
                      ('per_base_unaligned','Number of unaligned reads per base per exon'),
                      ('per_base_not_primary_alignment','Number of not primary alignment reads per base per exon'),
                      ('per_base_mate_unaligned','Number of reads with mate unaligned per base per exon'),
                      ('per_base_mate_on_other_chr','Number of reads with mate on another chromosome per base per exon'),
                      ('per_base_mate_gt_10kb_away','Number of reads with mate >10kb away per base per exon'),
                      ('per_base_soft_clipped','Number of soft clipped reads per base per exon'),
                      ('per_base_hard_clipped','Number of hard clipped reads per base per exon')
                      ]

    output_pdf = output_root + ".gene_exon_coverage_info.pdf"
    with PdfPages(output_pdf) as pdf:
        for value,value_description in values_to_plot_with_descriptions:
            big_df_arr = []
            for gene in genes_with_exons:
                gene_val_arr = [np.nan]*max_gene_exon
                for exon_idx, exon in enumerate(genes_with_exons[gene]):
                    gene_val_arr[exon_idx] = gene_info[gene][exon_idx][value]
                big_df_arr.append(gene_val_arr)

            big_df = pd.DataFrame(big_df_arr, columns=['Exon' + str(x) for x in range(1,max_gene_exon+1)])
            big_df.index = [x[0] + ' ' + x[1] for x in genes_with_exons]

            plt.figure(figsize=(8, 6))
            ax = sns.heatmap(big_df, cmap='coolwarm', mask=big_df.isnull())
            ax.set_title(value_description)
            plt.tight_layout()
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
    print('Plotted to ' + output_pdf)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some files.')
    parser.add_argument('--output_root', type=str, required=True, help='Root for report files and plots')
    parser.add_argument('--bam_input', type=str, required=True, help='Path to the BAM input file')
    parser.add_argument('--ref_gene_file', type=str, required=True, help='Path to the refGene.txt file')
    parser.add_argument('--gene_exon_list_file', type=str, required=True, help='Path to the gene exon list file')

    parser.add_argument('chr_index', type=int, nargs='?', default=2, help='The index of the chr column in the refGene.txt file')
    parser.add_argument('strand_index', type=int, nargs='?', default=3, help='The index of the strand column in the refGene.txt file')
    parser.add_argument('tx_start_index', type=int, nargs='?', default=4, help='The index of the txStart column in the refGene.txt file')
    parser.add_argument('tx_end_index', type=int, nargs='?', default=5, help='The index of the txEnd column in the refGene.txt file')
    parser.add_argument('cds_start_index', type=int, nargs='?', default=6, help='The index of the cdsStart column in the refGene.txt file')
    parser.add_argument('cds_end_index', type=int, nargs='?', default=7, help='The index of the cdsEnd column in the refGene.txt file')
    parser.add_argument('exon_index', type=int, nargs='?', default=9, help='The index of the exon column in the refGene.txt file')
    parser.add_argument('intron_index', type=int, nargs='?', default=10, help='The index of the intron column in the refGene.txt file')
    parser.add_argument('gene_id_index', type=int, nargs='?', default=12, help='The gene id index column in the refGene.txt file')
    parser.add_argument('exon_status_index', type=int, nargs='?', default=15, help='The index of the exon_status column in the refGene.txt file')
    args = parser.parse_args()

    try:
        output = subprocess.check_output(["samtools", "--version"], stderr=subprocess.STDOUT)
        print("Samtools is installed and functional.")
    except subprocess.CalledProcessError as e:
        print("Samtools is not installed or not functional.")
        print("Error details: ", e.output)

    gene_list = get_genes_from_file(args.gene_exon_list_file)

    genes_with_exons = get_exons_from_refgene(gene_list, args.ref_gene_file, args.chr_index, args.strand_index, args.tx_start_index, args.tx_end_index, args.cds_start_index, args.cds_end_index, args.exon_index, args.intron_index, args.gene_id_index, args.exon_status_index)

    gene_info = {}
    for gene_idx,gene in enumerate(genes_with_exons):
        print('Processing gene ' + gene[0] + ' (' + str(gene_idx+1) + ' of ' + str(len(genes_with_exons)) + ' with ' + str(len(genes_with_exons[gene])) + ' exons)')
        gene_exon_coverage_info = []
        for exon in genes_with_exons[gene]:
            exon_coverage = get_exon_coverage(args.bam_input, exon)
            gene_exon_coverage_info.append(exon_coverage)
        gene_info[gene] = gene_exon_coverage_info


    write_gene_exon_coverage_info_to_file(genes_with_exons, gene_info, args.output_root)
    plot_gene_exon_coverage_info(genes_with_exons, gene_info, args.output_root)

    print('Finished processing ' + str(len(gene_list)) + ' genes.')
