# Genomic Region Coverage

## Requirements
This script requires:
 - samtools
 - python libraries:
   - matplotlib
   - seaborn

## Usage
```
python get_exon_coverage.py [-h] --output_root OUTPUT_ROOT --bam_input BAM_INPUT --ref_gene_file REF_GENE_FILE
                            --gene_exon_list_file GENE_EXON_LIST_FILE
```

- `--output_root OUTPUT_ROOT`
                        Root for output report files and plots
  
- `--bam_input BAM_INPUT`
                        Path to the BAM input file. BAI file is assumed to exist and be in the standard location.
  
- `--ref_gene_file REF_GENE_FILE`
                        Path to the refGene.txt file containing gene and exon information, downloadable from here: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
  
- `--gene_exon_list_file GENE_EXON_LIST_FILE`
                        Path to the gene exon list file. This is a file with each gene of interest on each line.

## Output
Statistics are reported for the following metrics:

- `min_per_nucleotide_depth`: Minimum per-nucleotide read depth per exon
- `mean_per_nucleotide_depth`: Mean per-nucleotide read depth per exon
- `median_per_nucleotide_depth`: Median per-nucleotide read depth per exon
- `max_per_nucleotide_depth`: Maximum per-nucleotide read depth per exon
- `bq30_min_per_nucleotide_depth`: Minimum per-nucleotide read depth for bases with quality >30 per exon
- `bq30_mean_per_nucleotide_depth`: Mean per-nucleotide read depth for bases with quality >30 per exon
- `bq30_median_per_nucleotide_depth`: Median per-nucleotide read depth for bases with quality >30 per exon
- `bq30_max_per_nucleotide_depth`: Maximum per-nucleotide read depth for bases with quality >30 per exon
- `mq20_min_per_nucleotide_depth`: Minimum per-nucleotide read depth for reads with mapping quality >20 per exon
- `mq20_mean_per_nucleotide_depth`: Mean per-nucleotide read depth for reads with mapping quality >20 per exon
- `mq20_median_per_nucleotide_depth`: Median per-nucleotide read depth for reads with mapping quality >20 per exon
- `mq20_max_per_nucleotide_depth`: Maximum per-nucleotide read depth for reads with mapping quality >20 per exon
- `total_reads_processed`: Total reads per exon
- `count_reads_unaligned`: Count of reads unaligned per exon
- `count_reads_not_primary_alignment`: Count of reads not primary alignment per exon
- `count_reads_mate_unaligned`: Count of reads with mate unaligned per exon
- `count_reads_mate_on_other_chr`: Count of reads with mate on other chromosome per exon
- `count_reads_mate_gt_10kb_away`: Count of reads with mate >10kb away per exon
- `count_reads_soft_clipped`: Count of soft clipped reads per exon
- `count_reads_hard_clipped`: Count of hard clipped reads per exon
- `count_reads_fw_strand`: Count of reads on forward strand per exon
- `pct_reads_unaligned`: Percent of unaligned reads per exon
- `pct_reads_not_primary_alignment`: Percent of not primary alignment reads per exon
- `pct_reads_mate_unaligned`: Percent of reads with mate unaligned per exon
- `pct_reads_mate_on_other_chr`: Percent of reads with mate on another chromosome per exon
- `pct_reads_mate_gt_10kb_away`: Percent of reads with mate >10kb away per exon
- `pct_reads_soft_clipped`: Percent of soft clipped reads per exon
- `pct_reads_hard_clipped`: Percent of hard clipped reads per exon
- `pct_reads_fw_strand`: Percent of reads on forward strand per exon
- `per_base_unaligned`: Number of unaligned reads per base per exon
- `per_base_not_primary_alignment`: Number of not primary alignment reads per base per exon
- `per_base_mate_unaligned`: Number of reads with mate unaligned per base per exon
- `per_base_mate_on_other_chr`: Number of reads with mate on another chromosome per base per exon
- `per_base_mate_gt_10kb_away`: Number of reads with mate >10kb away per base per exon
- `per_base_soft_clipped`: Number of soft clipped reads per base per exon
- `per_base_hard_clipped`: Number of hard clipped reads per base per exon
