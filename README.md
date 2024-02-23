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

