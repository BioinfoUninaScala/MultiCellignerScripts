# MultiCellignerScripts

R code for generating and analyzing single-omics and multiomics alignments between tumor samples and cancer cell line models.

------------------------------------------------------------------------
# Data modalities

The framework supports:

- mRNA expression (RNA-seq)
- DNA methylation (HM450K and RRBS)
- Mutational signatures (WES)

and their combinations:

- DNA methylation + mRNA expression
- DNA methylation + Mutational signatures
- mRNA expression + Mutational signatures
- DNA methylation + mRNA expression + Mutational signatures

# Repository structure

```text
DATA/
├── Annotation/
│   ├── ann_multiomics_v9.rds
│   ├── fullsample_ann_expression.csv
│   └── sample_info_CCLE.csv
│
├── Expression/
│   ├── hgnc_complete_set
│   ├── CCLE_expression_full.csv
│   ├── combined_mat.rds
│   └── TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv 
│
├── Methylation/
│   ├── CCLE_RRBS_TSS_1kb_20180614.txt 
│   ├── combined_mat_meth.rds
│   ├── GDC-PANCAN_meth450.tsv.gz
│   └── id_tcga_pancan
│
└─── Mutational_process/
    ├── CCLE_mutations.csv
    ├── combined_mat_mut.rds
    └── mut_pancancer.rds

Large matrices are not distributed with the repository and must be downloaded separately.

GDC-PANCAN_meth450.tsv.gz
https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz

CCLE_RRBS_TSS_1kb_20180614.txt 
https://depmap.org/portal/data_page/?tab=allData

TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv 
https://xenabrowser.net/datapages/?dataset=TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv&host=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

CCLE_expression_full.csv
https://depmap.org/portal/data_page/?tab=allData

CCLE_mutations.csv
https://depmap.org/portal/data_page/?tab=allData
