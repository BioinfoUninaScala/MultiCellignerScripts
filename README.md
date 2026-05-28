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
```text
DATA/
├── Annotation/
│   ├── ann_multiomics_v9.rds
│   ├── fullsample_ann_expression.csv
│   └── sample_info_CCLE.csv
│
├── Expression/
│   ├── 22258446
│   ├── 25712759
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


