# MultiCellignerScripts

R code for generating and analyzing single-omics and multiomics alignments between tumor samples and cancer cell line models.

------------------------------------------------------------------------

<p align="center">

<img src="inst/www/MultiCelligner_Logo_2.png" width="200"/>

</p>

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

------------------------------------------------------------------------

<p align="center">

<img src="inst/www/graphicalAbstract.png" width="1000"/>

</p>

------------------------------------------------------------------------

# Repository structure

```text
DATA/
├── Annotation/
├── Expression/
├── Methylation/
└── Mutational_process/
DATA/
├── Annotation/
├── Expression/
├── Methylation/
└── Mutational_process/

Large matrices are not distributed with the repository and must be downloaded separately.


