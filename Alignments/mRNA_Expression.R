
#### Data availability ####
# Large input matrices are not distributed with the repository.
# Please download the required files and place them in the corresponding folders:
#
# DATA/
# ├── Expression/
# │   ├── hgnc_complete_set
# │   ├── TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv
# │   └── CCLE_expression_full.csv
# └── Annotation/
#     └── fullsample_ann_expression.csv



library(tidyverse)
library(Seurat)

source(utils)
source(global_parameters)

###### Get combined_mat of mRNA Expression data ###### 


hgnc_file <- file.path(
  "DATA", "Expression", "hgnc_complete_set"
)

if (!file.exists(hgnc_file)) {
  stop(
    "Missing file: hgnc_complete_set\n",
    "Please place the file in DATA/Expression/"
  )
}

hgnc.complete.set <- read.delim(hgnc_file)



annotation_file <- file.path(
  "DATA", "Annotation", "fullsample_ann_expression.csv"
)

if (!file.exists(annotation_file)) {
  stop(
    "Missing file: fullsample_ann_expression.csv\n",
    "Please place the file in DATA/Annotation/"
  )
}

fullsample_ann_expression <- read.csv(annotation_file)



tcga_expression_file <- file.path(
  "DATA", "Expression",
  "TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv"
)

if (!file.exists(tcga_expression_file)) {
  stop(
    "Missing file: TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv\n",
    "Please place the file in DATA/Expression/"
  )
}

TCGA_mat <- read.delim(tcga_expression_file)



ccle_expression_file <- file.path(
  "DATA", "Expression", "CCLE_expression_full.csv"
)

if (!file.exists(ccle_expression_file)) {
  stop(
    "Missing file: CCLE_expression_full.csv\n",
    "Please place the file in DATA/Expression/"
  )
}

CCLE_mat <- read_csv(ccle_expression_file)

###### Get combined_mat of mRNA Expression data ###### 

TCGA_mat_1 <- TCGA_mat %>%
  as.data.frame() %>%
  tibble::column_to_rownames('Gene') %>%
  as.matrix() %>%
  t()

common_genes <- intersect(colnames(TCGA_mat_1), hgnc$symbol)
TCGA_mat_1 <- TCGA_mat_1[,common_genes]
hgnc.complete.set <- filter(hgnc, symbol %in% common_genes)
hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$symbol),]
rownames(hgnc.complete.set) <- hgnc.complete.set$symbol
hgnc.complete.set <- hgnc.complete.set[common_genes,]
colnames(TCGA_mat_1) <- hgnc.complete.set$ensembl_gene_id

a <- CCLE_mat[,1]
a <- as.vector(a)
a <- a[["...1"]]

CCLE_mat_1 <- as.matrix(CCLE_mat)
rownames(CCLE_mat_1) <- a
CCLE_mat_1 <- CCLE_mat_1[,-1]

colnames(CCLE_mat_1) <- stringr::str_match(colnames(CCLE_mat_1), '\\((.+)\\)')[,2]

func_genes <- dplyr::filter(hgnc.complete.set, !locus_group %in% c('non-coding RNA', 'pseudogene'))$ensembl_gene_id
genes_used <- intersect(colnames(TCGA_mat_1), colnames(CCLE_mat_1))
genes_used <- intersect(genes_used, func_genes)

TCGA_mat_1 <- TCGA_mat_1[,genes_used]
CCLE_mat_1 <- CCLE_mat_1[,genes_used]

sum(!rownames(TCGA_mat_1) %in% fullsample_ann_expression$sampleID)

b <- gsub(pattern = "\\.", replacement = "-", x = rownames(TCGA_mat_1))
rownames(TCGA_mat_1) <- b

TCGA_ann <- dplyr::filter(fullsample_ann_expression, type=='tumor')
CCLE_ann <- dplyr::filter(fullsample_ann_expression, type=='CL')

CCLE_mat_2 <- as.matrix(apply(CCLE_mat_1, 2, as.double))
rownames(CCLE_mat_2) <- rownames(CCLE_mat_1)

common_genes <- intersect(colnames(TCGA_mat_1), colnames(CCLE_mat_2))

hgnc.complete.set <- hgnc.complete.set %>%
  dplyr::select(Gene = ensembl_gene_id, Symbol = symbol) %>%
  filter(Gene %in% common_genes)
hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$Gene), ]
rownames(hgnc.complete.set) <- hgnc.complete.set$Gene
hgnc.complete.set <- hgnc.complete.set[common_genes, ]

gene_stats <- data.frame(
  Tumor_SD = apply(TCGA_mat_1, 2, sd, na.rm=T),
  CCLE_SD = apply(CCLE_mat_2, 2, sd, na.rm=T),
  Tumor_mean = colMeans(TCGA_mat_1, na.rm=T),
  CCLE_mean = colMeans(CCLE_mat_2, na.rm=T),
  Gene = common_genes,
  stringsAsFactors = F) %>% 
  dplyr::mutate(max_SD = pmax(Tumor_SD, CCLE_SD, na.rm=T)) #add avg and max SD per gene

gene_stats <- left_join(hgnc.complete.set, gene_stats, by = "Gene")

comb_ann <- rbind(
  TCGA_ann %>% dplyr::select(sampleID, lineage, subtype) %>%
    dplyr::mutate(type = "tumor"),
  CCLE_ann %>% dplyr::select(sampleID, lineage, subtype) %>%
    dplyr::mutate(type = "CL")
)

TCGA_obj <- create_Seurat_object(TCGA_mat_1, TCGA_ann, type = "tumor")
CCLE_obj <- create_Seurat_object(CCLE_mat_2, CCLE_ann, type = "CL")

TCGA_obj <- cluster_data(TCGA_obj)
CCLE_obj <- cluster_data(CCLE_obj)

tumor_DE_genes <- find_differentially_expressed_genes(TCGA_obj)
CL_DE_genes <- find_differentially_expressed_genes(CCLE_obj)

DE_genes <- full_join(tumor_DE_genes, CL_DE_genes, by = "Gene", suffix = c("_tumor", "_CL")) %>%
  mutate(
    tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
    CL_rank = dplyr::dense_rank(-gene_stat_CL),
    best_rank = pmin(tumor_rank, CL_rank, na.rm = T)
  ) %>%
  dplyr::left_join(gene_stats, by = "Gene")

### take genes that are ranked in the top 1000 from either dataset, used for finding mutual nearest neighbors
DE_gene_set <- DE_genes %>%
  dplyr::filter(best_rank < MultiCelligner_parameters$top_DE_genes_per) %>%
  .[["Gene"]]

### calculating all cPCs in order to empower reproducibility

cov_diff_eig <- run_cPCA(TCGA_obj, CCLE_obj, pc_dims = NULL)
cur_vecs <- cov_diff_eig$rotation[, MultiCelligner_parameters$remove_cPCA_dims, drop = FALSE]
rownames(cur_vecs) <- colnames(TCGA_mat_1)

TCGA_cor <- resid(lm(t(TCGA_mat_1) ~ 0 + cur_vecs)) %>% t()
CCLE_cor <- resid(lm(t(CCLE_mat_2) ~ 0 + cur_vecs)) %>% t()

mnn_res <- run_MNN(CCLE_cor, TCGA_cor,
                   k1 = MultiCelligner_parameters$mnn_k_tumor, k2 = MultiCelligner_parameters$mnn_k_CL, ndist = MultiCelligner_parameters$mnn_ndist,
                   subset_genes = DE_gene_set
)

combined_mat <- rbind(mnn_res$corrected, CCLE_cor)

saveRDS(combined_mat, "DATA/Expression/combined_mat.rds")

