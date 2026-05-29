
#### Data availability ####
# Large input matrices are not distributed with the repository.
# Please download the required files and place them in the corresponding folders:
#
# DATA/
# ├── Methylation/
# │   ├── CCLE_RRBS_TSS_1kb_20180614.txt
# │   ├── GDC-PANCAN_meth450.tsv.gz
# │   └── id_tcga_pancan
# └── Annotation/
#     ├── sample_info_CCLE.csv
#     ├── final_TCGAprojects_ExpMethAnnot.RDS
#     └── ann_multiomics_v9.rds

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GenomicRanges)
library(tidyverse)
library(impute)
library(magrittr)
library(Seurat)
library(data.table)
library(stringr)

source(utils)
source(global_parameters)

#### Pre-process TCGA Data ####
# load CCLE gene promoter coordinates

ccle_rrbs_file <- file.path(
  "DATA", "Methylation", "CCLE_RRBS_TSS_1kb_20180614.txt"
)

if (!file.exists(ccle_rrbs_file)) {
  stop(
    "Missing file: CCLE_RRBS_TSS_1kb_20180614.txt\n",
    "Please place the file in DATA/Methylation/"
  )
}

CCLE_RRBS_1kb <- read_delim(
  ccle_rrbs_file,
  delim = " ",
  escape_double = FALSE,
  trim_ws = TRUE
)



##### Load TCGA Data #####

tcga_meth_file <- file.path(
  "DATA", "Methylation", "GDC-PANCAN_meth450.tsv.gz"
)

if (!file.exists(tcga_meth_file)) {
  stop(
    "Missing file: GDC-PANCAN_meth450.tsv.gz\n",
    "Please place the file in DATA/Methylation/"
  )
}

pancancer_meth <- fread(
  tcga_meth_file,
  header = FALSE
)


id_tcga_file <- file.path(
  "DATA", "Methylation", "id_tcga_pancan"
)

if (!file.exists(id_tcga_file)) {
  stop(
    "Missing file: id_tcga_pancan\n",
    "Please place the file in DATA/Methylation/"
  )
}

id_tcga_pancan <- read_table(
  id_tcga_file,
  col_names = FALSE
)

##### Load CCLE annotation #####

ccle_annot_file <- file.path(
  "DATA", "Annotation", "sample_info_CCLE.csv"
)

if (!file.exists(ccle_annot_file)) {
  stop(
    "Missing file: sample_info_CCLE.csv\n",
    "Please place the file in DATA/Annotation/"
  )
}

sample_info_CCLE <- read_csv(
  ccle_annot_file,
  col_types = cols(depmap_public_comments = col_character())
)



##### Load multiomics annotation #####

ann_multiomics_file <- file.path(
  "DATA", "Annotation", "ann_multiomics_v9.rds"
)

if (!file.exists(ann_multiomics_file)) {
  stop(
    "Missing file: ann_multiomics_v9.rds\n",
    "Please place the file in DATA/Annotation/"
  )
}

ann_multiomics_v9 <- readRDS(ann_multiomics_file)

##### Load TCGA projects annotation #####

tcga_projects_file <- file.path(
  "DATA", "Annotation", "final_TCGAprojects_ExpMethAnnot.RDS"
)

if (!file.exists(tcga_projects_file)) {
  stop(
    "Missing file: final_TCGAprojects_ExpMethAnnot.RDS\n",
    "Please place the file in DATA/Annotation/"
  )
}

sample_annot <- readRDS(tcga_projects_file)


#### Pre-process TCGA Data ####
# load CCLE gene promoter coordinates
# load cpg coordinates

loc <- data("Locations")
loc <- as.data.frame(Locations)
loc <- loc %>% mutate('cpg_id' = rownames(loc))
loc_range <- GRanges(seqnames = loc$chr, ranges = IRanges(start = loc$pos, end = loc$pos), strand = loc$strand)

rrbs_loc <- CCLE_RRBS_1kb %>%
  select("gene", "chr", "start" = fpos, "stop" = tpos, "strand")
rrbs_loc_range <- as(rrbs_loc, "GRanges")

over <- findOverlaps(loc_range, rrbs_loc_range)

ann_cpg <- cbind(loc[queryHits(over),], rrbs_loc[subjectHits(over),]) %>%
  select("cpg_id","gene") %>% distinct %>% as_tibble

ann_cpg %>% nrow #52,119

name_pancancer <- c("cpg_id", id_tcga_pancan$X1)
colnames(pancancer_meth) <- name_pancancer
pancancer_meth[1:5, 1:5]

##### Remove cpg with no signal #####
cpg_to_remove <- rowSums(is.na(pancancer_meth)) == ncol(pancancer_meth)-1
sum(cpg_to_remove) #5363

new_pancancer_meth <- pancancer_meth[!cpg_to_remove, ]
nrow(pancancer_meth) #51106
nrow(new_pancancer_meth) #45743

##### Annotate TCGA cpg with corresponding genes #####
all(new_pancancer_meth$cpg_id %in% ann_cpg$cpg_id) #TRUE
pancancer_meth_ann <- right_join(ann_cpg, new_pancancer_meth, by = "cpg_id", multiple = "all")
dim(pancancer_meth_ann) # 46644  9738
pancancer_meth_ann[1:5, 1:5]
pancancer_meth_ann %>%  group_by(gene) %>% dplyr::summarise(n = n() ) %>% arrange(desc(n))
# A tibble: 13,934 × 2
# gene          n
# <chr>     <int>
#   1 GNAS         51
# 2 GATA3-AS1    27
# 3 IGF2         27
# 4 BLCAP        25
# 5 BRD2         23

f_pancancer_meth_ann <- pancancer_meth_ann
dim(f_pancancer_meth_ann) #46644  9738
f_pancancer_meth_ann[1:5,1:5]

length(unique(f_pancancer_meth_ann$gene)) #13934

##### From TCGA cpg matrix to TCGA gene promoter matrix #####

###### 1) summarize at gene level ######
cpg_meth_mat <- as.matrix(f_pancancer_meth_ann[, -c(1,2)])
gene_vec <- f_pancancer_meth_ann$gene
split_idx <- split(seq_along(gene_vec), gene_vec)

meth_means_mat <- sapply(split_idx, function(idx) {
  colMeans(cpg_meth_mat[idx, , drop = FALSE], na.rm = TRUE)
})
dim(meth_means_mat) #9736  13934
meth_means_mat[1:5,1:5]

t_meth_means_mat <- t(meth_means_mat)
t_meth_means_mat[1:5,1:5]


###### 2) filtering by sample_type + remove duplicated ######
sample_names <- colnames(t_meth_means_mat)

sample_annot <- readRDS(paste0('DATA/Annotation/final_TCGAprojects_ExpMethAnnot.RDS'))
dim(sample_annot) #11895     9
all(colnames(pancancer_meth_ann)[-c(1:2)] %in% sample_names)
sample_annot$sample_type %>% table

# sample_type_to_remove <- c('Metastatic', 'Additional Metastatic', 'Solid Tissue Normal', 'Recurrent Tumor', 'Additional - New Primary')
sample_type_to_remove <- 'Solid Tissue Normal'
meth_sample_annot <- sample_annot %>%
  filter(!sample_type %in% sample_type_to_remove) %>%
  filter(sample.submitter_id %in% sample_names)
nrow(meth_sample_annot) # 8990

length(meth_sample_annot %>% pull(cases.submitter_id) %>% unique) # 8868
length(meth_sample_annot %>% pull(sampleID) %>% unique) # 8942

(duplicated_sample_df <- meth_sample_annot %>% group_by(sampleID) %>% summarize(n=n()) %>% filter(n > 1) %>% arrange(desc(n)))
meth_sample_annot %>% filter(sampleID %in% duplicated_sample_df$sampleID) %>% pull(sample_type) %>% table
meth_sample_annot %>% filter(sampleID %in% duplicated_sample_df$sampleID) %>% group_by(sampleID) %>% summarise(n = n()) %>% arrange(desc(n))
meth_sample_annot %>% filter(sampleID == 'TCGA-44-2656-01')
meth_sample_annot %>% filter(sampleID %in% duplicated_sample_df$sampleID) %>% group_by(sampleID) %>% summarise(n = n()) %>% pull(n) %>% unique

meth_sample_annot %>%
  mutate(last_letter = substr(sample.submitter_id, nchar(sample_names), nchar(sample_names) )) %>%
  pull(last_letter) %>% table
meth_sample_annot %>% filter(sampleID %in% duplicated_sample_df$sampleID) %>%
  mutate(last_letter = substr(sample.submitter_id, nchar(sample_names), nchar(sample_names) )) %>%
  pull(last_letter) %>% table

no_duplicated_sampleID <- meth_sample_annot %>% group_by(sampleID) %>% summarize(n=n()) %>% filter(n == 1) %>% pull(sampleID)
no_duplicated <- meth_sample_annot %>% filter(sampleID %in% no_duplicated_sampleID) %>% pull(sample.submitter_id)
to_keep_from_duplicated <- meth_sample_annot %>% filter(sampleID %in% duplicated_sample_df$sampleID) %>%
  mutate(last_letter = substr(sample.submitter_id, nchar(sample_names), nchar(sample_names) )) %>% filter(last_letter == 'A') %>% pull(sample.submitter_id)

sample_to_keep <- c(no_duplicated, to_keep_from_duplicated)
length(sample_to_keep) #8942

ff_pancancer_meth_ann <- t_meth_means_mat[, sample_to_keep]
dim(ff_pancancer_meth_ann) #13934  8942
ff_pancancer_meth_ann[1:5,1:5]

final_meth_sample_annot <- meth_sample_annot %>% filter(sample.submitter_id %in% colnames(ff_pancancer_meth_ann))
dim(final_meth_sample_annot)

###### 3) Transpose matrix ######

fff_pancancer_meth_ann <- t(ff_pancancer_meth_ann)
ff_pancancer_meth_ann[1:5,1:5]
fff_pancancer_meth_ann[1:5,1:5]

dim(fff_pancancer_meth_ann) #8942 13934


###### 4) remove NA > th samples ######
TCGAsample_th <- 0.3*ncol(fff_pancancer_meth_ann) #percentuale sui geni
filtSample_methTCGA_mat <- fff_pancancer_meth_ann[rowSums(is.na(fff_pancancer_meth_ann)) < TCGAsample_th,]
dim(filtSample_methTCGA_mat) #8942 13934
filtSample_methTCGA_mat[1:5,1:5]


###### 5)  remove NA > th genes ######
filtSample_methTCGA_mat <- filtSample_methTCGA_mat[,colSums(is.na(filtSample_methTCGA_mat)) < (0.5*nrow(filtSample_methTCGA_mat))]
dim(filtSample_methTCGA_mat) #  8942 13934


###### 6) remove lineages ######
lin_to_remove <- c("adrenal", "cervix", "germ cell", "thymus",  "bile duct", 'nerve', 'eye', 'soft tissue')
final_meth_sample_annot %>% pull(lineage) %>% sort %>% table
final_meth_sample_annot %>% filter(lineage %in% lin_to_remove) %>% pull(lineage ) %>% table
final_meth_sample_annot %>% filter(!lineage %in% lin_to_remove) %>% pull(lineage ) %>% table
final_meth_sample_annot %>% filter(!lineage %in% lin_to_remove) %>% pull(lineage ) %>% unique() %>% length #18

lin_TCGAsample_to_remove <- final_meth_sample_annot %>% filter(lineage %in% lin_to_remove) %>%
  pull(sample.submitter_id )
length(lin_TCGAsample_to_remove) #1234
f_filtSample_methTCGA_mat <- filtSample_methTCGA_mat[!rownames(filtSample_methTCGA_mat) %in% lin_TCGAsample_to_remove, ]
dim(f_filtSample_methTCGA_mat) #7708 13934


###### 7) Impute NA by lineage ######

final_tcga_annot <- final_meth_sample_annot %>%
  filter(sample.submitter_id %in% rownames(f_filtSample_methTCGA_mat))
nrow(final_tcga_annot)
final_tcga_annot %>% pull(lineage ) %>% sort %>% table
final_tcga_annot %>% pull(lineage ) %>% unique %>% length #18
final_tcga_annot %>% pull(sample_type) %>% unique


pre_TCGA_meth_impute <- imputeNA_by_lineage(mat = f_filtSample_methTCGA_mat,
                                        annot = final_tcga_annot, sample_column = 'sample.submitter_id')
dim(pre_TCGA_meth_impute) #13934  7708

sum(is.na(f_filtSample_methTCGA_mat)) #19842
sum(is.na(pre_TCGA_meth_impute)) #0

sample_info_CCLE$lineage <- tolower(gsub("_", " ", sample_info_CCLE$lineage))
unique(sample_info_CCLE$lineage) %>% sort

sample_info_CCLE %>% pull(primary_or_metastasis) %>% sort %>% table
sample_info_CCLE %>% pull(lineage) %>% sort %>% table

###### 1) summarize at gene level ######

nrow(CCLE_RRBS_1kb) == length(unique(CCLE_RRBS_1kb$gene)) #FALSE
CCLE_RRBS_1kb %>% group_by(gene) %>% dplyr::summarise(n = n()) %>% arrange(desc(n))
CCLE_meth <- CCLE_RRBS_1kb %>% select(-c('TSS_id', 'chr', 'fpos', 'tpos', 'strand', 'avg_coverage'))

CCLE_meth[1:5,1:5]

CCLE_meth_mat <- as.matrix(CCLE_meth[, -1])
gene_vec <- CCLE_meth$gene
split_idx <- split(seq_along(gene_vec), gene_vec)
CCLE_meth_means_mat <- sapply(split_idx, function(idx) {
  colMeans(CCLE_meth_mat[idx, , drop = FALSE], na.rm = TRUE)
})
dim(CCLE_meth_means_mat) #843 16494
CCLE_meth_means_mat[1:5, 1:5]

###### 2) rename ccle with depmap IDs ######
rownames(CCLE_meth_means_mat)[which(rownames(CCLE_meth_means_mat) == 'TT_OESOPHAGUS')] <- 'TDOTT_OESOPHAGUS'
rownames(CCLE_meth_means_mat) <- sub("\\_.*", "",rownames(CCLE_meth_means_mat))
CCLE_meth_means_mat[1:5, 1:5]

all(rownames(CCLE_meth_means_mat) %in% sample_info_CCLE$stripped_cell_line_name)
rownames(CCLE_meth_means_mat)[!rownames(CCLE_meth_means_mat) %in% sample_info_CCLE$stripped_cell_line_name] # "D341MED"
rownames(CCLE_meth_means_mat)[which(rownames(CCLE_meth_means_mat) == "D341MED")] <- "D341"
all(rownames(CCLE_meth_means_mat) %in% sample_info_CCLE$stripped_cell_line_name)

length(rownames(CCLE_meth_means_mat)) == length(unique(rownames(CCLE_meth_means_mat)))

sum(is.na(CCLE_meth_means_mat)) #523,840
sum(is.na(CCLE_meth_means_mat))/(ncol(CCLE_meth_means_mat) * nrow(CCLE_meth_means_mat))*100


na_remove_CCLE <- apply(CCLE_meth_means_mat, 1, function (x) all(is.na(x)))
sum(na_remove_CCLE) #0


renamed_CCLE_meth_means_mat <- CCLE_meth_means_mat %>% as.data.frame %>% rownames_to_column('stripped_cell_line_name' ) %>%
  left_join(., sample_info_CCLE %>% select('stripped_cell_line_name', 'DepMap_ID'), by ='stripped_cell_line_name') %>%
  select(-stripped_cell_line_name) %>% select(DepMap_ID, everything()) %>%  column_to_rownames('DepMap_ID')
dim(renamed_CCLE_meth_means_mat) #843 16494
renamed_CCLE_meth_means_mat[1:5, 1:5]


###### 3) remove ccle with many NAs ######

CCLEsample_th <- 0.1*ncol(renamed_CCLE_meth_means_mat)
filtSample_CCLE_meth_means_mat <- renamed_CCLE_meth_means_mat[rowSums(is.na(renamed_CCLE_meth_means_mat)) < CCLEsample_th,]
dim(filtSample_CCLE_meth_means_mat)  # 798  16,494

filtSample_CCLE_meth_means_mat <- filtSample_CCLE_meth_means_mat[,colSums(is.na(filtSample_CCLE_meth_means_mat)) < (0.1*nrow(filtSample_CCLE_meth_means_mat))]
dim(filtSample_CCLE_meth_means_mat)  # 798 15083


###### 4) Filtering by lineage ######
ccle_lineage_to_remove <- c('adrenal cortex', 'bile duct', 'bone', 'cervix',
                            'embryo', 'engineered', 'engineered blood', 'engineered bone',
                            'engineered breast engineered central nervous system', 'engineered kidney',
                            'engineered lung', 'engineered ovary', 'engineered prostate', 'eye',
                            'fibroblast', 'peripheral nervous system', 'soft tissue', 'unknown' )

ccle_lin_to_remove_names <- sample_info_CCLE %>% filter(lineage %in% ccle_lineage_to_remove) %>% pull(DepMap_ID)
length(ccle_lin_to_remove_names) #376

f_filtSample_CCLE_meth_means_mat <- filtSample_CCLE_meth_means_mat[!rownames(filtSample_CCLE_meth_means_mat) %in% ccle_lin_to_remove_names,]
dim(f_filtSample_CCLE_meth_means_mat) #731 15083

f_filtSample_CCLE_meth_means_mat <- as.matrix(f_filtSample_CCLE_meth_means_mat)
f_filtSample_CCLE_meth_means_mat[1:5, 1:5]

###### 5) Re-annotations #####
final_ccle_annot <- sample_info_CCLE %>%
  filter(DepMap_ID %in% rownames(f_filtSample_CCLE_meth_means_mat)) %>%
  mutate(lineage = ifelse(lineage == 'plasma cell', 'lymphocyte', lineage)
         )
unique(final_ccle_annot$lineage) %>% length
unique(final_ccle_annot$lineage) %>% sort
final_ccle_annot$lineage %>% sort %>% table
final_ccle_annot %>% pull(primary_or_metastasis) %>% sort %>% table

all.equal((final_tcga_annot$lineage %>% unique %>% sort),(final_ccle_annot$lineage %>% unique %>% sort))

###### 6) NA imputation #####

all(rownames(f_filtSample_CCLE_meth_means_mat) %in% final_ccle_annot$DepMap_ID)
pre_CCLE_meth_impute <- imputeNA_by_lineage(mat = f_filtSample_CCLE_meth_means_mat,
                                            annot = final_ccle_annot, sample_column = 'DepMap_ID')
dim(pre_CCLE_meth_impute) #  15083   731

nrow(pre_CCLE_meth_impute) == length(unique(rownames(pre_CCLE_meth_impute)))

#### Alignment ####

# ##### Intersect TCGA and CCLE at gene level #####
genes <- intersect(rownames(pre_TCGA_meth_impute), rownames(pre_CCLE_meth_impute))
length(genes) #12932

TCGA_meth_impute <- t(pre_TCGA_meth_impute[genes,])
CCLE_meth_impute <- t(pre_CCLE_meth_impute[genes,])

all(final_tcga_annot$sample.submitter_id %in% rownames(TCGA_meth_impute)) #TRUE
all(final_ccle_annot$DepMap_ID %in% rownames(CCLE_meth_impute)) #TRUE
all(rownames(TCGA_meth_impute) %in% final_tcga_annot$sample.submitter_id) #TRUE
all(rownames(CCLE_meth_impute) %in% final_ccle_annot$DepMap_ID) #TRUE
nrow(TCGA_meth_impute) == nrow(final_tcga_annot)
nrow(CCLE_meth_impute) == nrow(final_ccle_annot)

### Compute gene stats
common_genes_meth <- colnames(TCGA_meth_impute)

gene_stats_meth <- data.frame(
  Tumor_SD = apply(TCGA_meth_impute, 2, sd, na.rm=T),
  CCLE_SD = apply(CCLE_meth_impute, 2, sd, na.rm=T),
  Tumor_mean = colMeans(TCGA_meth_impute, na.rm=T),
  CCLE_mean = colMeans(CCLE_meth_impute, na.rm=T),
  Gene = common_genes_meth,
  stringsAsFactors = F) %>%
  dplyr::mutate(max_SD = pmax(Tumor_SD, CCLE_SD, na.rm=T)) #add avg and max SD per gene
dim(gene_stats_meth) # 12932   6

# the sampleID column must have the same order of the corresponding matrix's rownames
TCGA_ann <- final_tcga_annot %>% mutate(type='tumor', sampleID = sample.submitter_id) %>% select(sampleID, lineage, type) %>% as.data.frame()
CCLE_ann <- final_ccle_annot %>% mutate(type='CL', sampleID = DepMap_ID) %>% select(sampleID, lineage, type) %>% as.data.frame()

TCGA_ann <- TCGA_ann[match(rownames(TCGA_meth_impute), TCGA_ann$sampleID), ]
CCLE_ann <- CCLE_ann[match(rownames(CCLE_meth_impute), CCLE_ann$sampleID), ]

all.equal(rownames(TCGA_meth_impute), TCGA_ann$sampleID)
all.equal(rownames(CCLE_meth_impute), CCLE_ann$sampleID)

##### Create Seurat object: center the data + PCA + UMAP #####
TCGA_obj_meth <- create_Seurat_object(TCGA_meth_impute, TCGA_ann, type='tumor')
CCLE_obj_meth <- create_Seurat_object(CCLE_meth_impute, CCLE_ann, type='CL')

##### Find clusters with Louvain on PCA #####
TCGA_obj_meth <- cluster_data(TCGA_obj_meth) #Number of communities: 54
CCLE_obj_meth <- cluster_data(CCLE_obj_meth) #Number of communities: 25

##### DEA: Differential gene expresison analysis between clusters #####
tumor_DE_genes_meth <- find_differentially_expressed_genes(TCGA_obj_meth)
CL_DE_genes_meth <- find_differentially_expressed_genes(CCLE_obj_meth)

DE_genes_meth <- full_join(tumor_DE_genes_meth, CL_DE_genes_meth, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
  mutate(
    tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
    CL_rank = dplyr::dense_rank(-gene_stat_CL),
    best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
  dplyr::left_join(gene_stats_meth, by = 'Gene')

# take genes that are ranked in the top 1000 from either dataset, used for finding mutual nearest neighbors
DE_gene_set_meth <- DE_genes_meth %>%
  dplyr::filter(best_rank < MultiCelligner_parameters$top_DE_genes_per) %>%
  .[['Gene']]

length(unique(DE_gene_set_meth)) #1854

######################

comb_pre_mat <- rbind(TCGA_meth_impute, CCLE_meth_impute)
dim(comb_pre_mat)
comb_pre_mat <- t(comb_pre_mat)
dim(comb_pre_mat)

qnrm_comb <- preprocessCore::normalize.quantiles(comb_pre_mat)
colnames(qnrm_comb) <- colnames(comb_pre_mat)
rownames(qnrm_comb) <- rownames(comb_pre_mat)

CCLE_nrm <- qnrm_comb[,grepl(pattern = "ACH-00", x = colnames(qnrm_comb))]
TCGA_nrm <- qnrm_comb[,grepl(pattern = "TCGA", x = colnames(qnrm_comb))]

dim(CCLE_nrm)
dim(TCGA_nrm)

CCLE_nrm <- t(CCLE_nrm)
TCGA_nrm <- t(TCGA_nrm)

rownames(TCGA_nrm) <- str_remove(rownames(TCGA_nrm),"[A-Z]$")

TCGA_nrm <- TCGA_nrm[rownames(TCGA_nrm) %in% ann_multiomics_v9$sampleID,]
CCLE_nrm <- CCLE_nrm[rownames(CCLE_nrm) %in% ann_multiomics_v9$sampleID,]

### No cPCs correction before MNN
### MNN parameters associated with the max prop_agree in dist_all file

grid_all_max_y <- grid_all[which.max(grid_all$prop_agree_weigh_dist.y),]
grid_all_max_x <- grid_all[which.max(grid_all$prop_agree_weigh_dist.x),]

mnn_res <- run_MNN(CCLE_cor = CCLE_nrm, TCGA_cor = TCGA_nrm,
                   k1 = MultiCelligner_parameters$meth_mnn_k_tumor, k2 = MultiCelligner_parameters$meth_mnn_k_CL, ndist = MultiCelligner_parameters$ndist,
                   subset_genes = DE_gene_set_meth)

combined_mat_meth <- rbind(mnn_res$corrected,CCLE_nrm) 

saveRDS(combined_mat_meth, "DATA/Methylation/combined_mat_meth.rds")


