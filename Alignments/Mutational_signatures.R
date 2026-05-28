
#### Data availability ####
# Large input matrices are not distributed with the repository.
# Please download the required files and place them in the corresponding folders:
#
# DATA/
# ├── Mutational_process/
# │   ├── CCLE_mutations.csv
# │   └── mut_pancancer.rds
# └── Annotation/
#     └── ann_multiomics_v9.rds



library(tidyverse)
library(sigminer)
library(maftools)
library(foreach)
library(snow)
library(doParallel)

source(utils)
source(global_parameters)

###### Get combined_mat for Mutational Signatures data ###### 

###### Load Data ######

ccle_mut_file <- file.path(
  "DATA", "Mutational_process", "CCLE_mutations.csv"
)

if (!file.exists(ccle_mut_file)) {
  stop(
    "Missing file: CCLE_mutations.csv\n",
    "Please place the file in DATA/Mutational_process/"
  )
}

CCLE_mutations <- read.csv(ccle_mut_file)



mut_pancancer_file <- file.path(
  "DATA", "Mutational_process", "mut_pancancer.rds"
)

if (!file.exists(mut_pancancer_file)) {
  stop(
    "Missing file: mut_pancancer.rds\n",
    "Please place the file in DATA/Mutational_process/"
  )
}

mut_pancancer <- readRDS(mut_pancancer_file)



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

mut_ccle_maf_type_1 <- CCLE_mutations %>%  dplyr::select(DepMap_ID, Chromosome, Start_position, End_position,Reference_Allele,
                                                         Alternate_Allele,Hugo_Symbol, Variant_Classification,
                                                         Variant_Type, Variant_annotation)

colnames(mut_ccle_maf_type_1)[c(1,6,3,4)] <- c("Tumor_Sample_Barcode", "Tumor_Seq_Allele2","Start_Position","End_Position")

### FILTER MUTATIONS: keep only SNP for both CCLE and TCGA samples
mut_ccle_maf_type_2 <- mut_ccle_maf_type_1 %>% dplyr::filter(Variant_Type %in% "SNP")

mut_pancancer_2 <- mut_pancancer %>% dplyr::filter(effect %in% c("RNA","Missense_Mutation", "5'Flank", "3'Flank", "3'UTR",
                                                                 "Nonsense_Mutation", "Splice_Site", "Intron", "5'UTR", "Nonstop_Mutation", 
                                                                 "3'Flank", "Translation_Start_Site"))

### FILTER MUTATIONS: keep only SNP for both CCLE and TCGA samples
mut_pancancer_3 <- mut_pancancer_2 %>% mutate("Variant_Type" = rep(x = "SNP"))

halfmaf <- mut_pancancer_3 %>% dplyr::rename(Hugo_Symbol = gene, Chromosome = chr, 
                                             Start_position = start, End_position = end, Tumor_Seq_Allele2 = alt, Reference_Allele = reference,
                                             Variant_Classification = effect, Variant_Type = Variant_Type, Tumor_Sample_Barcode = sample)

mut_pancancer_maf <- maftools::read.maf(maf = halfmaf)

####################### Mutational Signatures REFERENCE FITTING ####################### 

###### Classifying SBS records into 96 components ###### 
mut_ccle_maf <- sigminer::read_maf(maf = mut_ccle_maf_type_1,
                                   verbose = TRUE)

pancancer_tally <- sigminer::sig_tally(mut_pancancer_maf,
                                       ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                                       use_syn = TRUE)

ccle_tally <- sigminer::sig_tally(mut_ccle_maf,
                                  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                                  use_syn = TRUE)

###### Compute exposure of pre-defined signatures for each TCGA-CCLE samples ###### 

mut_fit_sign_ccle <- get_fitting_signature(ccle_tally$nmf_matrix, cores = 30)
mut_fit_sign_pancancer <- get_fitting_signature(pancancer_tally$nmf_matrix, cores = 50)

mut_fit_sign_pancancer_v1 <- mut_fit_sign_pancancer[,colnames(mut_fit_sign_pancancer) %in% ann_multiomics_v9$sampleID]
mut_fit_sign_ccle_v1 <- mut_fit_sign_ccle[,colnames(mut_fit_sign_ccle) %in% ann_multiomics_v9$sampleID]

mut_fit_sign_pancancer_v1 <- t(mut_fit_sign_pancancer_v1)
mut_fit_sign_ccle_v1 <- t(mut_fit_sign_ccle_v1)

combined_mat_mut <- rbind(mut_fit_sign_ccle_v1, mut_fit_sign_pancancer_v1)

saveRDS(combined_mat_mut, "DATA/Mutational_process/combined_mat_mut.rds")









