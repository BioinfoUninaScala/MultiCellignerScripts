
library(tidyverse)
library(SNFtool)
library(uwot)
library(MoNETA)

combined_mat <- readRDS("/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/Expression/combined_mat.rds")
combined_mat_meth <- readRDS("/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/Methylation/combined_mat_meth.rds")
combined_mat_mut <- readRDS("/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/Mutational_process/combined_mat_mut.rds")

################# MOFA DATA INTEGRATION ################# 

#### Intersect samples among combined matrices 

sample_exp_meth <- dplyr::intersect(rownames(combined_mat), rownames(combined_mat_meth))
sample_exp_mut <- dplyr::intersect(rownames(combined_mat), rownames(combined_mat_mut))
sample_meth_mut <- dplyr::intersect(rownames(combined_mat_meth), rownames(combined_mat_mut))
sample_all <- dplyr::intersect(rownames(combined_mat_mut), sample_exp_meth)

#### Subset samples from the original combined matrices

mat_exp_meth <- combined_mat[rownames(combined_mat) %in% sample_exp_meth,] 
mat_exp_mut <- combined_mat[rownames(combined_mat) %in% sample_exp_mut,] 
mat_meth_exp <- combined_mat_meth[rownames(combined_mat_meth) %in% sample_exp_meth,] 
mat_meth_mut <- combined_mat_meth[rownames(combined_mat_meth) %in% sample_meth_mut,] 
mat_mut_exp <- combined_mat_mut[rownames(combined_mat_mut) %in% sample_exp_mut,] 
mat_mut_meth <- combined_mat_mut[rownames(combined_mat_mut) %in% sample_meth_mut,] 

mat_all_exp <-combined_mat[rownames(combined_mat) %in% sample_all,] 
mat_all_meth <- combined_mat_meth[rownames(combined_mat_meth) %in% sample_all,] 
mat_all_mut <- combined_mat_mut[rownames(combined_mat_mut) %in% sample_all,] 

#### Compute pairwise squared Euclidean distances between tumor/CCL sample.

dist_exp_meth <- SNFtool::dist2(mat_exp_meth, mat_exp_meth)
dist_exp_mut <- SNFtool::dist2(mat_exp_mut, mat_exp_mut)
dist_meth_exp <- SNFtool::dist2(mat_meth_exp, mat_meth_exp) 
dist_meth_mut <- SNFtool::dist2(mat_meth_mut, mat_meth_mut)
dist_mut_exp <- SNFtool::dist2(mat_mut_exp, mat_mut_exp)
dist_mut_meth <- SNFtool::dist2(mat_mut_meth, mat_mut_meth)

dist_all_exp <- SNFtool::dist2(mat_all_exp, mat_all_exp)
dist_all_meth <- SNFtool::dist2(mat_mut_meth, mat_mut_meth)
dist_all_mut <- SNFtool::dist2(mat_all_mut, mat_all_mut)

#### Compute the affinity matrices that represents the neighborhood graph of each tumor/CCL sample.

aff_exp_meth <- SNFtool::affinityMatrix(dist_exp_meth)
aff_exp_mut <- SNFtool::affinityMatrix(dist_exp_mut)
aff_meth_exp <- SNFtool::affinityMatrix(dist_meth_exp) 
aff_meth_mut <- SNFtool::affinityMatrix(dist_meth_mut)
aff_mut_exp <- SNFtool::affinityMatrix(dist_mut_exp)
aff_mut_meth <- SNFtool::affinityMatrix(dist_mut_meth)

aff_exp_meth <- aff_exp_meth[rownames(aff_meth_exp),]
aff_exp_meth <- aff_exp_meth[,colnames(aff_meth_exp)]

aff_exp_mut <- aff_exp_mut[rownames(aff_mut_exp),]
aff_exp_mut <- aff_exp_mut[,colnames(aff_mut_exp)]

aff_meth_mut <- aff_meth_mut[rownames(aff_mut_meth),]
aff_meth_mut <- aff_meth_mut[,colnames(aff_mut_meth)]

aff_all_exp <- SNFtool::affinityMatrix(dist_all_exp)
aff_all_meth <- SNFtool::affinityMatrix(dist_all_meth)
aff_all_mut <- SNFtool::affinityMatrix(dist_all_mut)

aff_all_meth <- aff_all_meth[rownames(aff_all_exp),]
aff_all_meth <- aff_all_meth[,colnames(aff_all_exp)]
aff_all_mut <- aff_all_mut[rownames(aff_all_exp),]
aff_all_mut <- aff_all_mut[,colnames(aff_all_exp)]

#### Fuse single network together to construct the results multiomics matrices 

SNF_exp_meth <-  SNFtool::SNF(list(aff_exp_meth,aff_meth_exp))
SNF_exp_mut <- SNFtool::SNF(list(aff_exp_mut,aff_mut_exp))
SNF_meth_mut <- SNFtool::SNF(list(aff_meth_mut,aff_mut_meth))
SNF_all <- SNFtool::SNF(list(aff_all_exp,aff_all_meth,aff_all_mut))

saveRDS(SNF_all, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/SNF_all.rds')
saveRDS(SNF_exp_meth, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/SNF_exp_meth.rds')
saveRDS(SNF_exp_mut, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/SNF_exp_mut.rds')
saveRDS(SNF_meth_mut, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/SNF_meth_mut.rds')

#### UMAP #### 

snf_umap_all <- uwot::umap(SNF_all, metric = "cosine")
snf_umap_exp_meth <- uwot::umap(SNF_exp_meth, metric = "cosine")
snf_umap_exp_mut <- uwot::umap(SNF_exp_mut, metric = "cosine")
snf_umap_meth_mut <- uwot::umap(SNF_meth_mut, metric = "cosine")

saveRDS(snf_umap_all, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/snf_umap_all.rds')
saveRDS(snf_umap_exp_meth, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/snf_umap_exp_meth.rds')
saveRDS(snf_umap_exp_mut, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/snf_umap_exp_mut.rds')
saveRDS(snf_umap_meth_mut, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/snf_umap_meth_mut.rds')

#### tSNE #### 

snf_tsne_exp_meth <- get_tsne_embedding(SNF_exp_meth, 
                                        embedding_size = 2, 
                                        perplexity = 70, 
                                        max_iter = 20000 , 
                                        num_threads = 80)

snf_tsne_all <- get_tsne_embedding(SNF_all, 
                                   embedding_size = 2, 
                                   perplexity = 70, 
                                   max_iter = 20000 , 
                                   num_threads = 80)

snf_tsne_exp_mut <- get_tsne_embedding(SNF_exp_mut, 
                                       embedding_size = 2, 
                                       perplexity = 70, 
                                       max_iter = 20000 , 
                                       num_threads = 80)

snf_tsne_meth_mut <- get_tsne_embedding(SNF_meth_mut, 
                                        embedding_size = 2, 
                                        perplexity = 70, 
                                        max_iter = 20000 , 
                                        num_threads = 80)

saveRDS(snf_tsne_exp_meth, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/snf_tsne_exp_meth.rds')
saveRDS(snf_tsne_all, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/snf_tsne_all.rds')
saveRDS(snf_tsne_exp_mut, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/snf_tsne_exp_mut.rds')
saveRDS(snf_tsne_meth_mut, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/SNF/snf_tsne_meth_mut.rds')

