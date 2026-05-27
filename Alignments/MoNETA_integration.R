
library(tidyverse)
library(MoNETA)

source(utils)
source(global_parameters)

combined_mat <- readRDS("/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/Expression/combined_mat.rds")
combined_mat_meth <- readRDS("/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/Methylation/combined_mat_meth.rds")
combined_mat_mut <- readRDS("/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/Mutational_process/combined_mat_mut.rds")

################# MoNETA DATA INTEGRATION ################# 

##### generate a weighted adjacency list for each omics matrices combination ##### 

#### DNA methylation - mRNA expression - Mutational signatures 
MoNETA_exp_meth_mut <- list(exp_data = k_star_net(t(combined_mat),
                                                  sparsity = .7,
                                                  distFun = "Euclidean",
                                                  cores = 60,
                                                  knn = MultiCelligner_parameters$MoNETA_knn,
                                                  MAX_ASSOC = MultiCelligner_parameters$MoNETA_MAX_ASSOC),
                            meth_data = k_star_net(t(combined_mat_meth),
                                                   sparsity = .7,
                                                   distFun = "Euclidean",
                                                   cores = 60,
                                                   knn = MultiCelligner_parameters$MoNETA_knn,
                                                   MAX_ASSOC = MultiCelligner_parameters$MoNETA_MAX_ASSOC),
                            mut_data = k_star_net(t(combined_mat_mut),
                                                  sparsity = .7,
                                                  distFun = "Euclidean",
                                                  cores = 60,
                                                  knn = MultiCelligner_parameters$MoNETA_knn,
                                                  MAX_ASSOC = MultiCelligner_parameters$MoNETA_MAX_ASSOC))

#### DNA methylation - mRNA expression 
MoNETA_exp_meth <- list(exp_data = k_star_net(t(combined_mat),
                                              sparsity = .7,
                                              distFun = "Euclidean",
                                              cores = 60,
                                              knn = MultiCelligner_parameters$MoNETA_knn,
                                              MAX_ASSOC = MultiCelligner_parameters$MoNETA_MAX_ASSOC),
                        meth_data = k_star_net(t(combined_mat_meth),
                                               sparsity = .7,
                                               distFun = "Euclidean",
                                               cores = 60,
                                               knn = MultiCelligner_parameters$MoNETA_knn,
                                               MAX_ASSOC = MultiCelligner_parameters$MoNETA_MAX_ASSOC))

#### mRNA expression - Mutational signatures 
MoNETA_exp_mut <- list(exp_data = k_star_net(t(combined_mat),
                                             sparsity = .7,
                                             distFun = "Euclidean",
                                             cores = 60,
                                             knn = MultiCelligner_parameters$MoNETA_knn,
                                             MAX_ASSOC = MultiCelligner_parameters$MoNETA_MAX_ASSOC),
                       mut_data = k_star_net(t(combined_mat_mut),
                                             sparsity = .7,
                                             distFun = "Euclidean",
                                             cores = 60,
                                             knn = MultiCelligner_parameters$MoNETA_knn,
                                             MAX_ASSOC = MultiCelligner_parameters$MoNETA_MAX_ASSOC))

#### DNA methylation - Mutational signatures 
MoNETA_meth_mut <- list(meth_data = k_star_net(t(combined_mat_meth),
                                               sparsity = .7,
                                               distFun = "Euclidean",
                                               cores = 60,
                                               knn = MultiCelligner_parameters$MoNETA_knn,
                                               MAX_ASSOC = MultiCelligner_parameters$MoNETA_MAX_ASSOC),
                        mut_data = k_star_net(t(combined_mat_mut),
                                              sparsity = .7,
                                              distFun = "Euclidean",
                                              cores = 60,
                                              knn = MultiCelligner_parameters$MoNETA_knn,
                                              MAX_ASSOC = MultiCelligner_parameters$MoNETA_MAX_ASSOC))


##### Get the multiplex network by combining the list of single-omics networks ##### 

#### DNA methylation - mRNA expression - Mutational signatures 
multiplex_exp_meth_mut <- create_multiplex(MoNETA_exp_meth_mut)

#### DNA methylation - mRNA expression 
multiplex_exp_meth <- create_multiplex(MoNETA_exp_meth)

#### mRNA expression - Mutational signatures 
multiplex_exp_mut <- create_multiplex(MoNETA_exp_mut)

#### DNA methylation - Mutational signatures 
multiplex_meth_mut <- create_multiplex(MoNETA_meth_mut)

##### Get the probabilities of transitioning between different omics layers for each multiplex ##### 

#### DNA methylation - mRNA expression 
layer_transition_1  <-  create_layer_transition_matrix(MoNETA_exp_meth)

#### mRNA expression - Mutational signatures 
layer_transition_2  <-  create_layer_transition_matrix(MoNETA_exp_mut)

#### DNA methylation - Mutational signatures 
layer_transition_3  <-  create_layer_transition_matrix(MoNETA_meth_mut)

#### DNA methylation - mRNA expression - Mutational signatures 
layer_transition_4 <- create_layer_transition_matrix(MoNETA_exp_meth_mut)

##### Compute the probability distribution of reaching other nodes in each multiplex by a random walk process ##### 

#### DNA methylation - mRNA expression - Mutational signatures 
RWR_mat_exp_meth_mut <- gen_sim_mat_M(network = multiplex_exp_meth_mut,
                                      tau = tau_multiomics_3, restart = 0.1,
                                      layer_transition = layer_transition_4,
                                      jump_neighborhood = F,
                                      weighted_multiplex = F,
                                      cores = 80)

#### DNA methylation - mRNA expression 
RWR_mat_exp_meth <- gen_sim_mat_M(network = multiplex_exp_meth,
                                  tau = tau_multiomics_1, restart = 0.1,
                                  layer_transition = layer_transition_1,
                                  jump_neighborhood = F,
                                  weighted_multiplex = F,
                                  cores = 80)

#### mRNA expression - Mutational signatures 
RWR_mat_exp_mut <- gen_sim_mat_M(network = multiplex_exp_mut,
                                 tau = tau_multiomics_2, restart = 0.1,
                                 layer_transition = layer_transition_2,
                                 jump_neighborhood = F,
                                 weighted_multiplex = F,
                                 cores = 80)

#### DNA methylation - Mutational signatures 
RWR_mat_meth_mut <- gen_sim_mat_M(network = multiplex_meth_mut,
                                  tau = tau_multiomics_2, restart = 0.1,
                                  layer_transition = layer_transition_3,
                                  jump_neighborhood = F,
                                  weighted_multiplex = F,
                                  cores = 80)

saveRDS(RWR_mat_exp_meth ,'/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/RWR_mat_exp_meth.rds')
saveRDS(RWR_mat_exp_mut ,'/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/RWR_mat_exp_mut.rds')
saveRDS(RWR_mat_meth_mut ,'/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/RWR_mat_meth_mut.rds')
saveRDS(RWR_mat_exp_meth_mut ,'/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/RWR_mat_exp_meth_mut.rds')

################# Dimensionally reductions ################# 

#### Embedding #### 

emb_exp_meth_1 <- MoNETA::get_embedding(RWR_mat_exp_meth, embedding_size = 70, cores = 80)
emb_exp_mut_1 <- MoNETA::get_embedding(RWR_mat_exp_mut, embedding_size = 70, cores = 80)
emb_meth_mut_1 <- MoNETA::get_embedding(RWR_mat_meth_mut, embedding_size = 70, cores = 80)
emb_exp_meth_mut_1 <- MoNETA::get_embedding(RWR_mat_exp_meth_mut, embedding_size = 70, cores = 80)

saveRDS(emb_exp_meth_1,'/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/emb_exp_meth_1.rds')
saveRDS(emb_exp_mut_1,'/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/emb_exp_mut_1.rds')
saveRDS(emb_meth_mut_1,'/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/emb_meth_mut_1.rds')
saveRDS(emb_exp_meth_mut_1,'/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/emb_exp_meth_mut_1.rds')

#### UMAP #### 

umap_exp_meth_1 <- get_parallel_umap_embedding(emb_exp_meth_1,
                                               embedding_size = 2,
                                               n_threads = 80,
                                               metric = 'euclidean')

umap_exp_mut_1 <- get_parallel_umap_embedding(emb_exp_mut_1,
                                              embedding_size = 2,
                                              n_threads = 80,
                                              metric = 'euclidean')

umap_meth_mut_1 <- get_parallel_umap_embedding(emb_meth_mut_1,
                                               embedding_size = 2,
                                               n_threads = 80,
                                               metric = 'euclidean')

umap_exp_meth_mut_1 <- get_parallel_umap_embedding(emb_exp_meth_mut_1,
                                                   embedding_size = 2,
                                                   n_threads = 80,
                                                   metric = 'euclidean')

#### tSNE #### 

tsne_exp_meth <- get_tsne_embedding(RWR_mat_exp_meth, 
                                    embedding_size = 2, 
                                    perplexity = 70, 
                                    max_iter = 20000 , 
                                    num_threads = 80)

tsne_exp_meth_mut <- get_tsne_embedding(RWR_mat_exp_meth_mut, 
                                        embedding_size = 2, 
                                        perplexity = 70, 
                                        max_iter = 20000 , 
                                        num_threads = 80)

tsne_exp_mut <- get_tsne_embedding(RWR_mat_exp_mut, 
                                   embedding_size = 2, 
                                   perplexity = 70, 
                                   max_iter = 20000 , 
                                   num_threads = 80)

tsne_meth_mut <- get_tsne_embedding(RWR_mat_meth_mut, 
                                    embedding_size = 2, 
                                    perplexity = 70, 
                                    max_iter = 20000 , 
                                    num_threads = 80)

saveRDS(umap_exp_meth_1, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/umap_exp_meth.rds')
saveRDS(umap_exp_mut_1, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/umap_exp_mut.rds')
saveRDS(umap_meth_mut_1, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/umap_meth_mut.rds')
saveRDS(umap_exp_meth_mut_1, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/umap_exp_meth_mut.rds')

saveRDS(tsne_exp_meth, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/tsne_exp_meth.rds')
saveRDS(tsne_exp_meth_mut, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/tsne_exp_meth_mut.rds')
saveRDS(tsne_exp_mut, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/tsne_exp_mut.rds')
saveRDS(tsne_meth_mut, '/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MoNETA/tsne_meth_mut.rds')











