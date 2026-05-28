
#### Data availability ####
# Large input matrices are not distributed with the repository.
# Please download the required files and place them in the corresponding folders:
#
# DATA/
# ├── Annotation/
# │   └── ann_multiomics_v9.rds
# ├── Expression/
# │   └── combined_mat.rds
# ├── Methylation/
# │   └── combined_mat_meth.rds
# ├── Mutational_process/
# │   └── combined_mat_mut.rds
# └── Annotation/
#     └── ann_multiomics_v9.rds



library(tidyverse)
library(MOFA2)



###### Load annotation ######

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



###### Load combined matrices ######

expression_combined_file <- file.path(
  "DATA", "Expression", "combined_mat.rds"
)

if (!file.exists(expression_combined_file)) {
  stop(
    "Missing file: combined_mat.rds\n",
    "Please place the file in DATA/Expression/"
  )
}

combined_mat <- readRDS(expression_combined_file)



methylation_combined_file <- file.path(
  "DATA", "Methylation", "combined_mat_meth.rds"
)

if (!file.exists(methylation_combined_file)) {
  stop(
    "Missing file: combined_mat_meth.rds\n",
    "Please place the file in DATA/Methylation/"
  )
}

combined_mat_meth <- readRDS(methylation_combined_file)



mutation_combined_file <- file.path(
  "DATA", "Mutational_process", "combined_mat_mut.rds"
)

if (!file.exists(mutation_combined_file)) {
  stop(
    "Missing file: combined_mat_mut.rds\n",
    "Please place the file in DATA/Mutational_process/"
  )
}

combined_mat_mut <- readRDS(mutation_combined_file)

################# MOFA DATA INTEGRATION ################# 

#### Intersect samples among combined matrices 

omics_sample <- intersect(rownames(combined_mat), rownames(combined_mat_meth))
omics_sample_1 <- intersect(omics_sample, rownames(combined_mat_mut))

mofa_mat_exp <- combined_mat[which(rownames(combined_mat) %in% omics_sample_1),]
mofa_mat_meth <- combined_mat_meth[rownames(combined_mat_meth) %in% omics_sample_1,]
mofa_mat_mut <- combined_mat_mut[rownames(combined_mat_mut) %in% omics_sample_1,]

#### Subset samples from the original combined matrices; to use only for multiomics combination with Mutational signatures layer

mofa_mat_exp <- mofa_mat_exp[omics_sample_1,]
mofa_mat_meth <- mofa_mat_meth[omics_sample_1,]
mofa_mat_mut <- mofa_mat_mut[omics_sample_1,]

#### Subset samples from the original combined matrices; to use only for multiomics combination of DNA methylation and mRNA expression layers

mofa_mat_meth_exp <- combined_mat_meth[rownames(combined_mat_meth) %in% omics_sample,]
mofa_mat_exp_meth <- combined_mat[rownames(combined_mat) %in% omics_sample,]

mofa_mat_meth_exp <- mofa_mat_meth_exp[omics_sample,]
mofa_mat_exp_meth <- mofa_mat_exp_meth[omics_sample,]

#### Create different lists, one for each layer combination #### 

MOFA_exp_meth <- list(exp = t(mofa_mat_exp_meth),
                      meth = t(mofa_mat_meth_exp))

MOFA_all <- list(exp = t(mofa_mat_exp),
                 meth = t(mofa_mat_meth),
                 mut = t(mofa_mat_mut))

MOFA_exp_mut <- list(exp = t(mofa_mat_exp),
                     mut = t(mofa_mat_mut))

MOFA_meth_mut <- list(meth = t(mofa_mat_meth),
                      mut = t(mofa_mat_mut))

#### Create the final MOFA_list which contain lists of each layers combination #### 

MOFA_list <- list(
  MOFA_exp_mut = MOFA_exp_mut,
  MOFA_meth_mut = MOFA_meth_mut,
  MOFA_all = MOFA_all,
  MOFA_exp_meth = MOFA_exp_meth
)

################# Run MOFA: for each iteration of the for, MOFA will be executed on a single list of layers combination ################# 

MOFA_results <- list() 

for (j in names(MOFA_list)) {
  
  MOFA_obj <- MOFA_list[[j]]
  
  MOFAobject <- create_mofa(MOFA_obj)
  
  data_opts <- get_default_data_options(MOFAobject)
  
  model_opts <- get_default_model_options(MOFAobject)
  
  train_opts <- get_default_training_options(MOFAobject)
  
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  setwd("/DATA/SCRATCH/scala/celligner/1_multiCellignerFinal/DATA/multiomics/MOFA")
  outfile = file.path(getwd(),"MOFA_multiomics_model.hdf5")
  
  #### in parallel
  Sys.setenv(OMP_NUM_THREADS=20)
  Sys.setenv(RMKL_NUM_THREADS=20)
  
  MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
  
  model <- load_model('MOFA_multiomics_model.hdf5')
  
  #### subset with used samples the annotation file 
  if(j == MOFA_exp_meth) {
    
    mofa_ann <- ann_multiomics_v9[ann_multiomics_v9$sampleID %in% omics_sample, ]
    
  } else {
    
    mofa_ann <- ann_multiomics_v9[ann_multiomics_v9$sampleID %in% omics_sample_1, ]
  
  }
  
  sample_metadata <- mofa_ann %>% dplyr::select(sample=sampleID, dplyr::everything())
  
  #### check
  all(samples_names(MOFAobject)[[1]] == sample_metadata$sample)
  
  samples_metadata(model) <- sample_metadata
  
  model <- run_tsne(model)
  model <- run_umap(model)
  
  MOFA <- model@expectations$Z$group1 
  
  mofa_tsne <- model@dim_red$TSNE[,-1] %>% t()
  mofa_umap <- model@dim_red$UMAP[,-1] %>% t()
  
  MOFA_results[[j]] <- list(
    model = model,
    tsne = mofa_tsne,
    umap = mofa_umap,
    MOFA_mat = MOFA
  )
  
}

