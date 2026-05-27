
################# MultiCelligner global parameters ################# 

MultiCelligner_parameters <- list(
  # tcga_DNAmethylation_matrix = ,
  # ccle_DNAmethylation_matrix = ,
  # tcga_mRNAExpression_matrix = ,
  # ccle_mRNAExpression_matrix = ,
  # tcga_Mut_sign_matrix = ,
  # ccle_Mut_sign_matrix = ,
  meth_mnn_k_CL = 10, # number of nearest neighbors of tumors in the cell line DNA methylation data
  meth_mnn_k_tumor = 45, # number of nearest neighbors of cell lines in the tumor DNA methylation data
  exp_mnn_k_CL = 5, # number of nearest neighbors of tumors in the cell line mRNA expression data
  exp_mnn_k_tumor = 50, # number of nearest neighbors of cell lines in the tumor mRNA expression data
  mnn_ndist = 3, # ndist parameter used for MNN (threshold beyond which neighbours are to be ignored when computing correction vectors)
  top_DE_genes_per = 1000, # differentially expressed genes with a rank better than this is in the cell line or tumor data
  # are used to identify mutual nearest neighbors in the MNN alignment step
  remove_cPCA_dims = c(1,2,3,4), # which cPCA dimensions to regress out of the data
  mod_clust_res = 5, # resolution parameter used for clustering the data
  n_PC_dims = 70, # number of PCs to use for dimensionality reduction
  reduction.use = 'umap', # 2D projection used for plotting
  fast_cPCA = 10, # to run fast cPCA (approximate the cPCA eigenvectors instead of calculating all) set this to a value >= 4
  MoNETA_knn = 5, # maximum number of neighbors to be considered as candidate neighbors for each node in 5he MoNETA single omic similarity network construction 
  MoNETA_MAX_ASSOC = 5, # number of maximum incoming edges that a node can have in the MoNETA single omic similarity network 
  MoNETA_tau_meth_exp =  c(0.5,0.5), # tau values specify likelihood of swapping in that omic layers for derived the final patient multiomics similarity matrix 
  MoNETA_tau_y_mut = c(0.7,0.3), # tau values specify likelihood of swapping in that omic layers for derived the final patient multiomics similarity matrix; y indicate both meth and exp 
  MoNETA_tau_meth_exp_mut = c(0.4,0.4,0.2) # tau values specify likelihood of swapping in that omic layers for derived the final patient multiomics similarity matrix 

)
