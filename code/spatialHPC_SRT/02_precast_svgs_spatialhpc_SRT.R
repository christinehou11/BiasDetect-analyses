suppressPackageStartupMessages({
  library(Seurat)
  library(PRECAST)
  library(dplyr)
  library(here)
})
set.seed(123)

#build spe
counts = rbind(
    read.csv(here("processed-data","spatialHPC_SRT",
                "spe-hpc_sub4_svgs-only_counts-1.csv"), 
    row.names=1),
    read.csv(here("processed-data","spatialHPC_SRT",
                "spe-hpc_sub4_svgs-only_counts-2.csv"), 
    row.names=1))
colnames(counts) = gsub("\\.","-",colnames(counts))
cdata = read.csv(
    here("processed-data","spatialHPC_SRT",
        "spe-hpc_sub4_svgs-only_colData.csv"), 
    row.names=1)
rdata = read.csv(
    here("processed-data","spatialHPC_SRT",
        "spe-hpc_sub4_svgs-only_rowData.csv"), 
    row.names=1)

spe = SpatialExperiment(
    assay = list("counts"=counts), 
    colData = cdata, rowData = rdata,
    spatialCoordsNames = c("array_row", "array_col"))
spe
# class: SpatialExperiment 
# dim: 2082 18945 
# metadata(0):
#   assays(1): counts
# rownames(2082): ENSG00000002330 ENSG00000002586 ... ENSG00000286214
# ENSG00000286961
# rowData names(6): source type ... gene_name gene_type
# colnames(18945): AAACAACGAATAGTTC-1_V11L05-333_B1
# AAACAAGTATCTCCCA-1_V11L05-333_B1 ... TTGTTTGTATTACACG-1_V11L05-336_A1
# TTGTTTGTGTAAATTC-1_V11L05-336_A1
# colData names(10): sample_id slide ... expr_chrM expr_chrM_ratio
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : array_row array_col
# imgData names(0):

#reformat spe to seurat list
l2 = unique(spe$sample_id)
names(l2) = l2
l2 = lapply(l2, function(x) spe[,colData(spe)$sample_id==x])

srt.sets = lapply(l2, function(x) {
    colnames(counts(x)) <- rownames(colData(x))
    colData(x)$col <- spatialCoords(x)[,"array_col"]
    colData(x)$row <- spatialCoords(x)[,"array_row"]
    count <- counts(x)
    a1 <- CreateAssayObject(count, assay = "RNA", 
                            min.features = 0, min.cells = 0)
    CreateSeuratObject(a1, meta.data = as.data.frame(colData(x)))
})
srt.sets
# $`Br3942_V11L05-333_B1`
# An object of class Seurat 
# 2082 features across 4985 samples within 1 assay 
# Active assay: RNA (2082 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`Br3942_V11L05-333_D1`
# An object of class Seurat 
# 2082 features across 4938 samples within 1 assay 
# Active assay: RNA (2082 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`Br8325_V11L05-335_D1`
# An object of class Seurat 
# 2082 features across 4483 samples within 1 assay 
# Active assay: RNA (2082 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`Br8667_V11L05-336_A1`
# An object of class Seurat 
# 2082 features across 4539 samples within 1 assay 
# Active assay: RNA (2082 features, 0 variable features)
# 2 layers present: counts, data

#run precast
preobj <- CreatePRECASTObject(seuList = srt.sets,
                            customGenelist=rownames(spe),
                            premin.spots=0, premin.features=0, 
                            postmin.spots=0, postmin.features=0)
PRECASTObj <- AddAdjList(preobj, platform = "Visium")
PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, 
                            verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj <- PRECAST(PRECASTObj, K=7)

#consolidate/ reformat results
PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
# An object of class PRECASTObj 
# with 4 datasets and  18945 spots in total, with spots for each dataset:  4985 4938 4483 4539 
# 2081 common variable genes selected

seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
seuInt
# An object of class Seurat 
# 2081 features across 18945 samples within 1 assay 
# Active assay: PRE_CAST (2081 features, 0 variable features)
# 2 layers present: counts, data
# 2 dimensional reductions calculated: PRECAST, position

# -----------
# save object
# -----------
save(spe,
    file = here("processed-data","spatialHPC_SRT","spe.rda"))

write.csv(seuInt@meta.data, 
        here("processed-data","spatialHPC_SRT",
        "seuInt-hpc_k-7_svgs_metadata.csv"), 
        row.names=T)

# -----------
# session information
# -----------
sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS 15.0.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
#   [1] PRECAST_1.6.5               gtools_3.9.5               
# [3] Seurat_5.1.0                SeuratObject_5.0.2         
# [5] sp_2.1-4                    tibble_3.2.1               
# [7] tidyr_1.3.1                 dplyr_1.1.4                
# [9] SpatialExperiment_1.15.1    SingleCellExperiment_1.27.2
# [11] SummarizedExperiment_1.35.4 Biobase_2.65.1             
# [13] GenomicRanges_1.57.2        GenomeInfoDb_1.41.2        
# [15] IRanges_2.39.2              S4Vectors_0.43.2           
# [17] BiocGenerics_0.51.3         MatrixGenerics_1.17.0      
# [19] matrixStats_1.4.1           here_1.0.1                 
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.4.1           later_1.3.2            
# [4] bitops_1.0-9            polyclip_1.10-7         fastDummies_1.7.4      
# [7] lifecycle_1.0.4         rstatix_0.7.2           rprojroot_2.0.4        
# [10] globals_0.16.3          lattice_0.22-6          MASS_7.3-61            
# [13] backports_1.5.0         magrittr_2.0.3          plotly_4.10.4          
# [16] rmarkdown_2.28          yaml_2.3.10             httpuv_1.6.15          
# [19] sctransform_0.4.1       spam_2.11-0             spatstat.sparse_3.1-0  
# [22] reticulate_1.39.0       cowplot_1.1.3           pbapply_1.7-2          
# [25] RColorBrewer_1.1-3      abind_1.4-8             zlibbioc_1.51.1        
# [28] Rtsne_0.17              purrr_1.0.2             RCurl_1.98-1.16        
# [31] GenomeInfoDbData_1.2.13 ggrepel_0.9.6           irlba_2.3.5.1          
# [34] listenv_0.9.1           spatstat.utils_3.1-0    goftest_1.2-3          
# [37] RSpectra_0.16-2         spatstat.random_3.3-2   fitdistrplus_1.2-1     
# [40] parallelly_1.38.0       leiden_0.4.3.1          codetools_0.2-20       
# [43] DelayedArray_0.31.14    scuttle_1.15.4          tidyselect_1.2.1       
# [46] UCSC.utils_1.1.0        farver_2.1.2            viridis_0.6.5          
# [49] ScaledMatrix_1.13.0     spatstat.explore_3.3-2  jsonlite_1.8.9         
# [52] BiocNeighbors_1.99.2    Formula_1.2-5           progressr_0.14.0       
# [55] ggridges_0.5.6          survival_3.7-0          scater_1.33.4          
# [58] tools_4.4.1             ica_1.0-3               Rcpp_1.0.13            
# [61] glue_1.8.0              gridExtra_2.3           SparseArray_1.5.44     
# [64] xfun_0.48               ggthemes_5.1.0          withr_3.0.1            
# [67] fastmap_1.2.0           fansi_1.0.6             digest_0.6.37          
# [70] rsvd_1.0.5              R6_2.5.1                mime_0.12              
# [73] colorspace_2.1-1        scattermore_1.2         tensor_1.5             
# [76] spatstat.data_3.1-2     utf8_1.2.4              generics_0.1.3         
# [79] data.table_1.16.2       httr_1.4.7              htmlwidgets_1.6.4      
# [82] S4Arrays_1.5.11         uwot_0.2.2              pkgconfig_2.0.3        
# [85] gtable_0.3.5            lmtest_0.9-40           XVector_0.45.0         
# [88] htmltools_0.5.8.1       carData_3.0-5           dotCall64_1.2          
# [91] scales_1.3.0            png_0.1-8               spatstat.univar_3.0-1  
# [94] knitr_1.48              rstudioapi_0.16.0       reshape2_1.4.4         
# [97] rjson_0.2.23            nlme_3.1-166            zoo_1.8-12             
# [100] stringr_1.5.1           KernSmooth_2.23-24      vipor_0.4.7            
# [103] miniUI_0.1.1.1          GiRaF_1.0.1             pillar_1.9.0           
# [106] grid_4.4.1              vctrs_0.6.5             RANN_2.6.2             
# [109] ggpubr_0.6.0            promises_1.3.0          car_3.1-3              
# [112] BiocSingular_1.21.4     DR.SC_3.4               beachmat_2.21.6        
# [115] xtable_1.8-4            cluster_2.1.6           beeswarm_0.4.0         
# [118] evaluate_1.0.1          magick_2.8.5            cli_3.6.3              
# [121] compiler_4.4.1          rlang_1.1.4             crayon_1.5.3           
# [124] ggsignif_0.6.4          future.apply_1.11.2     mclust_6.1.1           
# [127] plyr_1.8.9              ggbeeswarm_0.7.2        stringi_1.8.4          
# [130] viridisLite_0.4.2       deldir_2.0-4            BiocParallel_1.39.0    
# [133] munsell_0.5.1           lazyeval_0.2.2          spatstat.geom_3.3-3    
# [136] CompQuadForm_1.4.3      Matrix_1.7-0            RcppHNSW_0.6.0         
# [139] patchwork_1.3.0         future_1.34.0           ggplot2_3.5.1          
# [142] shiny_1.9.1             ROCR_1.0-11             broom_1.0.7            
# [145] igraph_2.0.3  
