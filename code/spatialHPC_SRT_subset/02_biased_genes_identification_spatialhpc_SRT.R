suppressPackageStartupMessages({
    library(BiasDetect)
    library(dplyr)
    library(here)
    library(SummarizedExperiment)
})

# ---------
# load data
# ---------
load(here("processed-data","spatialHPC_SRT","spe_sub4.rda"))

# ---------
# biased gene identification
# ---------
SVGs <- rowData(spe_sub4)$gene_id
batch_df <- BiasDetect::featureSelect(spe_sub4, 
                                    batch_effect = "sample_id", VGs = SVGs)
head(batch_df)
# gene gene_name dev_default rank_default dev_batch rank_batch
# 1 ENSG00000002330       BAD    15830.43         1310  15785.12       1289
# 2 ENSG00000002586      CD99    18837.21          865  18616.10        858
# 3 ENSG00000004059      ARF5    15033.32         1443  14985.57       1417
# 4 ENSG00000004660    CAMKK1    14586.54         1506  14479.27       1487
# 5 ENSG00000004779   NDUFAB1    19859.78          756  19469.72        764
# 6 ENSG00000005022   SLC25A5    16842.03         1150  16739.61       1131

# d_diff     nSD_dev r_diff   nSD_rank
# 1 0.002870288 -0.31052071    -21 -0.4914075
# 2 0.011877487 -0.13869673     -7 -0.1638025
# 3 0.003186151 -0.30449522    -26 -0.6084093
# 4 0.007408683 -0.22394495    -19 -0.4446068
# 5 0.020034113  0.01690146      8  0.1872029
# 6 0.006118529 -0.24855631    -19 -0.4446068

# list bias as data frame
bias <- BiasDetect::biasDetect(batch_df = batch_df, nSD_dev = 10, nSD_rank = 5)

# display bias in plots
bias_plots <- BiasDetect::biasDetect(batch_df = batch_df, 
                        nSD_dev = 10, nSD_rank = 5, visual = TRUE)
# deviance
png(here("plots", "spatialHPC_SRT","biased_iden_dev.png"), 
      width=5, height=5, units="in", res=300)
p1_dev <- bias_plots$deviance
p1_dev
dev.off()

# rank
png(here("plots", "spatialHPC_SRT","biased_iden_rank.png"), 
    width=5, height=5, units="in", res=300)
p2_rank <- bias_plots$rank
p2_rank
dev.off()

# -----------
# save object
# -----------
write.csv(batch_df, 
    here("processed-data","spatialHPC_SRT",
        "hpc_bindev_default-sample_svgs-only.csv"), 
    row.names=FALSE)

save(bias,
    file = here("processed-data","spatialHPC_SRT","hpc_biased_features.rda"))

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
#   [1] BiasDetect_0.99.0           PRECAST_1.6.5              
# [3] gtools_3.9.5                Seurat_5.1.0               
# [5] SeuratObject_5.0.2          sp_2.1-4                   
# [7] tibble_3.2.1                tidyr_1.3.1                
# [9] dplyr_1.1.4                 SpatialExperiment_1.15.1   
# [11] SingleCellExperiment_1.27.2 SummarizedExperiment_1.35.4
# [13] Biobase_2.65.1              GenomicRanges_1.57.2       
# [15] GenomeInfoDb_1.41.2         IRanges_2.39.2             
# [17] S4Vectors_0.43.2            BiocGenerics_0.51.3        
# [19] MatrixGenerics_1.17.0       matrixStats_1.4.1          
# [21] here_1.0.1                 
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.4.1           later_1.3.2            
# [4] bitops_1.0-9            polyclip_1.10-7         fastDummies_1.7.4      
# [7] lifecycle_1.0.4         rstatix_0.7.2           rprojroot_2.0.4        
# [10] processx_3.8.4          globals_0.16.3          lattice_0.22-6         
# [13] MASS_7.3-61             backports_1.5.0         magrittr_2.0.3         
# [16] plotly_4.10.4           rmarkdown_2.28          remotes_2.5.0          
# [19] yaml_2.3.10             httpuv_1.6.15           sctransform_0.4.1      
# [22] spam_2.11-0             pkgbuild_1.4.4          spatstat.sparse_3.1-0  
# [25] reticulate_1.39.0       cowplot_1.1.3           pbapply_1.7-2          
# [28] RColorBrewer_1.1-3      abind_1.4-8             zlibbioc_1.51.1        
# [31] Rtsne_0.17              purrr_1.0.2             RCurl_1.98-1.16        
# [34] GenomeInfoDbData_1.2.13 ggrepel_0.9.6           scry_1.17.0            
# [37] irlba_2.3.5.1           listenv_0.9.1           spatstat.utils_3.1-0   
# [40] goftest_1.2-3           RSpectra_0.16-2         spatstat.random_3.3-2  
# [43] fitdistrplus_1.2-1      parallelly_1.38.0       leiden_0.4.3.1         
# [46] codetools_0.2-20        DelayedArray_0.31.14    scuttle_1.15.4         
# [49] tidyselect_1.2.1        UCSC.utils_1.1.0        farver_2.1.2           
# [52] viridis_0.6.5           ScaledMatrix_1.13.0     spatstat.explore_3.3-2 
# [55] jsonlite_1.8.9          BiocNeighbors_1.99.2    Formula_1.2-5          
# [58] progressr_0.14.0        ggridges_0.5.6          survival_3.7-0         
# [61] scater_1.33.4           tools_4.4.1             ica_1.0-3              
# [64] Rcpp_1.0.13             glue_1.8.0              gridExtra_2.3          
# [67] SparseArray_1.5.44      xfun_0.48               ggthemes_5.1.0         
# [70] withr_3.0.1             BiocManager_1.30.25     fastmap_1.2.0          
# [73] fansi_1.0.6             callr_3.7.6             digest_0.6.37          
# [76] rsvd_1.0.5              R6_2.5.1                mime_0.12              
# [79] colorspace_2.1-1        scattermore_1.2         tensor_1.5             
# [82] spatstat.data_3.1-2     utf8_1.2.4              generics_0.1.3         
# [85] data.table_1.16.2       httr_1.4.7              htmlwidgets_1.6.4      
# [88] S4Arrays_1.5.11         uwot_0.2.2              pkgconfig_2.0.3        
# [91] gtable_0.3.5            lmtest_0.9-40           XVector_0.45.0         
# [94] htmltools_0.5.8.1       carData_3.0-5           dotCall64_1.2          
# [97] scales_1.3.0            png_0.1-8               spatstat.univar_3.0-1  
# [100] knitr_1.48              rstudioapi_0.16.0       reshape2_1.4.4         
# [103] rjson_0.2.23            curl_5.2.3              nlme_3.1-166           
# [106] zoo_1.8-12              stringr_1.5.1           KernSmooth_2.23-24     
# [109] vipor_0.4.7             miniUI_0.1.1.1          GiRaF_1.0.1            
# [112] desc_1.4.3              pillar_1.9.0            grid_4.4.1             
# [115] vctrs_0.6.5             RANN_2.6.2              ggpubr_0.6.0           
# [118] promises_1.3.0          car_3.1-3               BiocSingular_1.21.4    
# [121] DR.SC_3.4               beachmat_2.21.6         xtable_1.8-4           
# [124] cluster_2.1.6           beeswarm_0.4.0          evaluate_1.0.1         
# [127] magick_2.8.5            cli_3.6.3               compiler_4.4.1         
# [130] rlang_1.1.4             crayon_1.5.3            ggsignif_0.6.4         
# [133] future.apply_1.11.2     labeling_0.4.3          mclust_6.1.1           
# [136] ps_1.8.0                plyr_1.8.9              ggbeeswarm_0.7.2       
# [139] stringi_1.8.4           viridisLite_0.4.2       deldir_2.0-4           
# [142] BiocParallel_1.39.0     munsell_0.5.1           lazyeval_0.2.2         
# [145] spatstat.geom_3.3-3     CompQuadForm_1.4.3      Matrix_1.7-0           
# [148] RcppHNSW_0.6.0          patchwork_1.3.0         future_1.34.0          
# [151] ggplot2_3.5.1           shiny_1.9.1             ROCR_1.0-11            
# [154] broom_1.0.7             igraph_2.0.3 