suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(scater)
    library(ggspavis)
    library(scuttle)
    library(PRECAST)
})

# ---------
# load data
# ---------
load(here("processed-data","spatialHPC_SRT","spe_sub6.rda"))
clusters1_sub6 <- read.csv(
    here("processed-data","spatialHPC_SRT","seuInt_sub6-hpc_k-7_svgs_metadata.csv"), 
    row.names=1)
clusters2_sub6 <- read.csv(
    here("processed-data","spatialHPC_SRT",
        "seuInt_sub6-hpc_k-7_svgs-no-bias_metadata.csv"), 
    row.names=1)

# ---------
# add PRECAST results to colData
# ---------
identical(rownames(clusters1_sub6), rownames(colData(spe_sub6)))
spe_sub6$precast_k7 <- clusters1_sub6$cluster
spe_sub6$precast_k7_ordered <- factor(spe_sub6$precast_k7, levels=c(7,2,3,1,6,5,4), 
    labels=c("WM","WM (2)","SR/SL","CA1","CA3","DG GCL","DG ML"))

identical(rownames(clusters2_sub6), rownames(colData(spe_sub6)))
spe_sub6$precast_k7_nobias <- clusters2$cluster
spe_sub6$precast_k7_nobias_ordered <- factor(spe_sub6$precast_k7_nobias, 
    levels=c(1,2,7,5,6,4,3),
    labels=c("WM","SR/SL","CA1","CA1 (2)","CA3","DG GCL","DG ML"))

# ---------
# heatmap to justify cluster annotations
# ---------
spe_sub6 <- logNormCounts(spe_sub6)
markers <- c("MBP","GFAP","SPARCL1","FIBCD1",
            "COL5A2","KCNQ5","CARTPT","PCDH8","CALB1")

png(here("plots", "spatialHPC_SRT","cluster_heatmap_with_bias_sub6.png"), 
    width=5, height=5, units="in", res=300)
p1_heat <- plotGroupedHeatmap(spe_sub6, features = markers, 
                    swap_rownames="gene_name", 
                    group="precast_k7_ordered",
                    scale=TRUE, center=TRUE, 
                    cluster_rows=FALSE, cluster_cols=FALSE)
p1_heat
dev.off()

png(here("plots", "spatialHPC_SRT","cluster_heatmap_no_bias.png"), 
    width=5, height=5, units="in", res=300)
p2_heat <- plotGroupedHeatmap(spe, features = markers, 
                    swap_rownames="gene_name",
                    group="precast_k7_nobias_ordered",
                    scale=TRUE, center=TRUE, 
                    cluster_rows=FALSE, cluster_cols=FALSE)
p2_heat
dev.off()

# ---------
# plot cluster results
# ---------
l2_sub6 = unique(spe_sub6$sample_id)
names(l2_sub6) = l2_sub6
l2_sub6 = lapply(l2_sub6, function(x) spe_sub6[,colData(spe_sub6)$sample_id==x])

col.pal1 = c("#1f77b4FF","#aec7e8FF","#ffbb78FF",
            "#2ca02cFF","#ff7f0eFF","#d62728FF","#ff9896FF")
col.pal2 = c("#1f77b4FF","#ffbb78FF","#2ca02cFF",
            "#98df8aFF","#ff7f0eFF","#d62728FF","#ff9896FF")

png(here("plots", "spatialHPC_SRT","precast_cluster_bias_sub6.png"), 
    width=5, height=5, units="in", res=300)
c1 <- lapply(seq_along(l2_sub6), function(x) {
  plotSpots(l2_sub6[[x]], annotate="precast_k7_ordered", point_size=.3)+
    labs(color="clus")+
    scale_color_manual(values=col.pal1)+
    theme(plot.title=element_text(size=8))
})
drawFigs(c1, layout.dim = c(2,3), common.legend = TRUE, 
        legend.position = "right", align = "h")
dev.off()

png(here("plots", "spatialHPC_SRT","precast_cluster_no_bias.png"), 
    width=5, height=5, units="in", res=300)
c2 <- lapply(seq_along(l2_sub6), function(x) {
  plotSpots(l2_sub6[[x]], annotate="precast_k7_nobias_ordered", point_size=.3)+
    labs(color="clus")+
    scale_color_manual(values=col.pal2)+
    theme(plot.title=element_text(size=8))
})
drawFigs(c2, layout.dim = c(2, 3), common.legend = TRUE, 
        legend.position = "right", align = "h")
dev.off()

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
#   [1] ggspavis_1.11.0             scater_1.33.4              
# [3] ggplot2_3.5.1               scuttle_1.15.4             
# [5] BiasDetect_0.99.0           PRECAST_1.6.5              
# [7] gtools_3.9.5                Seurat_5.1.0               
# [9] SeuratObject_5.0.2          sp_2.1-4                   
# [11] tibble_3.2.1                tidyr_1.3.1                
# [13] dplyr_1.1.4                 SpatialExperiment_1.15.1   
# [15] SingleCellExperiment_1.27.2 SummarizedExperiment_1.35.4
# [17] Biobase_2.65.1              GenomicRanges_1.57.2       
# [19] GenomeInfoDb_1.41.2         IRanges_2.39.2             
# [21] S4Vectors_0.43.2            BiocGenerics_0.51.3        
# [23] MatrixGenerics_1.17.0       matrixStats_1.4.1          
# [25] here_1.0.1                 
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.4.1           later_1.3.2            
# [4] bitops_1.0-9            polyclip_1.10-7         fastDummies_1.7.4      
# [7] lifecycle_1.0.4         rstatix_0.7.2           rprojroot_2.0.4        
# [10] processx_3.8.4          globals_0.16.3          lattice_0.22-6         
# [13] MASS_7.3-61             backports_1.5.0         magrittr_2.0.3         
# [16] plotly_4.10.4           rmarkdown_2.28          remotes_2.5.0          
# [19] yaml_2.3.10             httpuv_1.6.15           sctransform_0.4.1      
# [22] ggside_0.3.1            spam_2.11-0             pkgbuild_1.4.4         
# [25] spatstat.sparse_3.1-0   reticulate_1.39.0       cowplot_1.1.3          
# [28] pbapply_1.7-2           RColorBrewer_1.1-3      abind_1.4-8            
# [31] zlibbioc_1.51.1         Rtsne_0.17              purrr_1.0.2            
# [34] RCurl_1.98-1.16         GenomeInfoDbData_1.2.13 ggrepel_0.9.6          
# [37] scry_1.17.0             irlba_2.3.5.1           listenv_0.9.1          
# [40] spatstat.utils_3.1-0    pheatmap_1.0.12         goftest_1.2-3          
# [43] RSpectra_0.16-2         spatstat.random_3.3-2   fitdistrplus_1.2-1     
# [46] parallelly_1.38.0       leiden_0.4.3.1          codetools_0.2-20       
# [49] DelayedArray_0.31.14    tidyselect_1.2.1        UCSC.utils_1.1.0       
# [52] farver_2.1.2            viridis_0.6.5           ScaledMatrix_1.13.0    
# [55] spatstat.explore_3.3-2  jsonlite_1.8.9          BiocNeighbors_1.99.2   
# [58] Formula_1.2-5           progressr_0.14.0        ggridges_0.5.6         
# [61] survival_3.7-0          tools_4.4.1             ica_1.0-3              
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
# [151] shiny_1.9.1             ROCR_1.0-11             broom_1.0.7            
# [154] igraph_2.0.3 