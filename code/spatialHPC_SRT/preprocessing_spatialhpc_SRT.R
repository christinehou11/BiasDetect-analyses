library(here)
library(SpatialExperiment)
library(dplyr)
library(tidyr)
library(tibble)

# ---------
# load data
# ---------

load(here("raw-data","spatialHPC_SRT","spe_precast_HE_domain.rda"))

fix_order = distinct(as.data.frame(colData(spe)), 
                slide, array, brnum, sample_id, position, sex) %>% 
            arrange(slide, array)
sub4 = fix_order$sample_id[c(14,16, 20,21)]
spe_sub4 = spe[,spe$sample_id %in% sub4]

#load and reformat nnSVG results
load(here("raw-data","spatialHPC_SRT","nnSVG_outs_HE_only.rda"))

res_df = pivot_longer(
            rownames_to_column(as.data.frame(res_ranks), var="gene_id"), 
            colnames(res_ranks), names_to="sample_id", values_to="rank", 
            values_drop_na=T)

# filter to only the top 2k sig features in the 4 samples we're using
res_df2 = filter(res_df, 
            sample_id %in% c("V11L05-333_B1","V11L05-333_D1",
                            "V11L05-335_D1","V11L05-336_A1"),
            rank<=2000)

# further filter to only features that are in the top 2k of >1 sample
svgs = group_by(res_df2, gene_id) %>% 
        tally() %>% 
        filter(n>1) 
nrow(svgs) # 2082

spe_sub4 = spe_sub4[svgs$gene_id,]
dim(spe_sub4) # 2082 18945

# -----------
# save object
# -----------

write.csv(as.data.frame(as.matrix(counts(spe_sub4)[1:1000,])), 
        here("processed-data","spatialHPC_SRT",
            "spe-hpc_sub4_svgs-only_counts-1.csv"),
        row.names=T)

write.csv(as.data.frame(as.matrix(counts(spe_sub4)[1001:2082,])), 
      here("processed-data","spatialHPC_SRT",
          "spe-hpc_sub4_svgs-only_counts-2.csv"),
      row.names=T)

write.csv(colData(spe_sub4)[,c("sample_id","slide","brnum","array",
                              "sex","in_tissue","array_row","array_col",
                              "sum_umi","sum_gene","expr_chrM",
                              "expr_chrM_ratio")],
      here("processed-data","spatialHPC_SRT",
          "spe-hpc_sub4_svgs-only_colData.csv"),
      row.names=T)

write.csv(rowData(spe_sub4), 
      here("processed-data","spatialHPC_SRT","spe-hpc_sub4_svgs-only_rowData.csv"),
      row.names=T)

write.csv(spatialCoords(spe_sub4), 
      here("processed-data","spatialHPC_SRT",
          "spe-hpc_sub4_svgs-only_spatialCoords.csv"),
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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tibble_3.2.1                tidyr_1.3.1                 dplyr_1.1.4                
# [4] SpatialExperiment_1.15.1    SingleCellExperiment_1.27.2 SummarizedExperiment_1.35.4
# [7] Biobase_2.65.1              GenomicRanges_1.57.2        GenomeInfoDb_1.41.2        
# [10] IRanges_2.39.2              S4Vectors_0.43.2            BiocGenerics_0.51.3        
# [13] MatrixGenerics_1.17.0       matrixStats_1.4.1           here_1.0.1                 
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.4              generics_0.1.3          SparseArray_1.5.44     
# [4] lattice_0.22-6          digest_0.6.37           magrittr_2.0.3         
# [7] evaluate_1.0.1          grid_4.4.1              fastmap_1.2.0          
# [10] rprojroot_2.0.4         jsonlite_1.8.9          Matrix_1.7-0           
# [13] httr_1.4.7              purrr_1.0.2             fansi_1.0.6            
# [16] UCSC.utils_1.1.0        abind_1.4-8             cli_3.6.3              
# [19] rlang_1.1.4             crayon_1.5.3            XVector_0.45.0         
# [22] withr_3.0.1             DelayedArray_0.31.14    yaml_2.3.10            
# [25] S4Arrays_1.5.11         tools_4.4.1             GenomeInfoDbData_1.2.13
# [28] vctrs_0.6.5             R6_2.5.1                lifecycle_1.0.4        
# [31] magick_2.8.5            zlibbioc_1.51.1         pkgconfig_2.0.3        
# [34] pillar_1.9.0            glue_1.8.0              Rcpp_1.0.13            
# [37] xfun_0.48               tidyselect_1.2.1        rstudioapi_0.16.0      
# [40] knitr_1.48              rjson_0.2.23            htmltools_0.5.8.1      
# [43] rmarkdown_2.28          compiler_4.4.1  