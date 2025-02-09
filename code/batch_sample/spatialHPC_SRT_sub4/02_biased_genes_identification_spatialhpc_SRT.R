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
#            gene gene_name dev_default rank_default dev_batch rank_batch
# ENSG00000002330       BAD    15830.43         1310  15785.12       1289
# ENSG00000002586      CD99    18837.21          865  18616.10        858
# ENSG00000004059      ARF5    15033.32         1443  14985.57       1417
# ENSG00000004660    CAMKK1    14586.54         1506  14479.27       1487
# ENSG00000004779   NDUFAB1    19859.78          756  19469.72        764
# ENSG00000005022   SLC25A5    16842.03         1150  16739.61       1131

# d_diff     nSD_dev r_diff   nSD_rank
# 0.002870288 -0.31052071    -21 -0.4914075
# 0.011877487 -0.13869673     -7 -0.1638025
# 0.003186151 -0.30449522    -26 -0.6084093
# 0.007408683 -0.22394495    -19 -0.4446068
# 0.020034113  0.01690146      8  0.1872029
# 0.006118529 -0.24855631    -19 -0.4446068

# list bias as data frame
bias <- BiasDetect::biasDetect(batch_df = batch_df, nSD_dev = 10, nSD_rank = 5)

# display bias in plots
bias_plots <- BiasDetect::biasDetect(batch_df = batch_df, 
                        nSD_dev = 10, nSD_rank = 5, visual = TRUE)
# deviance
png(here("plots", "spatialHPC_SRT","biased_iden_dev_sub4.png"), 
      width=5, height=5, units="in", res=300)
p1_dev <- bias_plots$deviance
p1_dev
dev.off()

# rank
png(here("plots", "spatialHPC_SRT","biased_iden_rank_sub4.png"), 
    width=5, height=5, units="in", res=300)
p2_rank <- bias_plots$rank
p2_rank
dev.off()

# -----------
# save object
# -----------
write.csv(batch_df, 
    here("processed-data","spatialHPC_SRT",
        "hpc_bindev_default-sample_svgs-only_sub4.csv"), 
    row.names=FALSE)

save(bias,
    file = here("processed-data","spatialHPC_SRT","hpc_biased_features_sub4.rda"))

