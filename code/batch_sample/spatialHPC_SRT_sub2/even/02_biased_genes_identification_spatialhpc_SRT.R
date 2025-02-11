suppressPackageStartupMessages({
    library(BiasDetect)
    library(dplyr)
    library(here)
    library(SummarizedExperiment)
})

# ---------
# load data
# ---------
load(here("processed-data","spatialHPC_SRT","spe_sub2.rda"))

# ---------
# biased gene identification
# ---------
SVGs <- rowData(spe_sub2)$gene_id
batch_df_sub2 <- BiasDetect::featureSelect(spe_sub2, 
                                    batch_effect = "sample_id", VGs = SVGs)
head(batch_df_sub2)
#            gene gene_name dev_default rank_default dev_batch rank_batch      d_diff 
# ENSG00000078369      GNB1   12534.709          503 12489.860        494 0.003590815
# ENSG00000187730     GABRD    8220.839         1201  8035.515       1216 0.023063058
# ENSG00000067606     PRKCZ    9847.529          909  9816.229        901 0.003188580
# ENSG00000158109    TPRG1L    8925.409         1072  8912.729       1053 0.001422636
# ENSG00000069424    KCNAB2   11927.047          578 11885.417        567 0.003502675 
# ENSG00000116254      CHD5    7778.070         1268  7724.372       1261 0.006951753
#    nSD_dev r_diff   nSD_rank
# -0.4554443     -9 -0.6683004
#  0.9450326     15  1.1138339
# -0.4843737     -8 -0.5940448
# -0.6113833    -19 -1.4108563
# -0.4617835    -11 -0.8168115
# -0.2137199     -7 -0.5197892

# Figure - nSD_dev
sd.interval = 4
batch_df_sub2$nSD.bin_dev = cut(abs(batch_df_sub2$nSD_dev), right=FALSE,
                                breaks=seq(0,max(batch_df_sub2$nSD_dev)+sd.interval, 
                                           by=sd.interval),
                                include.lowest=TRUE)

col.pal = RColorBrewer::brewer.pal(length(unique(batch_df_sub2$nSD.bin_dev)), "YlOrRd")
col.pal[1] = "grey"

png(here("plots", "spatialHPC_SRT","SVGs_nSD_dev_sub2.png"))
dev3 = ggplot(batch_df_sub2, aes(x=d_diff, fill=nSD.bin_dev))+
  geom_histogram(color="grey20", bins=50)+
  scale_fill_manual(values=col.pal)+
  labs(x="\u0394 deviance",
       fill="n abs(SD)", y="# SVGs", title="nSD bin width = 4")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=10^(0:4), labels=format(10^(0:4), scientific=F))+
  theme_bw()+
  theme(legend.position="inside",legend.position.inside=c(.8,.7),
        aspect.ratio=1)
dev3
dev.off()

# Figure - nSD_rank
sd.interval = 4
batch_df_sub2$nSD.bin_rank = cut(abs(batch_df_sub2$nSD_rank), right=FALSE,
                                 breaks=seq(0,max(batch_df_sub2$nSD_rank)+sd.interval, 
                                            by=sd.interval),
                                 include.lowest=TRUE)

col.pal3 = RColorBrewer::brewer.pal(length(unique(batch_df_sub2$nSD.bin_rank)), "YlOrRd")
col.pal3[1] = "grey"

png(here("plots", "spatialHPC_SRT","SVGs_nSD_rank_sub2.png"))
rank4 <- ggplot(batch_df_sub2, aes(x=r_diff, fill=nSD.bin_rank))+
  geom_histogram(color="grey20", bins=30)+
  scale_fill_manual(values=col.pal3)+
  labs(x="rank difference", fill="n abs(SD)", y="# genes", title="nSD bin width = 4")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=10^(0:4), labels=format(10^(0:4), scientific=F))+
  theme_bw()+theme(legend.position="inside",legend.position.inside=c(.8,.7))

rank4
dev.off()

# list bias as data frame
bias <- BiasDetect::biasDetect(batch_df = batch_df_sub2, nSD_dev = 4, nSD_rank = 4)

# display bias in plots
bias_plots <- BiasDetect::biasDetect(batch_df = batch_df_sub2, 
                        nSD_dev = 4, nSD_rank = 4, visual = TRUE)
# deviance
png(here("plots", "spatialHPC_SRT","biased_iden_dev_sub2.png"), 
      width=5, height=5, units="in", res=300)
p1_dev <- bias_plots$deviance
p1_dev
dev.off()

# rank
png(here("plots", "spatialHPC_SRT","biased_iden_rank_sub2.png"), 
    width=5, height=5, units="in", res=300)
p2_rank <- bias_plots$rank
p2_rank
dev.off()

# -----------
# save object
# -----------
write.csv(batch_df_sub2, 
    here("processed-data","spatialHPC_SRT",
        "hpc_bindev_default-sample_svgs-only_sub2.csv"), 
    row.names=FALSE)

save(bias,
    file = here("processed-data","spatialHPC_SRT","hpc_biased_features_sub2.rda"))

