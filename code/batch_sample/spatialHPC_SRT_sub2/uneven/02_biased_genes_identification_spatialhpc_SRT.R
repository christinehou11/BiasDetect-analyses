suppressPackageStartupMessages({
    library(BiasDetect)
    library(dplyr)
    library(here)
    library(SummarizedExperiment)
})

# ---------
# load data
# ---------
load(here("processed-data","spatialHPC_SRT","spe_sub2_2.rda"))

# ---------
# biased gene identification
# ---------
SVGs <- rowData(spe_sub2_2)$gene_id
batch_df_sub2_2 <- BiasDetect::featureSelect(spe_sub2_2, 
                                    batch_effect = "sample_id", VGs = SVGs)
head(batch_df_sub2_2)
#            gene gene_name dev_default rank_default dev_batch
# ENSG00000078369      GNB1    6068.537          531  6056.118
# ENSG00000187730     GABRD    4119.284         1189  4061.836
# ENSG00000067606     PRKCZ    4900.011          909  4874.738
# ENSG00000158109    TPRG1L    4336.524         1110  4321.293
# ENSG00000069424    KCNAB2    5774.298          615  5757.917
# ENSG00000116254      CHD5    3762.140         1292  3751.480
# rank_batch      d_diff    nSD_dev r_diff   nSD_rank
#        523 0.002050745 -0.4284529     -8 -0.7066865
#       1197 0.014143304  0.3600388      8  0.7066865
#        904 0.005184634 -0.2241086     -5 -0.4416790
#       1098 0.003524823 -0.3323360    -12 -1.0600297
#        602 0.002844943 -0.3766674    -13 -1.1483655
#       1288 0.002841432 -0.3768963     -4 -0.3533432

# Figure - nSD_dev
sd.interval = 4
batch_df_sub2_2$nSD.bin_dev = cut(abs(batch_df_sub2_2$nSD_dev), right=FALSE,
                                breaks=seq(0,max(batch_df_sub2_2$nSD_dev)+sd.interval, 
                                           by=sd.interval),
                                include.lowest=TRUE)

col.pal = RColorBrewer::brewer.pal(length(unique(batch_df_sub2_2$nSD.bin_dev)), "YlOrRd")
col.pal[1] = "grey"

png(here("plots", "spatialHPC_SRT","SVGs_nSD_dev_sub2_2.png"))
dev3 = ggplot(batch_df_sub2_2, aes(x=d_diff, fill=nSD.bin_dev))+
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
sd.interval = 5
batch_df_sub2_2$nSD.bin_rank = cut(abs(batch_df_sub2_2$nSD_rank), right=FALSE,
                                 breaks=seq(0,max(batch_df_sub2_2$nSD_rank)+sd.interval, 
                                            by=sd.interval),
                                 include.lowest=TRUE)

col.pal3 = RColorBrewer::brewer.pal(length(unique(batch_df_sub2_2$nSD.bin_rank)), "YlOrRd")
col.pal3[1] = "grey"

png(here("plots", "spatialHPC_SRT","SVGs_nSD_rank_sub2_2.png"))
rank4 <- ggplot(batch_df_sub2_2, aes(x=r_diff, fill=nSD.bin_rank))+
  geom_histogram(color="grey20", bins=30)+
  scale_fill_manual(values=col.pal3)+
  labs(x="rank difference", fill="n abs(SD)", y="# genes", title="nSD bin width = 5")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=10^(0:4), labels=format(10^(0:4), scientific=F))+
  theme_bw()+theme(legend.position="inside",legend.position.inside=c(.8,.7))

rank4
dev.off()

# list bias as data frame
bias_sub2_2 <- BiasDetect::biasDetect(batch_df = batch_df_sub2_2, nSD_dev = 4, nSD_rank = 5)

# display bias in plots
bias_plots <- BiasDetect::biasDetect(batch_df = batch_df_sub2_2, 
                        nSD_dev = 4, nSD_rank = 4, visual = TRUE)
# deviance
png(here("plots", "spatialHPC_SRT","biased_iden_dev_sub2_2.png"), 
      width=5, height=5, units="in", res=300)
p1_dev <- bias_plots$deviance
p1_dev
dev.off()

# rank
png(here("plots", "spatialHPC_SRT","biased_iden_rank_sub2_2.png"), 
    width=5, height=5, units="in", res=300)
p2_rank <- bias_plots$rank
p2_rank
dev.off()

# -----------
# save object
# -----------
write.csv(batch_df_sub2_2, 
    here("processed-data","spatialHPC_SRT",
        "hpc_bindev_default-sample_svgs-only_sub2_2.csv"), 
    row.names=FALSE)

save(bias_sub2_2,
    file = here("processed-data","spatialHPC_SRT","hpc_biased_features_sub2_2.rda"))

