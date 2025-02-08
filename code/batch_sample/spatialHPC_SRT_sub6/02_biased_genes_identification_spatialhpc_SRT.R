suppressPackageStartupMessages({
    library(BiasDetect)
    library(dplyr)
    library(here)
    library(SummarizedExperiment)
})

# ---------
# load data
# ---------
load(here("processed-data","spatialHPC_SRT","spe_sub6.rda"))

# ---------
# biased gene identification
# ---------
SVGs_sub6 <- rowData(spe_sub6)$gene_id
length(SVGs_sub6) # 2333

batch_df_sub6 <- BiasDetect::featureSelect(spe_sub6, 
            batch_effect = "sample_id", VGs = SVGs_sub6)
head(batch_df_sub6)
# gene gene_name dev_default rank_default dev_batch rank_batch      d_diff    nSD_dev
# 1 ENSG00000131584     ACAP3    15248.39         1549  15123.92       1515 0.008230092 -0.6120571
# 2 ENSG00000175756  AURKAIP1    19044.93         1024  19017.70        958 0.001432143 -0.7366072
# 3 ENSG00000242485    MRPL20    19666.65          951  19518.14        895 0.007608703 -0.6234420
# 4 ENSG00000179403      VWA1    12565.06         1816  12297.94       1811 0.021720397 -0.3648917
# 5 ENSG00000160075     SSU72    17752.21         1229  17578.88       1183 0.009860182 -0.5821911
# 6 ENSG00000078369      GNB1    24903.75          504  24073.06        499 0.034507129 -0.1306169
# r_diff   nSD_rank
# 1    -34 -0.7642247
# 2    -66 -1.4834950
# 3    -56 -1.2587230
# 4     -5 -0.1123860
# 5    -46 -1.0339511
# 6     -5 -0.1123860

# Figure
sd.interval = 5
batch_df_sub6$nSD.bin_dev = cut(abs(batch_df_sub6$nSD_dev), right=FALSE,
                           breaks=seq(0,max(batch_df_sub6$nSD_dev)+sd.interval, 
                                      by=sd.interval),
                           include.lowest=TRUE)

col.pal = RColorBrewer::brewer.pal(length(unique(batch_df_sub6$nSD.bin_dev)), "YlOrRd")
col.pal[1] = "grey"

png(here("plots", "spatialHPC_SRT","SVGs_nSD_dev_sub6.png"), 
    width=5, height=5, units="in", res=300)
dev3 = ggplot(batch_df_sub6, aes(x=d_diff, fill=nSD.bin_dev))+
  geom_histogram(color="grey20", bins=50)+
  scale_fill_manual(values=col.pal)+
  labs(x="\u0394 deviance",
       fill="n abs(SD)", y="# SVGs")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=10^(0:4), labels=format(10^(0:4), scientific=F))+
  theme_bw()+
  theme(legend.position="inside",legend.position.inside=c(.8,.7),
        aspect.ratio=1)
dev.off()


# list bias as data frame
bias_sub6 <- BiasDetect::biasDetect(batch_df = batch_df_sub6, nSD_dev = 5, nSD_rank = 3)

# display bias in plots
bias_plots_sub6 <- BiasDetect::biasDetect(batch_df = batch_df_sub6, 
                        nSD_dev = 10, nSD_rank = 5, visual = TRUE)
# deviance
png(here("plots", "spatialHPC_SRT","biased_iden_dev_sub6.png"), 
      width=5, height=5, units="in", res=300)
p1_dev <- bias_plots_sub6$deviance
p1_dev
dev.off()

# rank
png(here("plots", "spatialHPC_SRT","biased_iden_rank_sub6.png"), 
    width=5, height=5, units="in", res=300)
p2_rank <- bias_plots_sub6$rank
p2_rank
dev.off()

# -----------
# save object
# -----------
write.csv(batch_df_sub6, 
    here("processed-data","spatialHPC_SRT",
        "hpc_bindev_default-sample_svgs-only_sub6.csv"), 
    row.names=FALSE)

save(bias_sub6,
    file = here("processed-data","spatialHPC_SRT","hpc_biased_features_sub6.rda"))
