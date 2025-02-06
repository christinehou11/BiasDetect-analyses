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
load(here("processed-data","spatialHPC_SRT","spe_sub4.rda"))
clusters1_sub4 <- read.csv(
    here("processed-data","spatialHPC_SRT","seuInt_sub4-hpc_k-7_svgs_metadata.csv"), 
    row.names=1)
clusters2_sub4 <- read.csv(
    here("processed-data","spatialHPC_SRT",
        "seuInt_sub4-hpc_k-7_svgs-no-bias_metadata.csv"), 
    row.names=1)

# ---------
# add PRECAST results to colData
# ---------
identical(rownames(clusters1_sub4), rownames(colData(spe_sub4)))
spe_sub4$precast_k7 <- clusters1_sub4$cluster
spe_sub4$precast_k7_ordered <- factor(spe_sub4$precast_k7, levels=c(7,2,3,1,6,5,4), 
    labels=c("WM","WM (2)","SR/SL","CA1","CA3","DG GCL","DG ML"))

identical(rownames(clusters2_sub4), rownames(colData(spe_sub4)))
spe_sub4$precast_k7_nobias <- clusters2_sub4$cluster
spe_sub4$precast_k7_nobias_ordered <- factor(spe_sub4$precast_k7_nobias, 
    levels=c(1,2,7,5,6,4,3),
    labels=c("WM","SR/SL","CA1","CA1 (2)","CA3","DG GCL","DG ML"))

# ---------
# heatmap to justify cluster annotations
# ---------
spe_sub4 <- logNormCounts(spe_sub4)
markers <- c("MBP","GFAP","SPARCL1","FIBCD1",
            "COL5A2","KCNQ5","CARTPT","PCDH8","CALB1")

png(here("plots", "spatialHPC_SRT","cluster_sub4_heatmap_with_bias.png"), 
    width=5, height=5, units="in", res=300)
p1_heat <- plotGroupedHeatmap(spe_sub4, features = markers, 
                    swap_rownames="gene_name", 
                    group="precast_k7_ordered",
                    scale=TRUE, center=TRUE, 
                    cluster_rows=FALSE, cluster_cols=FALSE)
p1_heat
dev.off()

png(here("plots", "spatialHPC_SRT","cluster_sub4_heatmap_no_bias.png"), 
    width=5, height=5, units="in", res=300)
p2_heat <- plotGroupedHeatmap(spe_sub4, features = markers, 
                    swap_rownames="gene_name",
                    group="precast_k7_nobias_ordered",
                    scale=TRUE, center=TRUE, 
                    cluster_rows=FALSE, cluster_cols=FALSE)
p2_heat
dev.off()

# ---------
# plot cluster results
# ---------
l2_sub4 = unique(spe_sub4$sample_id)
names(l2_sub4) = l2_sub4
l2_sub4 = lapply(l2_sub4, function(x) spe_sub4[,colData(spe_sub4)$sample_id==x])

col.pal1 = c("#1f77b4FF","#aec7e8FF","#ffbb78FF",
            "#2ca02cFF","#ff7f0eFF","#d62728FF","#ff9896FF")
col.pal2 = c("#1f77b4FF","#ffbb78FF","#2ca02cFF",
            "#98df8aFF","#ff7f0eFF","#d62728FF","#ff9896FF")

png(here("plots", "spatialHPC_SRT","precast_sub4_cluster_bias.png"), 
    width=5, height=5, units="in", res=300)
c1 <- lapply(seq_along(l2_sub4), function(x) {
  plotSpots(l2_sub4[[x]], annotate="precast_k7_ordered", point_size=.3)+
    labs(color="clus")+
    scale_color_manual(values=col.pal1)+
    theme(plot.title=element_text(size=8))
})
drawFigs(c1, layout.dim = c(1, 4), common.legend = TRUE, 
        legend.position = "right", align = "h")
dev.off()

png(here("plots", "spatialHPC_SRT","precast_sub4_cluster_no_bias.png"), 
    width=5, height=5, units="in", res=300)
c2 <- lapply(seq_along(l2_sub4), function(x) {
  plotSpots(l2_sub4[[x]], annotate="precast_k7_nobias_ordered", point_size=.3)+
    labs(color="clus")+
    scale_color_manual(values=col.pal2)+
    theme(plot.title=element_text(size=8))
})
drawFigs(c2, layout.dim = c(1, 4), common.legend = TRUE, 
        legend.position = "right", align = "h")
dev.off()
