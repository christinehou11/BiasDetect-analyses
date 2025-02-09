suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(scater)
    library(ggspavis)
    library(scuttle)
    library(PRECAST)
    library(purrr)
})

# ---------
# load data
# ---------
load(here("processed-data","spatialHPC_SRT","spe_sub2.rda"))
clusters1_sub2 <- read.csv(
    here("processed-data","spatialHPC_SRT","seuInt_sub2-hpc_k-7_svgs_metadata.csv"), 
    row.names=1)
clusters2_sub2 <- read.csv(
    here("processed-data","spatialHPC_SRT",
        "seuInt_sub2-hpc_k-7_svgs-no-bias_metadata.csv"), 
    row.names=1)

# Label
cluster1_domain_map <- map_dfr(srt.sets.sub2, ~ .x@meta.data) %>% select(domain)
cluster1_domain_map <- merge(cluster1_domain_map, clusters1_sub2, by = "row.names", all = TRUE)
rownames(cluster1_domain_map) <- cluster1_domain_map$Row.names
cluster1_domain_map <- cluster1_domain_map %>% select(domain, cluster)

cluster1_domain_map %>% 
  group_by(cluster, domain) %>% 
  summarize(n = n(), .groups = "drop") %>% 
  group_by(cluster) %>% 
  filter(n == max(n))
# A tibble: 7 × 3
# Groups:   cluster [7]
# cluster domain     n
# <int> <fct>   <int>
#      1 ML        799
#      2 Choroid    94
#      3 CA2.4     884
#      4 GCL       502
#      5 SR.SLM   1469
#      6 CA1      1499
#      7 WM.1      781

cluster2_domain_map <- map_dfr(srt.sets.sub2, ~ .x@meta.data) %>% select(domain)
cluster2_domain_map <- merge(cluster2_domain_map, clusters2_sub2, by = "row.names", all = TRUE)
rownames(cluster2_domain_map) <- cluster2_domain_map$Row.names
cluster2_domain_map <- cluster2_domain_map %>% select(domain, cluster)

cluster2_domain_map %>% 
  group_by(cluster, domain) %>% 
  summarize(n = n(), .groups = "drop") %>% 
  group_by(cluster) %>% 
  filter(n == max(n))
# A tibble: 7 × 3
# Groups:   cluster [7]
# cluster domain   n
# <int> <fct>  <int>
#     1 ML       814
#     2 GABA     405
#     3 WM.1     752
#     4 GCL      500
#     5 CA1     1467
#     6 SR.SLM  1478
#     7 CA2.4   1047

# ---------
# add PRECAST results to colData
# ---------
identical(rownames(clusters1_sub2), rownames(colData(spe_sub2)))
spe_sub2$precast_k7 <- clusters1_sub2$cluster
spe_sub2$precast_k7_ordered <- factor(spe_sub2$precast_k7, levels=c(7,1,4,5,6,3,2), 
    labels=c('WM','DG ML','DG GCL',"SR/SL",'CA1','CA2','Choroid'))

identical(rownames(clusters2_sub2), rownames(colData(spe_sub2)))
spe_sub2$precast_k7_nobias <- clusters2_sub2$cluster
spe_sub2$precast_k7_nobias_ordered <- factor(spe_sub2$precast_k7_nobias, 
    levels=c(3,1,4,6,5,7,2),
    labels=c("WM",'DG ML','DG GCL',"SR/SL",'CA1','CA2','GABA'))

# ---------
# heatmap to justify cluster annotations
# ---------
spe_sub2 <- logNormCounts(spe_sub2)
markers <- c("MBP","GFAP","SPARCL1","FIBCD1",
            "COL5A2","KCNQ5","CARTPT","PCDH8","CALB1")

png(here("plots", "spatialHPC_SRT","cluster_sub2_heatmap_with_bias.png"), 
    width=5, height=5, units="in", res=300)
p1_heat <- plotGroupedHeatmap(spe_sub2, features = markers, 
                    swap_rownames="gene_name", 
                    group="precast_k7_ordered",
                    scale=TRUE, center=TRUE, 
                    cluster_rows=FALSE, cluster_cols=FALSE)
p1_heat
dev.off()

png(here("plots", "spatialHPC_SRT","cluster_sub2_heatmap_no_bias.png"), 
    width=5, height=5, units="in", res=300)
p2_heat <- plotGroupedHeatmap(spe_sub2, features = markers, 
                    swap_rownames="gene_name",
                    group="precast_k7_nobias_ordered",
                    scale=TRUE, center=TRUE, 
                    cluster_rows=FALSE, cluster_cols=FALSE)
p2_heat
dev.off()

# ---------
# plot cluster results
# ---------
l2_sub2 = unique(spe_sub2$sample_id)
names(l2_sub2) = l2_sub2
l2_sub2 = lapply(l2_sub2, function(x) spe_sub2[,colData(spe_sub2)$sample_id==x])

# CA1: 2ca02cFF
# CA2: 98df8aFF
# GABA: ff7f0eFF
# Choroid: AEC7E8FF
# DG ML: C5B0D5FF
# DG GCL: 17BECFFF
# SR/SL: d62728FF
# WM: 1f77b4FF
col.pal1 = c("#1f77b4FF", "#C5B0D5FF","#17BECFFF", 
             "#d62728FF","#2ca02cFF", "#98df8aFF", "#AEC7E8FF")
col.pal2 = c("#1f77b4FF", "#C5B0D5FF","#17BECFFF", 
             "#d62728FF","#2ca02cFF", "#98df8aFF", "#ff7f0eFF")

png(here("plots", "spatialHPC_SRT","precast_sub2_cluster_bias.png"), 
    width=5, height=5, units="in", res=300)
c1 <- lapply(seq_along(l2_sub2), function(x) {
  plotSpots(l2_sub2[[x]], annotate="precast_k7_ordered", point_size=.3)+
    labs(color="clus")+
    scale_color_manual(values=col.pal1)+
    theme(plot.title=element_text(size=8))
})
drawFigs(c1, layout.dim = c(1, 2), common.legend = TRUE, 
        legend.position = "right", align = "h")
dev.off()

png(here("plots", "spatialHPC_SRT","precast_sub2_cluster_no_bias.png"), 
    width=5, height=5, units="in", res=300)
c2 <- lapply(seq_along(l2_sub2), function(x) {
  plotSpots(l2_sub2[[x]], annotate="precast_k7_nobias_ordered", point_size=.3)+
    labs(color="clus")+
    scale_color_manual(values=col.pal2)+
    theme(plot.title=element_text(size=8))
})
drawFigs(c2, layout.dim = c(1, 2), common.legend = TRUE, 
        legend.position = "right", align = "h")
dev.off()
