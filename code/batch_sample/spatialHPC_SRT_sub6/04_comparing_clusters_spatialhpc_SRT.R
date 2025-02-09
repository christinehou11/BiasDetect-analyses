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
load(here("processed-data","spatialHPC_SRT","spe_sub6.rda"))
clusters1_sub6 <- read.csv(
    here("processed-data","spatialHPC_SRT","seuInt_sub6-hpc_k-7_svgs_metadata.csv"), 
    row.names=1)
clusters2_sub6 <- read.csv(
    here("processed-data","spatialHPC_SRT",
        "seuInt_sub6-hpc_k-7_svgs-no-bias_metadata.csv"), 
    row.names=1)

# Label
cluster1_domain_map <- map_dfr(srt.sets.sub6, ~ .x@meta.data) %>% select(domain)
cluster1_domain_map <- merge(cluster1_domain_map, clusters1_sub6, by = "row.names", all = TRUE)
rownames(cluster1_domain_map) <- cluster1_domain_map$Row.names
cluster1_domain_map <- cluster1_domain_map %>% select(domain, cluster)

cluster1_domain_map %>% 
  group_by(cluster, domain) %>% 
  summarize(n = n(), .groups = "drop") %>% 
  group_by(cluster) %>% 
  filter(n == max(n))
# A tibble: 6 × 3
# Groups:   cluster [6]
# cluster domain     n
# <int> <fct>  <int>
#     1 CA2.4   4374
#     2 SR.SLM  1838
#     3 ML      1299
#     4 GABA     957
#     5 SR.SLM  3030
#     6 CA1     2174

cluster2_domain_map <- map_dfr(srt.sets.sub6, ~ .x@meta.data) %>% select(domain)
cluster2_domain_map <- merge(cluster2_domain_map, clusters2_sub6, by = "row.names", all = TRUE)
rownames(cluster2_domain_map) <- cluster2_domain_map$Row.names
cluster2_domain_map <- cluster2_domain_map %>% select(domain, cluster)

cluster2_domain_map %>% 
  group_by(cluster, domain) %>% 
  summarize(n = n(), .groups = "drop") %>% 
  group_by(cluster) %>% 
  filter(n == max(n))
# A tibble: 6 × 3
# Groups:   cluster [6]
# cluster domain     n
# <int> <fct>  <int>
#     1 CA2.4   3328
#     2 WM.1    1603
#     3 SR.SLM  4062
#     4 GABA     689
#     5 CA1     2002
#     6 ML      1101

# ---------
# add PRECAST results to colData
# ---------
identical(rownames(clusters1_sub6), rownames(colData(spe_sub6)))
spe_sub6$precast_k6 <- clusters1_sub6$cluster
spe_sub6$precast_k6_ordered <- factor(spe_sub6$precast_k7, levels=c(6, 1,4,3,2,5), 
    labels=c("CA1", "CA2", "GABA","ML", "SR/SL", "SR/SL (2)"))

identical(rownames(clusters2_sub6), rownames(colData(spe_sub6)))
spe_sub6$precast_k6_nobias <- clusters2_sub6$cluster
spe_sub6$precast_k6_nobias_ordered <- factor(spe_sub6$precast_k7_nobias, 
    levels=c(2, 5, 1,4,6,3),
    labels=c("WM", "CA1", "CA2", "GABA",  "ML", "SR/SL"))

# ---------
# heatmap to justify cluster annotations
# ---------
spe_sub6 <- logNormCounts(spe_sub6)
markers <- c("MBP","GFAP","SPARCL1","FIBCD1",
            "COL5A2","KCNQ5","CARTPT","PCDH8","CALB1")
  # genes of interest for visualization in the heatmap

png(here("plots", "spatialHPC_SRT","cluster_heatmap_with_bias_sub6.png"), 
    width=5, height=5, units="in", res=300)
p1_heat <- plotGroupedHeatmap(spe_sub6, features = markers, 
                    swap_rownames="gene_name", 
                    group="precast_k6_ordered",
                    scale=TRUE, center=TRUE, 
                    cluster_rows=FALSE, cluster_cols=FALSE)
p1_heat
dev.off()

png(here("plots", "spatialHPC_SRT","cluster_heatmap_no_bias_sub6.png"), 
    width=5, height=5, units="in", res=300)
p2_heat <- plotGroupedHeatmap(spe_sub6, features = markers, 
                    swap_rownames="gene_name",
                    group="precast_k6_nobias_ordered",
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

# CA1: 2ca02cFF
# CA2: 98df8aFF
# GABA: ff7f0eFF
# ML: C5B0D5FF
# SR/SL: d62728FF
# SR/SL (2): ff9896FF
# WM: 1f77b4FF
col.pal1 = c("#2ca02cFF", "#98df8aFF","#ff7f0eFF", 
             "#C5B0D5FF","#d62728FF", "#ff9896FF")
col.pal2 = c("#1f77b4FF","#2ca02cFF","#98df8aFF", 
             "#ff7f0eFF","#C5B0D5FF","#d62728FF")

png(here("plots", "spatialHPC_SRT","precast_cluster_bias_sub6.png"), 
    width=5, height=5, units="in", res=300)
c1 <- lapply(seq_along(l2_sub6), function(x) {
  plotSpots(l2_sub6[[x]], annotate="precast_k6_ordered", point_size=.3)+
    labs(color="clus")+
    scale_color_manual(values=col.pal1)+
    theme(plot.title=element_text(size=8))
})
drawFigs(c1, layout.dim = c(2,3), common.legend = TRUE, 
        legend.position = "right", align = "h")
dev.off()

png(here("plots", "spatialHPC_SRT","precast_cluster_no_bias_sub6.png"), 
    width=5, height=5, units="in", res=300)
c2 <- lapply(seq_along(l2_sub6), function(x) {
  plotSpots(l2_sub6[[x]], annotate="precast_k6_nobias_ordered", point_size=.3)+
    labs(color="clus")+
    scale_color_manual(values=col.pal2)+
    theme(plot.title=element_text(size=8))
})
drawFigs(c2, layout.dim = c(2, 3), common.legend = TRUE, 
        legend.position = "right", align = "h")
dev.off()
