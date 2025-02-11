suppressPackageStartupMessages({
    library(Seurat)
    library(SummarizedExperiment)
    library(SpatialExperiment)
    library(PRECAST)
    library(dplyr)
    library(here)
    library(ExperimentHub)
    library(tidyr)
    library(tibble)
})
set.seed(123)

# ---------
# load data
# ---------
ehub <- ExperimentHub()
spe <- ehub[["EH9605"]]
spe
# class: SpatialExperiment 
# dim: 31483 150917 
# metadata(1): Obtained_from
# assays(2): counts logcounts
# rownames(31483): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type gene_search
# colnames(150917): AAACAACGAATAGTTC-1_V10B01-086_D1
# AAACAAGTATCTCCCA-1_V10B01-086_D1 ... TTGTTTCCATACAACT-1_Br2720_B1
# TTGTTTGTATTACACG-1_Br2720_B1
# colData names(150): sample_id in_tissue ... nmf99 nmf100
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

# ---------
# choose the subset (2 samples) from raw SRT data
# ---------
fix_order = distinct(as.data.frame(colData(spe)), slide, array, brnum, sample_id, position, sex) %>% 
  arrange(slide, array)
sub2 = fix_order$sample_id[c(14,16)]
spe_sub2 = spe[,spe$sample_id %in% sub2] # 31483, 9923

# ---------
# make the features ratio as 2:1 (2/3 from 14 and 1/3 from 16)
# ---------
cols_14 <- which(spe$sample_id == fix_order$sample_id[14]) # 4985
cols_16 <- which(spe$sample_id == fix_order$sample_id[16]) # 4938

num_14 <- round((2/3) * length(cols_14)) # 3323
num_16 <- round((1/3) * length(cols_16)) # 1646

selected_cols_14 <- sample(cols_14, num_14)
selected_cols_16 <- sample(cols_16, num_16)
selected_cols <- c(selected_cols_14, selected_cols_16) # 4969

spe_sub2_2<- spe[, selected_cols] # 31483  4969
spe_sub2_2
# class: SpatialExperiment 
# dim: 31483 4969 
# metadata(1): Obtained_from
# assays(2): counts logcounts
# rownames(31483): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type gene_search
# colnames(4969): GTTAACTATGTTGTCA-1_V11L05-333_B1 ACGTAGATTGCTGATG-1_V11L05-333_B1 ...
# TGGCTACACTCTACCT-1_V11L05-333_D1 GGTCGGTAATTAGACA-1_V11L05-333_D1
# colData names(150): sample_id in_tissue ... nmf99 nmf100
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

# ---------
# filtered the SVGs from selected samples
# ---------
load(here("raw-data","spatialHPC_SRT","nnSVG_outs_HE_only.rda"))
res_ranks

res_df_sub2_2 <- pivot_longer(
  rownames_to_column(as.data.frame(res_ranks), var<-"gene_id"), 
  colnames(res_ranks), names_to="sample_id", values_to="rank", 
  values_drop_na=T)

res_df2_sub2_2 <- filter(res_df_sub2_2,
    sample_id %in% c("V11L05-333_B1","V11L05-333_D1"),
        rank <= 2000) # top 2k sig features
nrow(res_df2_sub2_2) # 4000, 3

svgs_sub2_2 <- group_by(res_df2_sub2_2, gene_id) %>% 
  tally() %>% 
  filter(n>1) # >1 sample
nrow(svgs_sub2_2) # 1498

# ---------
# reformatted to SpatialExperiment object
# ---------
spe_sub2_2 <- spe_sub2_2[rowData(spe_sub2_2)$gene_id %in% svgs_sub2_2$gene_id,]
rownames(spe_sub2_2) <- rowData(spe_sub2_2)$gene_id
dim(spe_sub2_2) # 1498, 4969

# ---------
# reformat to seurat list
# ---------
l2_sub2_2 = unique(spe_sub2_2$sample_id)
names(l2_sub2_2) = l2_sub2_2
l2_sub2_2 = lapply(l2_sub2_2, function(x) spe_sub2_2[,colData(spe_sub2_2)$sample_id==x])

srt.sets.sub2_2 = lapply(l2_sub2_2, function(x) {
    colnames(counts(x)) <- rownames(colData(x))
    colData(x)$col <- spatialCoords(x)[,"pxl_col_in_fullres"]
    colData(x)$row <- spatialCoords(x)[,"pxl_row_in_fullres"]
    count <- counts(x)
    a1 <- CreateAssayObject(count, assay = "RNA", 
                            min.features = 0, min.cells = 0)
    CreateSeuratObject(a1, meta.data = as.data.frame(colData(x)))
})
srt.sets.sub2_2
# $`V11L05-333_B1`
# An object of class Seurat 
# 1498 features across 3323 samples within 1 assay 
# Active assay: RNA (1498 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`V11L05-333_D1`
# An object of class Seurat 
# 1498 features across 1646 samples within 1 assay 
# Active assay: RNA (1498 features, 0 variable features)
# 2 layers present: counts, data

# ---------
# run precast
# ---------
preobj_sub2_2 <- CreatePRECASTObject(seuList = srt.sets.sub2_2,
                            customGenelist=rownames(spe_sub2_2),
                            premin.spots=0, premin.features=0, 
                            postmin.spots=0, postmin.features=0)
PRECASTObj_sub2_2 <- AddAdjList(preobj_sub2_2, platform = "Visium")
PRECASTObj_sub2_2 <- AddParSetting(PRECASTObj_sub2_2, maxIter = 20, 
                            verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj_sub2_2 <- PRECAST(PRECASTObj_sub2_2, K=7)

#consolidate/ reformat results
PRECASTObj_sub2_2 <- SelectModel(PRECASTObj_sub2_2, criteria="MBIC")
PRECASTObj_sub2_2
# An object of class PRECASTObj 
# with 2 datasets and  4969 spots in total, with spots for each dataset:  3323 1646  
# 1498 common variable genes selected

seuInt_sub2_2 <- IntegrateSpaData(PRECASTObj_sub2_2, species = "Human")
seuInt_sub2_2
# An object of class Seurat 
# 1498 features across 4969 samples within 1 assay 
# Active assay: PRE_CAST (1498 features, 0 variable features)
# 2 layers present: counts, data
# 2 dimensional reductions calculated: PRECAST, position

# -----------
# save object
# -----------
save(spe_sub2_2,
    file = here("processed-data","spatialHPC_SRT","spe_sub2_2.rda"))

save(srt.sets.sub2_2,
    file = here("processed-data","spatialHPC_SRT","srt_sets_sub2_2.rda"))

save(PRECASTObj_sub2_2,
     file = here("processed-data","spatialHPC_SRT","PRECASTObj_sub2_2.rda"))

save(seuInt_sub2_2,
     file = here("processed-data","spatialHPC_SRT","seuInt_sub2_2.rda"))

write.csv(seuInt_sub2_2@meta.data, 
        here("processed-data","spatialHPC_SRT",
        "seuInt_sub2_2-hpc_k-7_svgs_metadata.csv"), 
        row.names=T)
