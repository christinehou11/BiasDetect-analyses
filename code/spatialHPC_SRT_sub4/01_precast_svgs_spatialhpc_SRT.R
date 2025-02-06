suppressPackageStartupMessages({
    library(Seurat)
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
# choose the subset (4 samples) from raw SRT data
# ---------
fix_order = distinct(as.data.frame(colData(spe)), slide, array, brnum, sample_id, position, sex) %>% 
  arrange(slide, array)
sub4 = fix_order$sample_id[c(14,16,
                             20,21)]
spe_sub4 = spe[,spe$sample_id %in% sub4] # 31483, 18945

# ---------
# filtered the SVGs from selected samples
# ---------
load(here("raw-data","spatialHPC_SRT","nnSVG_outs_HE_only.rda"))

res_df_sub4 <- pivot_longer(
  rownames_to_column(as.data.frame(res_ranks), var<-"gene_id"), 
  colnames(res_ranks), names_to="sample_id", values_to="rank", 
  values_drop_na=T)

res_df2_sub4 <- filter(res_df_sub4,
    sample_id %in% c("V11L05-333_B1","V11L05-333_D1","V11L05-335_D1","V11L05-336_A1"),
        rank <= 2000) # top 2k sig features
nrow(res_df2_sub4) # 7559, 3

svgs_sub4 <- group_by(res_df2_sub4, gene_id) %>% 
  tally() %>% 
  filter(n>1) # >1 sample

nrow(svgs_sub4) # 2082

# ---------
# reformatted to SpatialExperiment object
# ---------
spe_sub4 <- spe_sub4[rowData(spe_sub4)$gene_id %in% svgs_sub4$gene_id,]
rownames(spe_sub4) <- rowData(spe_sub4)$gene_id
dim(spe_sub4) # 2082, 18945

# ---------
# reformat to seurat list
# ---------
l2_sub4 = unique(spe_sub4$sample_id)
names(l2_sub4) = l2_sub4
l2_sub4 = lapply(l2_sub4, function(x) spe_sub4[,colData(spe_sub4)$sample_id==x])

srt.sets.sub4 = lapply(l2_sub4, function(x) {
    colnames(counts(x)) <- rownames(colData(x))
    colData(x)$col <- spatialCoords(x)[,"pxl_col_in_fullres"]
    colData(x)$row <- spatialCoords(x)[,"pxl_row_in_fullres"]
    count <- counts(x)
    a1 <- CreateAssayObject(count, assay = "RNA", 
                            min.features = 0, min.cells = 0)
    CreateSeuratObject(a1, meta.data = as.data.frame(colData(x)))
})
srt.sets.sub4
# $`V11L05-333_B1`
# An object of class Seurat 
# 2082 features across 4985 samples within 1 assay 
# Active assay: RNA (2082 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`V11L05-333_D1`
# An object of class Seurat 
# 2082 features across 4938 samples within 1 assay 
# Active assay: RNA (2082 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`V11L05-335_D1`
# An object of class Seurat 
# 2082 features across 4483 samples within 1 assay 
# Active assay: RNA (2082 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`V11L05-336_A1`
# An object of class Seurat 
# 2082 features across 4539 samples within 1 assay 
# Active assay: RNA (2082 features, 0 variable features)
# 2 layers present: counts, data

# ---------
# run precast
# ---------
preobj_sub4 <- CreatePRECASTObject(seuList = srt.sets.sub4,
                            customGenelist=rownames(spe_sub4),
                            premin.spots=0, premin.features=0, 
                            postmin.spots=0, postmin.features=0)
PRECASTObj_sub4 <- AddAdjList(preobj_sub4, platform = "Visium")
PRECASTObj_sub4 <- AddParSetting(PRECASTObj_sub4, maxIter = 20, 
                            verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj_sub4 <- PRECAST(PRECASTObj_sub4, K=7)

#consolidate/ reformat results
PRECASTObj_sub4 <- SelectModel(PRECASTObj_sub4, criteria="MBIC")
PRECASTObj_sub4
# An object of class PRECASTObj 
# with 4 datasets and  18945 spots in total, with spots for each dataset:  4985 4938 4483 4539 
# 2081 common variable genes selected

seuInt_sub4 <- IntegrateSpaData(PRECASTObj_sub4, species = "Human")
seuInt_sub4
# An object of class Seurat 
# 2081 features across 18945 samples within 1 assay 
# Active assay: PRE_CAST (2081 features, 0 variable features)
# 2 layers present: counts, data
# 2 dimensional reductions calculated: PRECAST, position

# -----------
# save object
# -----------
save(spe_sub4,
    file = here("processed-data","spatialHPC_SRT","spe_sub4.rda"))

save(srt.sets.sub4,
    file = here("processed-data","spatialHPC_SRT","srt_sets_sub4.rda"))

save(PRECASTObj_sub4,
     file = here("processed-data","spatialHPC_SRT","PRECASTObj_sub4.rda"))

save(seuInt_sub4,
     file = here("processed-data","spatialHPC_SRT","seuInt_sub4.rda"))

write.csv(seuInt_sub4@meta.data, 
        here("processed-data","spatialHPC_SRT",
        "seuInt_sub4-hpc_k-7_svgs_metadata.csv"), 
        row.names=T)
