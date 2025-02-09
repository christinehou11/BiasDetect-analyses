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
# choose the subset (6 samples) from raw SRT data
# ---------
fix_order = distinct(as.data.frame(colData(spe)), slide, array, brnum, sample_id, position, sex) %>% 
  arrange(slide, array)
sub6 = fix_order$sample_id[c(13,14,19,20,21,23)]
spe_sub6 = spe[,spe$sample_id %in% sub6] # 31483, 28306

# ---------
# filtered the SVGs from selected samples
# ---------
load(here("raw-data","spatialHPC_SRT","nnSVG_outs_HE_only.rda"))

res_df <- pivot_longer(
  rownames_to_column(as.data.frame(res_ranks), var<-"gene_id"), 
  colnames(res_ranks), names_to="sample_id", values_to="rank", 
  values_drop_na=T)

res_df_sub6 <- filter(res_df,
    sample_id %in% c("V11L05-333_A1","V11L05-333_B1",
                     "V11L05-335_C1","V11L05-335_D1",
                     "V11L05-336_A1","V11L05-336_C1"),
        rank <= 2000) # top 2k sig features
nrow(res_df_sub6) # 10914

svgs_sub6 <- group_by(res_df_sub6, gene_id) %>% 
  tally() %>% 
  filter(n>1) # >1 sample

nrow(svgs_sub6) # 2333

# ---------
# reformatted to SpatialExperiment object
# ---------
spe_sub6 <- spe_sub6[rowData(spe_sub6)$gene_id %in% svgs_sub6$gene_id,]
rownames(spe_sub6) <- rowData(spe_sub6)$gene_id
dim(spe_sub6) # 2333, 28306

# ---------
# reformat to seurat list
# ---------
l2_sub6 = unique(spe_sub6$sample_id)
names(l2_sub6) = l2_sub6
l2_sub6 = lapply(l2_sub6, function(x) spe_sub6[,colData(spe_sub6)$sample_id==x])

srt.sets.sub6 = lapply(l2_sub6, function(x) {
    colnames(counts(x)) <- rownames(colData(x))
    colData(x)$col <- spatialCoords(x)[,"pxl_col_in_fullres"]
    colData(x)$row <- spatialCoords(x)[,"pxl_row_in_fullres"]
    count <- counts(x)
    a1 <- CreateAssayObject(count, assay = "RNA", 
                            min.features = 0, min.cells = 0)
    CreateSeuratObject(a1, meta.data = as.data.frame(colData(x)))
})
srt.sets.sub6
# $`V11L05-333_A1`
# An object of class Seurat 
# 2333 features across 4969 samples within 1 assay 
# Active assay: RNA (2333 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`V11L05-333_B1`
# An object of class Seurat 
# 2333 features across 4985 samples within 1 assay 
# Active assay: RNA (2333 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`V11L05-335_C1`
# An object of class Seurat 
# 2333 features across 4887 samples within 1 assay 
# Active assay: RNA (2333 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`V11L05-335_D1`
# An object of class Seurat 
# 2333 features across 4483 samples within 1 assay 
# Active assay: RNA (2333 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`V11L05-336_A1`
# An object of class Seurat 
# 2333 features across 4539 samples within 1 assay 
# Active assay: RNA (2333 features, 0 variable features)
# 2 layers present: counts, data
# 
# $`V11L05-336_C1`
# An object of class Seurat 
# 2333 features across 4443 samples within 1 assay 
# Active assay: RNA (2333 features, 0 variable features)
# 2 layers present: counts, data

# ---------
# run precast
# ---------
preobj_sub6 <- CreatePRECASTObject(seuList = srt.sets.sub6,
                            customGenelist=rownames(spe_sub6),
                            premin.spots=0, premin.features=0, 
                            postmin.spots=0, postmin.features=0)
PRECASTObj_sub6 <- AddAdjList(preobj_sub6, platform = "Visium")
PRECASTObj_sub6 <- AddParSetting(PRECASTObj_sub6, maxIter = 20, 
                            verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj_sub6 <- PRECAST(PRECASTObj_sub6, K=6)

#consolidate/ reformat results
PRECASTObj_sub6 <- SelectModel(PRECASTObj_sub6, criteria="MBIC")
PRECASTObj_sub6
# An object of class PRECASTObj 
# with 6 datasets and 28306 spots in total, with spots for each dataset:  
# 4969 4985 4887 4483 4539 4443
# 2332 common variable genes selected

seuInt_sub6 <- IntegrateSpaData(PRECASTObj_sub6, species = "Human")
seuInt_sub6
# An object of class Seurat 
# 2332 features across 28306 samples within 1 assay 
# Active assay: PRE_CAST (2332 features, 0 variable features)
# 2 layers present: counts, data
# 2 dimensional reductions calculated: PRECAST, position

# -----------
# save object
# -----------
save(spe_sub6,
    file = here("processed-data","spatialHPC_SRT","spe_sub6.rda"))

save(srt.sets.sub6,
    file = here("processed-data","spatialHPC_SRT","srt_sets_sub6.rda"))

save(PRECASTObj_sub6,
     file = here("processed-data","spatialHPC_SRT","PRECASTObj_sub6.rda"))

save(seuInt_sub6,
     file = here("processed-data","spatialHPC_SRT","seuInt_sub6.rda"))

write.csv(seuInt_sub6@meta.data, 
        here("processed-data","spatialHPC_SRT",
        "seuInt_sub6-hpc_k-7_svgs_metadata.csv"), 
        row.names=T)
