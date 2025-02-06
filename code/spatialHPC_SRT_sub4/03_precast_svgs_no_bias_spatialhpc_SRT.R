suppressPackageStartupMessages({
    library(Seurat)
    library(PRECAST)
    library(dplyr)
    library(here)
})
set.seed(123)

# ---------
# load data
# ---------
load(here("processed-data","spatialHPC_SRT","spe_sub4.rda"))
load(here("processed-data","spatialHPC_SRT","srt_sets_sub4.rda"))
load(here("processed-data","spatialHPC_SRT","hpc_biased_features_sub4.rda"))

# ---------
# remove biased genes
# ---------
svgs_filt_sub4 = setdiff(rownames(spe_sub4), bias)
length(svgs_filt_sub4) #should be 2067

# ---------
# run precast
# ---------
preobj_filt_sub4 <- CreatePRECASTObject(seuList = srt.sets.sub4,
                              customGenelist=svgs_filt_sub4,
                              premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj_filt_sub4 <- AddAdjList(preobj_filt_sub4, platform = "Visium")
PRECASTObj_filt_sub4 <- AddParSetting(PRECASTObj_filt_sub4, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj_filt_sub4 <- PRECAST(PRECASTObj_filt_sub4, K=7)

#consolidate/ reformat results
PRECASTObj_filt_sub4 <- SelectModel(PRECASTObj_filt_sub4, criteria="MBIC")
# An object of class PRECASTObj 
# with 4 datasets and  18945 spots in total, with spots for each dataset:  4985 4938 4483 4539 
# 2067 common variable genes selected

seuInt_filt_sub4 <- IntegrateSpaData(PRECASTObj_filt_sub4, species = "Human")
# An object of class Seurat 
# 2067 features across 18945 samples within 1 assay 
# Active assay: PRE_CAST (2067 features, 0 variable features)
# 2 layers present: counts, data
# 2 dimensional reductions calculated: PRECAST, position

# ---------
# save object
# ---------
save(PRECASTObj_filt_sub4,
     file = here("processed-data","spatialHPC_SRT","PRECASTObj_filt_sub4.rda"))

save(seuInt_filt_sub4,
     file = here("processed-data","spatialHPC_SRT","seuInt_filt_sub4.rda"))

write.csv(seuInt_filt_sub4@meta.data, 
    here("processed-data","spatialHPC_SRT",
        "seuInt_sub4-hpc_k-7_svgs-no-bias_metadata.csv"), row.names=T)

