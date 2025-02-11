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
load(here("processed-data","spatialHPC_SRT","spe_sub2.rda"))
load(here("processed-data","spatialHPC_SRT","srt_sets_sub2.rda"))
load(here("processed-data","spatialHPC_SRT","hpc_biased_features_sub2.rda"))

# ---------
# remove biased genes
# ---------
svgs_filt_sub2 = setdiff(rownames(spe_sub2), bias)
length(svgs_filt_sub2) #should be 1485

# ---------
# run precast
# ---------
preobj_filt_sub2 <- CreatePRECASTObject(seuList = srt.sets.sub2,
                              customGenelist=svgs_filt_sub2,
                              premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj_filt_sub2 <- AddAdjList(preobj_filt_sub2, platform = "Visium")
PRECASTObj_filt_sub2 <- AddParSetting(PRECASTObj_filt_sub2, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj_filt_sub2 <- PRECAST(PRECASTObj_filt_sub2, K=7)

#consolidate/ reformat results
PRECASTObj_filt_sub2 <- SelectModel(PRECASTObj_filt_sub2, criteria="MBIC")
# An object of class PRECASTObj 
# with 2 datasets and  9923 spots in total, with spots for each dataset:  4985 4938 4483 4539 
# 1485 common variable genes selected

seuInt_filt_sub2 <- IntegrateSpaData(PRECASTObj_filt_sub2, species = "Human")
# An object of class Seurat 
# 1485 features across 9923 samples within 1 assay 
# Active assay: PRE_CAST (1485 features, 0 variable features)
# 2 layers present: counts, data
# 2 dimensional reductions calculated: PRECAST, position

# ---------
# save object
# ---------
save(PRECASTObj_filt_sub2,
     file = here("processed-data","spatialHPC_SRT","PRECASTObj_filt_sub2.rda"))

save(seuInt_filt_sub2,
     file = here("processed-data","spatialHPC_SRT","seuInt_filt_sub2.rda"))

write.csv(seuInt_filt_sub2@meta.data, 
    here("processed-data","spatialHPC_SRT",
        "seuInt_sub2-hpc_k-7_svgs-no-bias_metadata.csv"), row.names=T)

