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
load(here("processed-data","spatialHPC_SRT","spe_sub6.rda"))
load(here("processed-data","spatialHPC_SRT","srt_sets_sub6.rda"))
load(here("processed-data","spatialHPC_SRT","hpc_biased_features_sub6.rda"))

# ---------
# remove biased genes
# ---------
svgs_filt_sub6 = setdiff(rownames(spe_sub6), bias_sub6)
length(svgs_filt_sub6) #should be 2316

# ---------
# run precast
# ---------
preobj_filt_sub6 <- CreatePRECASTObject(seuList = srt.sets.sub6,
                              customGenelist=svgs_filt_sub6,
                              premin.spots=0, premin.features=0, 
                              postmin.spots=0, postmin.features=0)
PRECASTObj_filt_sub6 <- AddAdjList(preobj_filt_sub6, platform = "Visium")
PRECASTObj_filt_sub6 <- AddParSetting(PRECASTObj_filt_sub6, 
                                      maxIter = 20, verbose = TRUE, 
                                      Sigma_equal=FALSE, coreNum=12)
PRECASTObj_filt_sub6 <- PRECAST(PRECASTObj_filt_sub6, K=6)

#consolidate/ reformat results
PRECASTObj_filt_sub6 <- SelectModel(PRECASTObj_filt_sub6, criteria="MBIC")
# An object of class PRECASTObj 
# with 6 datasets and  28306 spots in total, with spots for each dataset:
# 4969 4985 4887 4483 4539 4443
# 2316 common variable genes selected

seuInt_filt_sub6 <- IntegrateSpaData(PRECASTObj_filt_sub6, species = "Human")
# An object of class Seurat 
# 2316 features across 28306 samples within 1 assay 
# Active assay: PRE_CAST (2316 features, 0 variable features)
# 2 layers present: counts, data
# 2 dimensional reductions calculated: PRECAST, position

# ---------
# save object
# ---------
save(PRECASTObj_filt_sub6,
     file = here("processed-data","spatialHPC_SRT","PRECASTObj_filt_sub6.rda"))

save(seuInt_filt_sub6,
     file = here("processed-data","spatialHPC_SRT","seuInt_filt_sub6.rda"))

write.csv(seuInt_filt_sub6@meta.data, 
    here("processed-data","spatialHPC_SRT",
        "seuInt_sub6-hpc_k-7_svgs-no-bias_metadata.csv"), row.names=T)
