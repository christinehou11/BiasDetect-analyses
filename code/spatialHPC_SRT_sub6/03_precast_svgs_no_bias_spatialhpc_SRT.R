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
load(here("processed-data","spatialHPC_SRT","spe.rda"))
load(here("processed-data","spatialHPC_SRT","srt_sets.rda"))
load(here("processed-data","spatialHPC_SRT","hpc_biased_features.rda"))

# ---------
# remove biased genes
# ---------
svgs_filt = setdiff(rownames(spe), bias)
length(svgs_filt) #should be 2067

# ---------
# run precast
# ---------
preobj_filt <- CreatePRECASTObject(seuList = srt.sets,
                              customGenelist=svgs_filt,
                              premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj_filt <- AddAdjList(preobj_filt, platform = "Visium")
PRECASTObj_filt <- AddParSetting(PRECASTObj_filt, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj_filt <- PRECAST(PRECASTObj_filt, K=7)

#consolidate/ reformat results
PRECASTObj_filt <- SelectModel(PRECASTObj_filt, criteria="MBIC")
# An object of class PRECASTObj 
# with 4 datasets and  18945 spots in total, with spots for each dataset:  4985 4938 4483 4539 
# 2067 common variable genes selected

seuInt_filt <- IntegrateSpaData(PRECASTObj_filt, species = "Human")
# An object of class Seurat 
# 2067 features across 18945 samples within 1 assay 
# Active assay: PRE_CAST (2067 features, 0 variable features)
# 2 layers present: counts, data
# 2 dimensional reductions calculated: PRECAST, position

# ---------
# save object
# ---------
write.csv(seuInt_filt@meta.data, 
    here("processed-data","spatialHPC_SRT",
        "seuInt-hpc_k-7_svgs-no-bias_metadata.csv"), row.names=T)

