##############################
#
# Load all used packages
#
# Written by Beáta Szeitz
# Project: SCLC proteomics
#
##############################

packages = c("ggplot2","ggbiplot","ggplotify", "ggpubr","lme4", "circlize","ComplexHeatmap", "patchwork","pvca","Biobase",
             "readxl","openxlsx", "edgeR", "clusterProfiler","ReactomePA", "org.Hs.eg.db","ggrepel", "RColorBrewer",
             "enrichplot", "gridExtra", "cowplot", "ConsensusClusterPlus","amap", "cluster", "S4Vectors", "mixOmics",
             "reshape2", "varhandle", "DTK", "car")

lapply(packages, function(x) {library(x, character.only = T)})
