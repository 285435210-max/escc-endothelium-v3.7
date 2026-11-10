suppressPackageStartupMessages({library(Seurat)})
OUT <- "D:/escc_sc/out"
obj <- readRDS(file.path(OUT, "endo_annotated.rds"))
obj <- subset(obj, subset = !grepl("^Contam", as.character(celltype)))
lig <- c("CXCL9","CXCL10","CXCL11","CXCL16","CCL5","CX3CL1","DLL1","DLL4","JAG1",
         "VEGFA","VEGFB","ICAM1","VCAM1","SELE","SELP","HLA-A","HLA-B","HLA-C",
         "HLA-DRA","CD274","PDCD1LG2","LGALS9","TNFRSF14","TNFSF9","CD70","CXCL12")
lig <- intersect(lig, rownames(obj))
av <- as.matrix(AverageExpression(obj, features = lig, group.by = "celltype", slot = "data")$RNA)
write.csv(round(av, 3), file.path(OUT, "ec_ligand_mean_by_celltype.csv"))
cat("ligands:", length(lig), " celltypes:", ncol(av), "\n")
print(round(av, 2))
