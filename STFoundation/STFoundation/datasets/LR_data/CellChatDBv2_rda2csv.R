## 2025.10.15 convert CellChatDBv2 rda to csv files
## CellChatDBv2 download from https://github.com/jinworks/CellChat/tree/main/data

#!/usr/bin/env Rscript

# --- Parse command line arguments ---
args <- commandArgs(trailingOnly = TRUE)
species <- ifelse(length(args) >= 1, args[[1]], "human")    # Default to "human"
rda_file <- ifelse(length(args) >= 2, args[[2]], paste0("CellChatDB.", species, ".rda"))

rda_path <- file.path("0_CellChatDBv2", rda_file)

cat("Species:", species, "\n")
cat("RDA file:", rda_file, "\n")
cat("RDA path:", rda_path, "\n")
cat("Current working directory:", getwd(), "\n")

# --- Load CellChatDB RDA file from ./0_CellChatDBv2/ ---
if (!file.exists(rda_path)) {
    stop(paste("File not found:", rda_path))
}
load(rda_path)

# --- Detect loaded object automatically ---
obj_name <- ls()[grepl("^CellChatDB\\.", ls())]
if (length(obj_name) == 0) stop("No CellChatDB.* object found in the RDA file!")
cat("Loaded object:", obj_name, "\n")
CellChatDB <- get(obj_name)

# --- Show structure ---
str(CellChatDB, max.level = 1)

# --- Set output prefix and save as CSV in ./LR_data/ (create if missing) ---
out_dir <- "STFoundation/datasets/LR_data"
if (!dir.exists(out_dir)) dir.create(out_dir)

database <- paste0(species, "-CellChatDBv2_")
write.csv(CellChatDB$interaction, file = file.path(out_dir, paste0(database, "interaction.csv")), row.names = FALSE)
write.csv(CellChatDB$complex,     file = file.path(out_dir, paste0(database, "complex.csv")),     row.names = FALSE)
write.csv(CellChatDB$cofactor,    file = file.path(out_dir, paste0(database, "cofactor.csv")),    row.names = FALSE)
write.csv(CellChatDB$geneInfo,    file = file.path(out_dir, paste0(database, "geneInfo.csv")),    row.names = FALSE)

cat("Exported CellChatDB tables to", out_dir, "\n")


## usage:
# 在含有 CellChatDB.human.rda 的目录下运行
# cd /mnt/lingyu/nfs_share2/Python/STPatch/KANsignal_GPU2/STFoundation
# Rscript STFoundation/datasets/LR_data/CellChatDBv2_rda2csv.R mouse 