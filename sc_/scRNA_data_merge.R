# This script was used to merge scRNA-seq samples
library(Seurat)
library(tools)
library(stringr)

setwd("/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/scRNA/")

untar("GSE202297_RAW.tar", exdir = "GSE202297_RAW2")
untar("GSE98969_RAW.tar", exdir = "GSE98969_RAW2")

dir1 <- "GSE202297_RAW2"
dir2 <- "GSE98969_RAW2"

files <- c(
  list.files(dir1, pattern = "*.txt.gz$", full.names = TRUE),
  list.files(dir2, pattern = "*.txt.gz$", full.names = TRUE)
)

sample_ids <- basename(files) |> 
  file_path_sans_ext() |> 
  str_remove("\\.txt$")

read_as_seurat <- function(file, sample_id) 
{
  expr <- read.delim(gzfile(file), row.names = 1, check.names = FALSE)
  expr <- as.matrix(expr)
  storage.mode(expr) <- "integer"
  
  colnames(expr) <- paste(sample_id, colnames(expr), sep = "_")
  
  seu <- CreateSeuratObject(counts = expr, project = sample_id)
  seu$sample <- sample_id
  
  if (grepl("GSE202297", file)) {
    seu$group <- "GSE202297"
  } else if (grepl("GSE98969", file)) {
    seu$group <- "GSE98969"
  }
  
  return(seu)
}

seurat_list <- mapply(read_as_seurat, files, sample_ids, SIMPLIFY = FALSE)

seurat_merged <- Reduce(function(x, y) merge(x, y), seurat_list)

saveRDS(seurat_merged, file = "merged__scRNAseq2.rds")


