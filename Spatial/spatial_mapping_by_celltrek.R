#remotes::install_version(package = 'Matrix', version = package_version('1.6-3'))
#remotes::install_version(package = 'SeuratObject', version = package_version('4.1.3'))
#remotes::install_version(package = 'Seurat', version = package_version('4.3.0'))

options(stringsAsFactors = F)
library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")
library("hdf5r")

# Run drug group
stdata<- Load10X_Spatial(
  data.dir = "/path/to/spatial/data/h5",  
  filename = "filtered_feature_bc_matrix.h5"
)
scdata <- readRDS("/path/to/single/cell/rds")

# Run CellTrek
stdata <- RenameCells(stdata, new.names=make.names(Cells(stdata)))
scdata <- RenameCells(scdata, new.names=make.names(Cells(scdata)))

stdata <- SCTransform(stdata, assay = "Spatial", verbose = FALSE)
stdata <- RunPCA(stdata, assay = "SCT", verbose = FALSE)
stdata <- FindNeighbors(stdata, reduction = "pca", dims = 1:30)
stdata <- FindClusters(stdata, verbose = FALSE)
stdata <- RunUMAP(stdata, reduction = "pca", dims = 1:30)

scdata <- NormalizeData(scdata)
scdata <- FindVariableFeatures(scdata)
scdata <- ScaleData(scdata)
scdata <- RunPCA(scdata)
scdata <- RunUMAP(scdata, dims = 1:20) 
scdata <- FindNeighbors(scdata, dims = 1:20)
scdata <- FindClusters(scdata, resolution = 0.5)
scdata$cell_type <- Idents(scdata)

## Visualize the scRNA-seq data
DimPlot(scdata, label = T, label.size = 4.5)
traint <- CellTrek::traint(st_data=stdata, sc_data=scdata, sc_assay='RNA', cell_names='cell_type')
DimPlot(traint, group.by = "type") 
celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data=scdata, sc_assay = 'RNA', reduction='pca', intp=T, intp_pnt=1000, intp_lin=F, nPCs=15, ntree=300, dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T, num.threads = 16)$celltrek
celltrek$cell_type <- factor(celltrek$seurat_clusters, levels=sort(unique(celltrek$seurat_clusters)))

# Get mapped cellular coordinates
coord <- celltrek@meta.data[,c("coord_x","coord_y")]
# Get mapped cellular barcode
spatial_ids <- celltrek@meta.data$id_raw
valid_spatial_ids <- intersect(spatial_ids, colnames(scdata))
spatial_cells <- scdata[, valid_spatial_ids]

spatial_cells@assays$RNA@data <- spatial_cells@assays$RNA@counts
write.csv(spatial_cells@assays$RNA@data, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/drug/expression_matrix_3.csv")

# Save coordinates
write.csv(coord, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/drug/coord_3.csv")





# Run model group
stdata<- Load10X_Spatial(
  data.dir = "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/model/",  
  filename = "filtered_feature_bc_matrix.h5"
)
scdata <- readRDS("/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/scRNA/merged__scRNAseq2.rds")

stdata <- RenameCells(stdata, new.names=make.names(Cells(stdata)))
scdata <- RenameCells(scdata, new.names=make.names(Cells(scdata)))


stdata <- SCTransform(stdata, assay = "Spatial", verbose = FALSE)
stdata <- RunPCA(stdata, assay = "SCT", verbose = FALSE)
stdata <- FindNeighbors(stdata, reduction = "pca", dims = 1:30)
stdata <- FindClusters(stdata, verbose = FALSE)
stdata <- RunUMAP(stdata, reduction = "pca", dims = 1:30)

SpatialDimPlot(stdata)

scdata <- NormalizeData(scdata)
scdata <- FindVariableFeatures(scdata)
scdata <- ScaleData(scdata)
scdata <- RunPCA(scdata)
scdata <- RunUMAP(scdata, dims = 1:20) 
scdata <- FindNeighbors(scdata, dims = 1:20)
scdata <- FindClusters(scdata, resolution = 0.5)
scdata$cell_type <- Idents(scdata)

DimPlot(scdata, label = T, label.size = 4.5)

traint <- CellTrek::traint(st_data=stdata, sc_data=scdata, sc_assay='RNA', cell_names='cell_type')
DimPlot(traint, group.by = "type") 

celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data=scdata, sc_assay = 'RNA', reduction='pca', intp=T, intp_pnt=1000, intp_lin=F, nPCs=15, ntree=300, dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T, num.threads = 16)$celltrek
celltrek$cell_type <- factor(celltrek$seurat_clusters, levels=sort(unique(celltrek$seurat_clusters)))


coord <- celltrek@meta.data[,c("coord_x","coord_y")]

spatial_ids <- celltrek@meta.data$id_raw
valid_spatial_ids <- intersect(spatial_ids, colnames(scdata))
spatial_cells <- scdata[, valid_spatial_ids]

dim(spatial_cells)
spatial_cells@assays$RNA@data <- spatial_cells@assays$RNA@counts
write.csv(spatial_cells@assays$RNA@data, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/model/expression_matrix_2.csv")
write.csv(coord, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/model/coord_2.csv")





# run control group
stdata<- Load10X_Spatial(
  data.dir = "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/control/",  
  filename = "filtered_feature_bc_matrix.h5"
)
scdata <- readRDS("/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/scRNA/merged__scRNAseq2.rds")

stdata <- RenameCells(stdata, new.names=make.names(Cells(stdata)))
scdata <- RenameCells(scdata, new.names=make.names(Cells(scdata)))


stdata <- SCTransform(stdata, assay = "Spatial", verbose = FALSE)
stdata <- RunPCA(stdata, assay = "SCT", verbose = FALSE)
stdata <- FindNeighbors(stdata, reduction = "pca", dims = 1:30)
stdata <- FindClusters(stdata, verbose = FALSE)
stdata <- RunUMAP(stdata, reduction = "pca", dims = 1:30)

SpatialDimPlot(stdata)

scdata <- NormalizeData(scdata)
scdata <- FindVariableFeatures(scdata)
scdata <- ScaleData(scdata)
scdata <- RunPCA(scdata)
scdata <- RunUMAP(scdata, dims = 1:20) 
scdata <- FindNeighbors(scdata, dims = 1:20)
scdata <- FindClusters(scdata, resolution = 0.5)
scdata$cell_type <- Idents(scdata)

DimPlot(scdata, label = T, label.size = 4.5)

traint <- CellTrek::traint(st_data=stdata, sc_data=scdata, sc_assay='RNA', cell_names='cell_type')
DimPlot(traint, group.by = "type") 

celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data=scdata, sc_assay = 'RNA', reduction='pca', intp=T, intp_pnt=1000, intp_lin=F, nPCs=15, ntree=300, dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T, num.threads = 16)$celltrek

celltrek$cell_type <- factor(celltrek$seurat_clusters, levels=sort(unique(celltrek$seurat_clusters)))


coord <- celltrek@meta.data[,c("coord_x","coord_y")]

spatial_ids <- celltrek@meta.data$id_raw

valid_spatial_ids <- intersect(spatial_ids, colnames(scdata))
spatial_cells <- scdata[, valid_spatial_ids]

dim(spatial_cells)


spatial_cells@assays$RNA@data <- spatial_cells@assays$RNA@counts
write.csv(spatial_cells@assays$RNA@data, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/control/expression_matrix_1.csv")
write.csv(coord, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/control/coord_1.csv")
