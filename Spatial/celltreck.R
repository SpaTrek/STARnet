# 这几个包的版本要换成这些版本，其他包的版本都不用动
#remotes::install_version(package = 'Matrix', version = package_version('1.6-3'))
#remotes::install_version(package = 'SeuratObject', version = package_version('4.1.3'))
#remotes::install_version(package = 'Seurat', version = package_version('4.3.0'))

# 设置多个包路径，路径的顺序决定了优先加载哪个路径
.libPaths(c("/slurm/home/yrd/fanlab/gaomeng/R/x86_64-pc-linux-gnu-library/4.3/", 
            "/slurm/home/yrd/fanlab/gaomeng/anaconda3/lib/R/library"))

options(stringsAsFactors = F)
library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")
library("hdf5r")
stdata<- Load10X_Spatial(
  data.dir = "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/drug/",  
  filename = "filtered_feature_bc_matrix.h5"
)
scdata <- readRDS("/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/scRNA/merged__scRNAseq2.rds")

#随后安装官网运行指令进行CellTrek分析
stdata <- RenameCells(stdata, new.names=make.names(Cells(stdata)))
scdata <- RenameCells(scdata, new.names=make.names(Cells(scdata)))


stdata <- SCTransform(stdata, assay = "Spatial", verbose = FALSE)
stdata <- RunPCA(stdata, assay = "SCT", verbose = FALSE)
stdata <- FindNeighbors(stdata, reduction = "pca", dims = 1:30)
stdata <- FindClusters(stdata, verbose = FALSE)
stdata <- RunUMAP(stdata, reduction = "pca", dims = 1:30)

# 现在重新画
SpatialDimPlot(stdata)

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
a<-c
traint <- CellTrek::traint(st_data=stdata, sc_data=scdata, sc_assay='RNA', cell_names='cell_type')
DimPlot(traint, group.by = "type") 
cat("🟢 正在开始 celltrek 映射...\n")
start_time <- Sys.time()
celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data=scdata, sc_assay = 'RNA', reduction='pca', intp=T, intp_pnt=1000, intp_lin=F, nPCs=15, ntree=300, dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T, num.threads = 16)$celltrek
cat("✅ celltrek 映射完成\n")
print(Sys.time() - start_time)
celltrek$cell_type <- factor(celltrek$seurat_clusters, levels=sort(unique(celltrek$seurat_clusters)))

#接卷积结果的可视化，可以不运行
#CellTrek::celltrek_vis(celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_type:id_new),
#                       celltrek@images$D1.AD12.L1@image, celltrek@images$D1.AD12.L1@scale.factors$lowres)


#获取解卷积的细胞以及其坐标
coord <- celltrek@meta.data[,c("coord_x","coord_y")]
# 获取解卷积出的细胞 ID（字符向量）
spatial_ids <- celltrek@meta.data$id_raw
# 保证这些细胞 ID 在 scdata 里确实存在
valid_spatial_ids <- intersect(spatial_ids, colnames(scdata))
# 仅提取 scdata 中存在的细胞
spatial_cells <- scdata[, valid_spatial_ids]

dim(spatial_cells)


spatial_cells@assays$RNA@data <- spatial_cells@assays$RNA@counts
write.csv(spatial_cells@assays$RNA@data, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/drug/spatial__adata2.csv")

#Seurat对象转为h5ad
#SaveH5Seurat(spatial_cells, filename = "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial//spatial_cells1.h5Seurat")
#Convert("/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/control/spatial_cells1.h5Seurat", dest = "h5ad")

#坐标存储
write.csv(coord, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/drug/coord2_.csv")


#model
stdata<- Load10X_Spatial(
  data.dir = "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/model/",  
  filename = "filtered_feature_bc_matrix.h5"
)
scdata <- readRDS("/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/scRNA/merged__scRNAseq2.rds")

#随后安装官网运行指令进行CellTrek分析
stdata <- RenameCells(stdata, new.names=make.names(Cells(stdata)))
scdata <- RenameCells(scdata, new.names=make.names(Cells(scdata)))


stdata <- SCTransform(stdata, assay = "Spatial", verbose = FALSE)
stdata <- RunPCA(stdata, assay = "SCT", verbose = FALSE)
stdata <- FindNeighbors(stdata, reduction = "pca", dims = 1:30)
stdata <- FindClusters(stdata, verbose = FALSE)
stdata <- RunUMAP(stdata, reduction = "pca", dims = 1:30)

# 现在重新画
SpatialDimPlot(stdata)

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
a<-c
traint <- CellTrek::traint(st_data=stdata, sc_data=scdata, sc_assay='RNA', cell_names='cell_type')
DimPlot(traint, group.by = "type") 
cat("🟢 正在开始 celltrek 映射...\n")
start_time <- Sys.time()
celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data=scdata, sc_assay = 'RNA', reduction='pca', intp=T, intp_pnt=1000, intp_lin=F, nPCs=15, ntree=300, dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T, num.threads = 16)$celltrek
cat("✅ celltrek 映射完成\n")
print(Sys.time() - start_time)
celltrek$cell_type <- factor(celltrek$seurat_clusters, levels=sort(unique(celltrek$seurat_clusters)))

#接卷积结果的可视化，可以不运行
#CellTrek::celltrek_vis(celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_type:id_new),
#                       celltrek@images$D1.AD12.L1@image, celltrek@images$D1.AD12.L1@scale.factors$lowres)


#获取解卷积的细胞以及其坐标
coord <- celltrek@meta.data[,c("coord_x","coord_y")]
# 获取解卷积出的细胞 ID（字符向量）
spatial_ids <- celltrek@meta.data$id_raw
# 保证这些细胞 ID 在 scdata 里确实存在
valid_spatial_ids <- intersect(spatial_ids, colnames(scdata))
# 仅提取 scdata 中存在的细胞
spatial_cells <- scdata[, valid_spatial_ids]

dim(spatial_cells)


spatial_cells@assays$RNA@data <- spatial_cells@assays$RNA@counts
write.csv(spatial_cells@assays$RNA@data, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/model/spatial__adata2.csv")

#Seurat对象转为h5ad
#SaveH5Seurat(spatial_cells, filename = "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial//spatial_cells1.h5Seurat")
#Convert("/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/control/spatial_cells1.h5Seurat", dest = "h5ad")

#坐标存储
write.csv(coord, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/model/coord2_.csv")
          
 #control
stdata<- Load10X_Spatial(
  data.dir = "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/control/",  
  filename = "filtered_feature_bc_matrix.h5"
)
scdata <- readRDS("/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/scRNA/merged__scRNAseq2.rds")

#随后安装官网运行指令进行CellTrek分析
stdata <- RenameCells(stdata, new.names=make.names(Cells(stdata)))
scdata <- RenameCells(scdata, new.names=make.names(Cells(scdata)))


stdata <- SCTransform(stdata, assay = "Spatial", verbose = FALSE)
stdata <- RunPCA(stdata, assay = "SCT", verbose = FALSE)
stdata <- FindNeighbors(stdata, reduction = "pca", dims = 1:30)
stdata <- FindClusters(stdata, verbose = FALSE)
stdata <- RunUMAP(stdata, reduction = "pca", dims = 1:30)

# 现在重新画
SpatialDimPlot(stdata)

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
a<-c
traint <- CellTrek::traint(st_data=stdata, sc_data=scdata, sc_assay='RNA', cell_names='cell_type')
DimPlot(traint, group.by = "type") 
cat("🟢 正在开始 celltrek 映射...\n")
start_time <- Sys.time()
celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data=scdata, sc_assay = 'RNA', reduction='pca', intp=T, intp_pnt=1000, intp_lin=F, nPCs=15, ntree=300, dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T, num.threads = 16)$celltrek
cat("✅ celltrek 映射完成\n")
print(Sys.time() - start_time)
celltrek$cell_type <- factor(celltrek$seurat_clusters, levels=sort(unique(celltrek$seurat_clusters)))

#接卷积结果的可视化，可以不运行
#CellTrek::celltrek_vis(celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_type:id_new),
#                       celltrek@images$D1.AD12.L1@image, celltrek@images$D1.AD12.L1@scale.factors$lowres)


#获取解卷积的细胞以及其坐标
coord <- celltrek@meta.data[,c("coord_x","coord_y")]
# 获取解卷积出的细胞 ID（字符向量）
spatial_ids <- celltrek@meta.data$id_raw
# 保证这些细胞 ID 在 scdata 里确实存在
valid_spatial_ids <- intersect(spatial_ids, colnames(scdata))
# 仅提取 scdata 中存在的细胞
spatial_cells <- scdata[, valid_spatial_ids]

dim(spatial_cells)


spatial_cells@assays$RNA@data <- spatial_cells@assays$RNA@counts
write.csv(spatial_cells@assays$RNA@data, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/control/spatial__adata2.csv")

#Seurat对象转为h5ad
#SaveH5Seurat(spatial_cells, filename = "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial//spatial_cells1.h5Seurat")
#Convert("/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/control/spatial_cells1.h5Seurat", dest = "h5ad")

#坐标存储
write.csv(coord, "/slurm/home/yrd/fanlab/gaomeng/NRI/spatial/control/coord2_.csv")
