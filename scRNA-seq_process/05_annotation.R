library(Seurat)
library(dplyr)

library(future)
options(future.fork.enable = TRUE) ## In case you're in RStudio
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^5)
plans()

Idents(alldata) <- 'scvi_snn_res.1'
markers <- FindAllMarkers(alldata, only.pos = TRUE)

# Festem
library(Festem)

myfestem <- function(obj,by){
  alldata1 <- obj
  alldata1[["RNA"]] <- as(object = alldata1[["RNA"]], Class = "Assay5")
  G <- length(unique(alldata1@meta.data[[by]]))
  alldata1 <- RunFestem(alldata1,batch = 'orig.ident',num.threads=80, prior = by,G = G)
  df.marker2 <- AllocateMarker(alldata1,VariableFeatures(alldata1),num_cores = 80,group_by=by)
  return(df.marker2)
}

df_festem <- myfestem(alldata,'scvi_snn_res.1')

alldata$celltype_level1 <- recode(alldata$scvi_snn_res.1,
                                '1'='Neutrophils',
                                '2'='B',
                                '3'='Neutrophils',
                                '4'='T',
                                '5'='T',
                                '6'='Endothelial',
                                '7'='Epithelial',
                                '8'='Macrophages',
                                '9'='B',
                                '10'='B',
                                '11'='T',
                                '12'='Monocytes',
                                '13'='Macrophages',
                                '14'='Epithelial',
                                '15'='Epithelial',
                                '16'='T',
                                '17'='NK',
                                '18'='T',
                                '19'='Epithelial',
                                '20'='Monocytes',
                                '21'='Macrophages',
                                '22'='AF',
                                '23'='Neutrophils',
                                '24'='Macrophages',
                                '25'='Endothelial',
                                '26'='AF',
                                '27'='DC',
                                '28'='AT2',
                                '29'='Epithelial',
                                '30'='B',
                                '31'='Neutrophils',
                                '32'='Club',
                                '33'='Plasma',
                                '34'='DC',
                                '35'='Club',
                                '36'='AF',
                                '37'='Monocytes',
                                '38'='Endothelial',
                                '39'='Macrophages',
                                '40'='AT1',
                                '41'='Neutrophils',
                                '42'='SMC',
                                '43'='DC',
                                '44'='B',
                                '45'='AF',
                                '46'='Macrophages',
                                '47'='Neutrophils',
                                '48'='B',
                                '49'='Epithelial',
                                '50'='Neutrophils')

qs::qsave(alldata,'alldata_filtered.qs')