# 10x PAH NucSeq - CellChat Processing
# Yvon Mbouamboua

# Setting parameters
options(future.globals.maxSize = 80000*1024^2)
setwd("/data/data_mbouamboua/projects/10x_htap_nucseq/analysis/")
set.seed(4743)

# Required packages
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(Seurat)
  library(scater)
  library(ggplot2)
  library(patchwork)
  library(NMF)
  library(CellChat)
  library(ggalluvial)
  library(ComplexHeatmap)
})

# Prepare data
htap <- readRDS("htap.filt.rds")
DefaultAssay(htap) <- "RNA"

cell_order <- c("Basal","MCC", "PNEC","Pre-TB-SC","TRB-SC","AT0", "AT2a","AT2b","AT1-AT2","AT1-imm","AT1",
                "aCap", "gCap", "ArtEC", "V-PulmEC", "V-SysEC","Lymph",
                "AdvFibro","AlvFibro","MyoFibro","ASMC","VSMC","Peri",
                "MacroATP10A","MacroCD68","Macro-prolif","DC1","DC2","Mono","Inter-Macro","CD4","CD8","T-prolif","NK","Mast","Plasma","pDC","B")

# Create a CTRL CellChat object
invisible(gc())
sub <- subset(htap, subset = group == "CTRL")
cells.include <- Cells(sub)[!is.na(sub$celltype)]
meta = data.frame(labels = as.factor(sub$celltype)[cells.include]) 
counts = GetAssayData(sub, assay = "RNA",slot = "data")[,cells.include] 
cellchat <- createCellChat(object = counts, meta = meta, group.by = "labels")
cellchat@idents <- factor(cellchat@idents,levels = cell_order)
cellchat <- setIdent(cellchat, ident.use = "labels") 
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use 
cellchat <- subsetData(cellchat)  # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
invisible(gc())
cellchat <- projectData(cellchat, PPI.human)
cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
cellchat <- computeCommunProb(cellchat)
#cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP") 
selectK(cellchat, pattern = "outgoing")
invisible(gc())
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, width = 4, height = 9)
invisible(gc())
selectK(cellchat, pattern = "incoming")
invisible(gc())
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, width = 4, height = 9)
invisible(gc())
saveRDS(cellchat, "cellchat.CTRL.rds")


# Create a HTAP CellChat object
sub <- subset(htap, subset = group == "PAH")
cells.include <- Cells(sub)[!is.na(sub$celltype)]
meta = data.frame(labels = as.factor(sub$celltype)[cells.include]) 
counts = GetAssayData(sub, assay = "RNA",slot = "data")[,cells.include] 
cellchat <- createCellChat(object = counts, meta = meta, group.by = "labels")
cellchat@idents <- factor(cellchat@idents,levels = cell_order)
cellchat <- setIdent(cellchat, ident.use = "labels") 
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use 
cellchat <- subsetData(cellchat)  
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
invisible(gc())
cellchat <- projectData(cellchat, PPI.human)
cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
cellchat <- computeCommunProb(cellchat)
#cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP") 
selectK(cellchat, pattern = "outgoing")
invisible(gc())
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, width = 4, height = 9)
invisible(gc())
selectK(cellchat, pattern = "incoming")
invisible(gc())
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, width = 4, height = 9)
invisible(gc())
saveRDS(cellchat, "cellchat.PAH.rds")

# Session info
utils::capture.output(devtools::session_info())


