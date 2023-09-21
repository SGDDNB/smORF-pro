## Preparation of cell image

library(png)
library(raster)
library(colorspace)
library(ggplot2)


cell_img=readPNG("Cell organelle.png")
cell_img=as.raster(cell_img)
#plot(cell_img)
#View(as.data.frame(table(cell_img)))

organelles=list()
organelles[["Cytoplasm"]]=which(cell_img=="#FFE791FF")
organelles[["Nucleus"]]=which(cell_img=="#BFD5FFFF")
organelles[["Extracellular"]]=which(cell_img=="#FFFF00FF")
organelles[["Cell Membrane"]]=which(cell_img=="#00E6FFFF")
organelles[["Mitochondrion"]]=which(cell_img=="#F4A76AFF")
organelles[["Mitochondrion"]]=c(organelles[["Mitochondrion"]],which(cell_img=="#F4A769FF"))
organelles[["Plastid"]]=which(cell_img=="")
organelles[["Endoplasmic reticulum"]]=which(cell_img=="#7D3705FF")
organelles[["Lyosome/Vacuole"]]=which(cell_img=="#D66370FF")
organelles[["Lyosome/Vacuole"]]=c(organelles[["Lyosome/Vacuole"]],which(cell_img=="#7300FFFF"))
organelles[["Golgi apparatus"]]=which(cell_img=="#8CC63FFF")
organelles[["Peroxisome"]]=which(cell_img=="#FF00FFFF")

saveRDS(organelles,"Organelles.rds")

# cell_img[organelles$Cytoplasm]="#00000000"
# cell_img[organelles$Nucleus]="#00000000"
# cell_img[organelles$Extracellular]="#00000000"
# cell_img[organelles$`Cell Membrane`]="#00000000"
# cell_img[organelles$Mitochondrion]="#00000000"
# cell_img[organelles$Plastid]="#00000000"
# cell_img[organelles$`Endoplasmic reticulum`]="#00000000"
# cell_img[organelles$`Lyosome/Vacuole`]="#00000000"
# cell_img[organelles$`Golgi apparatus`]="#00000000"
# cell_img[organelles$Peroxisome]="#00000000"
# plot(cell_img)


## Tree
library(ape)
library(stringr)
library(phylogram)
library(ggtree)
library(tidyverse)
library(rphylopic)
library(phytools)
library(deeptime)
library(ggimage)

data("vertebrate.tree")

Species_tree=read.tree("hg38.100way.commonNames.nh")
Species_tree$tip.label=gsub("00000"," ",Species_tree$tip.label)
Species_tree$tip.label=gsub("AAAAA","'",Species_tree$tip.label)

Species_group=data.frame(Groups=c("fish","Euarchontoglires","primates",
                                  "Laurasiatheria"),
                         uuid=c("ba7d2620-6e6d-4985-9de6-6951a03606b6",
                                "88a07585-846a-405d-9195-c15c010e7443",
                                "2d078b25-e6a0-4beb-a5d3-5d6f16be8ebf",
                                "f72f8145-1860-407b-838a-cff909a6bd69"))
Species_group$svg= lapply(Species_group$uuid, get_phylopic)
Species_group$x=c(-9,-5.5,-7,-8)
Species_group$y=c(11,98,85,71)
Species_group$size=c(1.5,2,2,1)

ggtree(Species_tree,layout="dendrogram",branch.length = "none")+geom_tiplab(align = T,)+
  hexpand(-0.8)+geom_tippoint(aes(1:199),color=c("black",rep("white",198)))
  add_phylopic(Species_group$svg,
               x=Species_group$x,
               y=Species_group$y,
               ysize=Species_group$size)+
  geom_highlight(node = c(100,186),fill="blue")+
  geom_highlight(node = 123,fill="green")+
  geom_highlight(node = 111,fill="orange")+
  geom_highlight(node = 135,fill="yellow")+
  geom_highlight(node = 160,fill="purple")+
  geom_nodepoint(aes(subset=node==111),size=5)

Species_tree$edge[,2]=""

phylopic_info=vertebrate_data
phylopic_info$species=1:11



data(vertebrate.tree)
vertebrate_data <- data.frame(species = vertebrate.tree$tip.label, uuid = NA)
# Try to get PhyloPic UUIDs for the species names
vertebrate_data$uuid <- sapply(vertebrate.tree$tip.label,
                               function(x) {
                                 tryCatch(get_uuid(x), error = function(e) NA)
                               })
vertebrate_data
vertebrate_data$uuid[vertebrate_data$species == "Myotis_lucifugus"] <-
  get_uuid("Vespertilioninae")
vertebrate_data$svg <- lapply(vertebrate_data$uuid, get_phylopic)


# Load pathways
pathways_to_use = c("c2.cp.kegg","c5.mf","c5.bp","c5.cc","h.all")
pathways_list = list()
path_gsea_annot="Coexpression/gsea_pathways/genename_annotations_downloaded_MSigDB_v6.2/"
for(pathway in pathways_to_use){
  pathway_anno = gmtPathways(paste(path_gsea_annot,pathway,".v6.2.symbols.gmt",sep=""))
  pathways_list[[pathway]] = pathway_anno
  Details=read.csv(paste(path_gsea_annot,pathway,".csv",sep = ""), header = FALSE)
  #  MitoMatrix[[pathway]]=as.data.frame(Details$V1)
}
P=pathways_list
P_to_run=c(P$c2.cp.kegg,P$c5.mf,P$c5.bp,P$c5.cc,P$h.all)
saveRDS(P_to_run,"P_ro_run.rds")





