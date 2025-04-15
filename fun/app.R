setwd("~/smORF-pro/")
#### Load libraries ####

# library(Seurat)
library(shiny)
library(bslib)
library(shinydashboard)
# library(data.table)
library(ggplot2)
# library(colorspace)
# library(gridExtra)
# library(reshape)
# library(ggmsa)
# library(png)
# library(raster)
# library(WGCNA)
# library(fgsea)
library(stringr)
library(shinythemes)
library(DT)
library(gganatogram)
library(seqinr)

#### Load data ####
Annotation=readRDS("data/Annotation_full.rds")#[1:100,]
df_DT=readRDS("data/Annotation_DT.rds")

RiboAllMerged=readRDS("data/Ribo_All_gsea_padj.rds")
RiboAllMergedNES=readRDS("data/Ribo_All_gsea_NES.rds")

gene_names=read.csv("data/gene_names.csv")

GO_cluster=readRDS("data/GO_clusters.rds")
markers=readRDS("data_preparation/Expression/Coexpression/Ribo/GO_clusters/Cluster_markers_All.rds")

P_to_run1=readRDS("data/P_to_run1.rds")
P_to_run2=readRDS("data/P_to_run2.rds")

load_smORF("ENCT00000053446_ncRNA_3")
GWAS_results=read.csv("data/GWAS_smORF_overlap.csv")[,-1]


ui=page_navbar(title="smORF-pro",
               id="navbar",theme = bs_theme(preset = "flatly"),
               nav_panel("Home",
                         HomeUI("Home"),
                         layout_columns(actionButton("smORF",label = "smORF"),
                                        actionButton("Explore",label = "Explore"),
                                        col_widths = c(-3,3,3,-3))),
               nav_panel("smORF",
                         smorfUI("smORF")),
               nav_panel("Explore",
                         ExploreUI("Explore"))
               )

server=function(input,output,session){
  observeEvent(input$smORF,{
    nav_select("navbar",selected = "smORF")
  })
  observeEvent(input$Explore,{
    nav_select("navbar",selected = "Explore")
  })
  HomeServer("Home")
  smorfServer("smORF")
  ExploreServer("Explore")
}

shinyApp(ui, server)







