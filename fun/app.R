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
library(NGLVieweR)

#### Load data ####
Annotation=readRDS("data/Annotation_full.rds")
df_DT=Annotation

source("fun/fun.R",local = T)
source("fun/modules.R",local = T)
source("fun/structure.R",local = T)

Explore_smORFs_ID=readRDS("data/Explore_smORF_ID_names_for_pathways.rds")

code_names=read.csv("data/code_names.csv")
gene_names=read.csv("data/gene_names.csv")

GO_cluster=readRDS("data/GO_clusters.rds")
markers=readRDS("data/Cluster_markers_All.rds")

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
                                        actionButton("Find",label = "Find"),
                                        col_widths = c(4,4,4))),
               nav_panel("smORF",
                         smorfUI("smORF")),
               nav_panel("Explore",
                         ExploreUI("Explore")),
               nav_panel("Find",
                         FindUI("Find"))
               )

server=function(input,output,session){
  observeEvent(input$smORF,{
    nav_select("navbar",selected = "smORF")
  })
  observeEvent(input$Explore,{
    nav_select("navbar",selected = "Explore")
  })
  observeEvent(input$Find,{
    nav_select("navbar",selected = "Find")
  })
  HomeServer("Home")
  smorfServer("smORF",session)
  ExploreServer("Explore",session)
  FindServer("Find",session)
}

shinyApp(ui, server)