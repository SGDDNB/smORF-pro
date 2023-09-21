library(shiny)
library(shinydashboard)
library(data.table)
library(ggplot2)
library(colorspace)
library(gridExtra)
library(reshape)
library(ggmsa)
library(png)
library(raster)
library(WGCNA)
library(fgsea)
library(stringr)

#####
setwd("/home/baptiste/Desktop/PhD/smorfs functionalization")

##############################

Annotation=read.csv("data/Annotation.csv")
cell_img=as.raster(readPNG("Cell organelle.png"))
organelles=readRDS("Organelles.rds")
P_to_run=readRDS("P_to_run.rds")
clean_pathways=readRDS("clean_pathways.rds")
All_genes_files=list.files("Coexpression/NES_per_gene/")
All_genes_files=str_sub(All_genes_files,1,nchar(All_genes_files)-4)

# page for smORF ID
ui=fluidPage(
  navbarPage(h4("smORF-pro"),
             header=SelectSmorfUI("smORF_ID"),
             HomeTabUI("Home"),
             SmorfSummaryUI("Summary"),
             ConservationUI("Conservation"),
             TissueSpecificityUI("Tissue"),
             CellTypeUI("CellType"),
             GeneOntologyUI("GO"),
             SignatureClusterUI("SignatureCluster"),
             AABasedAnalysisUI("AABasedAnalysis"),
             MutationUI("Mutation"),
             UTRCorrelationUI("UTRCorrelation")
  )
)

server=function(input,output,session){
  selected_smORF=reactive(SelectSmorfServer("smORF_ID"))
  HomeTabUI("Home")
  SmorfSummaryServer("Summary",selected_smORF())
  ConservationServer("Conservation",selected_smORF())
  TissueSpecificityServer("Tissue",selected_smORF())
  CellTypeServer("CellType",selected_smORF())
  GeneOntologyServer("GO",selected_smORF())
  SignatureClusterServer("SignatureCluster",selected_smORF())
  AABasedAnalysisServer("AABasedAnalysis",selected_smORF())
  MutationServer("Mutation",selected_smORF())
  UTRCorrelationServer("UTRCorrelation",selected_smORF())
}

shinyApp(ui,server)


#shiny::runApp("/home/baptiste/smORF-pro/Main.R")







