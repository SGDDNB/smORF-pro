
# load libraries
library(ggplot2)
library(data.table)
library(stringr)
library(dplyr)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(shiny)
library(ggmsa)
library(Seurat)
library(gridExtra)
library(Peptides)
library(WGCNA)
library(fgsea)
library(Matrix)
library(reshape2)
library(presto)

# setwd("/home/baptiste/smORF-pro/")
setwd("/home/baptiste/Desktop/PhD/smorfs functionalization/")
# setwd("/Users/baptiste/Desktop/smorfs functionalization/")
# setwd("C:/Users/e0545037/Desktop/Baptiste/PhD/smorfs_functionalization/")

# Load functions

source("2.smORFs_func.R")

# Load data per tab
Annotation=read.csv("data/Annotation.csv")
#Annotation_new_smorfs=read.csv("Include new smORFs/new_smorfs_annotation.csv")
#colnames(Annotation_new_smorfs)=colnames(Annotation)
#Annotation=rbind(Annotation,Annotation_new_smorfs)
load("data/DEG_list.Rdata")
load("data/DEG_smORFs_summary.Rdata")
ConservationTable=read.csv("data/Conservation_table.csv")
interproscan=read.csv("data/Interproscan_smORFs.csv")
#interproscan_new_smorfs=fread("Include new smORFs/new_smorfs_interproscan.tsv")
#colnames(interproscan_new_smorfs)=colnames(interproscan)
#interproscan=rbind(interproscan,interproscan_new_smorfs)
load("data/tpm_data.Rdata")

#percentageID=fread("data/protein_percentage_ids.txt")
#percentageID[which(is.na(percentageID),arr.ind = T)]=0

# speciesNames=fread("data/species_names.txt")
# percentageID$V1=gsub('"',"",percentageID$V1)
# percentageID=as.data.frame(percentageID)
# rownames(percentageID)=percentageID$V1
# percentageID=percentageID[,-1]
# colnames(percentageID)=speciesNames$COMMON[match(colnames(percentageID),speciesNames$GENOME)]
# percentageID=percentageID[,c("Human","Rat","Mouse","Zebrafish")]

GTex_Annotation=fread("Coexpression/GTEx_Analysis_v8_Annotations_SampleAttributesDS.csv")
GTex=fread("Coexpression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
GTex$Name=gsub("\\..*","",GTex$Name)
GTex=as.data.frame(GTex)
GTex=GTex[!duplicated(GTex$Name),]
rownames(GTex)=GTex$Name


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


sc_obj=readRDS("Coexpression//score_brain-hippocampus_sc_obj.rds")
#sc_obj=readRDS("score_Cells-Culturedfibroblasts_GO_obj.rds")
markers=wilcoxauc(sc_obj)

# ui with tabs

ui <- navbarPage("smORF-pro",
                 tabPanel("Home",
                          sidebarLayout(
                            sidebarPanel(
                              h2("Documentation"),
                              p("smORF-pro is a shiny app developped on Rstudio. The code can be found on github at:",
                                style="text-align: justify;"),
                              a("www.github/SGDDNB/smORF-pro/",href="www.github/SGDDNB/smORF-pro/"),
                              br(),
                              br(),
                              p("(Link to paper)"),
                              br(),
                              br(),
                              div(img(src="github.png",height=100,width=100), style="text-align: center;"),
                            ),
                            mainPanel(
                              h1("Introducing smORF-pro"),
                              p(paste0("small Open Reading Frames (smORFs), are proteins that are shorter than 100 amino acids. ",
                              "Due to their small size, they have remained understudied until the development of technologies",
                              " like ribosome-profiling data. We believe that several thousands of them are expressed in human",
                              "and are yet to be characterised."),
                              style="text-align: justify;"),
                              br(),
                              p(paste0("In order to facilitate their function characterisation, we develop smORF-pro. ",
                                       "smORF-pro is in an in silico platform that gather deep analysis of ",
                                       "transcription profiles, coexpression patterns, conservations and protein motifs",
                                       "; with the purpose of helping our community to generate hypothesis on smORF functions."),
                                style="text-align: justify;"),
                              ))),
                 tabPanel("Differentially expressed smORFs",
                          plotOutput("p1",height = "100%"),
                          plotOutput("p2",height = "100%"),
                          sliderInput("threshold","Remove smORFs with no |LFC| above threshold:",
                                      value = 9,min = 0,max = 10,step = 0.5),
                          plotOutput("p3",height = "100%"),
                          selectInput("Experiment_name", label = "Experiment",
                                      choices = names(DEG)),
                          plotOutput("p4",height = "100%")),
                 tabPanel("Annotation",
                          dataTableOutput("Annotation"),
                          downloadButton("smorfs_download","Download database"),
                          plotOutput("p5",height = "100%"),
                          plotOutput("p6",height = "100%"),
                          textInput("smORF_ID","ORF ID"),
                          plotOutput("p7",height = "100%"),
                          plotOutput("p8",height = "100%")),
                 #tabPanel("Coexpression",
                          #selectInput("Tissue_GTex",label = "Tissue of interest",
                          #            choices=unique(GTex_Annotation$SMTSD)),
                          #selectInput("iORF_ID1",label = "smORF ID",
                          #            choices=Annotation$iORF_id[Annotation$Gene_id%in%GTex$Name]),
                          #dataTableOutput("CoExpression_i"),
                          #dataTableOutput("GSEA_i")),
                 tabPanel("Genetic mutations"),
                 tabPanel("Protein-Protein interaction"),
                 tabPanel("Gene Ontology",
                          selectInput("Tissue_GO",label="Select Tissue",
                                      choices=c("Brain-Hyppocampus")),
                          titlePanel("T-SNE visualization"),
                          plotOutput("TSNE_plot",width = "100%",height = "100%"),
                          selectInput("gene_go",label = "GeneID",
                                      choices = colnames(sc_obj)),
                          textOutput("text_GO_cluster"),
                          titlePanel("Cluster details"),
                          selectInput("cluster_nb",label = "Cluster number",
                                      choices = 0:24),
                          dataTableOutput("Cluster_head_table")),
                 tabPanel("Cell specificity"),
                 tabPanel("smORFs conservation",
                          dataTableOutput("ConservationTable")),
                 tabPanel("smORF card",
                          selectInput("iORF_ID",label = "smORF ID",
                                      choices=Annotation$iORF_id),
                          textOutput("text1"),
                          titlePanel("smORF conservation plot"),
                          plotOutput("pcard1",width = "100%",height = "100%"),
                          dataTableOutput("FewSpecies"),
                          titlePanel("smORF hydrophobicity"),
                          plotOutput("hydro",width = "100%",height = "100%"),
                          titlePanel("Annotation of smORF"),
                          dataTableOutput("Annotation_i"),
                          titlePanel("Interproscan detailed table"),
                          dataTableOutput("interproscan"),
                          titlePanel("RNA and RIBO TPM in tissues"),
                          plotOutput("smORFs_tpm_plots",height = "100%",width = "100%"),
                          titlePanel("Gene ontology signature"),
                          plotOutput("smORF_bubble_GO",height = "100%",width = "100%")))

server <- function(input, output, session) {
  # Plot how many DE smORFs per experiment
  output$p1=renderPlot(Plot_DE_smorfs_per_exp(DEG_summary$padj),height = 500,width = 900)

  # Plot smORFs DE in how many experiments
  output$p2=renderPlot(Plot_DE_smorfs_many_exp(DEG_summary$per_exp),height = 500,width = 900)

  #Plot heatmap
  output$p3=renderPlot({
    Thres=input$threshold
    Plot_heatmap(Threshold = Thres)},height = 800)

  # Plot volcano plot
  output$p4=renderPlot(Plot_volcano(input$Experiment_name),height = 500,width = 900)

  # Output the smORFs table
  output$Annotation=renderDataTable(Annotation)




  # Coexpression tab ####

  # output$CoExpression_i=renderDataTable({
  #   Tissue=input$Tissue_GTex
  #   ORF=input$iORF_ID1
  #   ORF_Gene_ID=Annotation$Gene_id[Annotation$iORF_id==ORF]
  #   GTex_counts=GTex[,colnames(GTex)%in%c("Name",GTex_Annotation$SAMPID[GTex_Annotation$SMTSD==Tissue])]
  #   rownames(GTex_counts)=make.names(names = GTex_counts$Name,unique = T)
  #   GTex_counts=GTex_counts[,-1]
  #   GTex_counts=GTex_counts[,GTex_counts[ORF_Gene_ID,]>1]
  #   coexpr=corAndPvalue(t(GTex_counts),t(GTex_counts[ORF_Gene_ID,]))
  #   coexpr=data.frame(Gene_ID=rownames(GTex_counts),cor=coexpr$cor,pval=coexpr$p)
  #   colnames(coexpr)=c("Gene_ID","Correlation","pvalue")
  #   coexpr
  # })

  output$GSEA_i=renderDataTable({
    # Tissue=input$Tissue_GTex
    # ORF=input$iORF_ID1
    # ORF_Gene_ID=Annotation$Gene_id[Annotation$iORF_id==ORF]
    # GTex_counts=GTex[,colnames(GTex)%in%c("Name",GTex_Annotation$SAMPID[GTex_Annotation$SMTSD==Tissue])]
    # rownames(GTex_counts)=make.names(names = GTex_counts$Name,unique = T)
    # GTex_counts=GTex_counts[,-1]
    # GTex_counts=GTex_counts[,GTex_counts[ORF_Gene_ID,]>1]
    # coexpr=corAndPvalue(t(GTex_counts),t(GTex_counts[ORF_Gene_ID,]))
    # coexpr=data.frame(Gene_ID=rownames(GTex_counts),cor=coexpr$cor,pval=coexpr$p)
    # colnames(coexpr)=c("Gene_ID","Correlation","pvalue")
    # coexpr

    coexpr=coexpr[order(coexpr$Correlation,decreasing = T),]
    ranked_list=data.frame(geneID=coexpr$Gene_ID,corr=coexpr$Correlation)
    ranked_list$geneName=NA
    ranked_list$geneName=Genes_ID_Names$geneName[match(ranked_list$geneID,Genes_ID_Names$geneID)]
    ranked_list=ranked_list[!is.na(ranked_list$corr),]
    Rank_forgsea=ranked_list$corr
    names(Rank_forgsea)=ranked_list$geneName
    fgseaRes <- fgsea(pathways = P_to_run, stats = Rank_forgsea, minSize=10, maxSize=1000)
    fgseaRes
  })

  # Gene ontology tab

  output$TSNE_plot=renderPlot({
    DimPlot(sc_obj,reduction="tsne",label = T)},height = 400,width = 600)

  output$text_GO_cluster=renderText({
    Gene=input$gene_go
    cluster_nb=Idents(sc_obj)[Gene]
    paste("This gene is part of cluster",cluster_nb,sep = " ")})

    output$Cluster_head_table=renderDataTable({
      cluster_nb=input$cluster_nb
      m_i=markers[which(markers$group==cluster_nb),]
      m_i=m_i[order(m_i$logFC,decreasing = T),]
      m_i[1:15,]
      })






  # Output Conservation table
  output$ConservationTable=renderDataTable(ConservationTable)

  #### smorf card tab ####
  observe({
    output$pcard1=renderPlot({
      ORF=input$iORF_ID
      Plot_conservation(SEP = ORF)},
      height = nrow(fread(paste0("data/Conservation_fasta/",input$iORF_ID,".fa")))*20,
      width = nchar(Annotation$Peptide.seq[which(Annotation$iORF_id==input$iORF_ID)])*40)
  })

  # output$FewSpecies=renderDataTable({
  #   ORF=input$iORF_ID
  #   percentageID[rownames(percentageID)==ORF,]
  # })

  # observe({
  #   output$hydro=renderPlot({
  #     ORF=input$iORF_ID
  #     plot_hydro(SEP = ORF)},
  #     height = 500,width = 700)
  # })

  output$text1=renderText({
    ORF=input$iORF_ID
    smORF_text(ORF)})

  output$Annotation_i=renderDataTable({
    ORF=input$iORF_ID
    Annotation[Annotation$iORF_id==ORF,]
  })

  output$interproscan=renderDataTable({
    ORF=input$iORF_ID
    interproscan[interproscan$iORF_ID==ORF,]
  })

  output$smORFs_tpm_plots=renderPlot({
    ORF=input$iORF_ID
    smORF_tpm_plot(ORF)},height = 500,width = 900)

  output$smORF_bubble_GO=renderPlot({
    ORF=input$iORF_ID
    Gene=Annotation$Gene_id[which(Annotation$iORF_id==ORF)]
    bubble_GO(Gene=Gene)},height = 900,width = 1200)




  # output$p5=renderPlot(ggplot(TMHMM_sum,aes(x=Var1,y=Freq,label=Freq))+
  #                        geom_bar(stat = "identity",fill="dodgerblue")+
  #                        theme_bw()+ylab("Number of smORFs")+xlab("Protein domain")+ggtitle("TMHMM prediction")+
  #                        geom_text(vjust=-0.5),height = 500,width = 500)
  # output$p6=renderPlot(ggplot(TargetP_sum,aes(x=Var1,y=Freq,label=Freq))+
  #                        geom_bar(stat = "identity",fill="dodgerblue")+
  #                        theme_bw()+ylab("Number of smORFs")+xlab("Prediction")+ggtitle("TargetP prediction")+
  #                        geom_text(vjust=-0.5),height = 500,width = 500)

}


shinyApp(ui, server)


#shiny::runApp("/home/baptiste/smORF-pro/Main_app/app.R")









