#### Modules


############################# Select smORF module

SelectSmorfUI=function(id){
  ns=NS(id)
  tagList(selectInput(ns("smORF_ID"),label = "Select a smORF ID",choices = Annotation$iORF_id)
  )
}

SelectSmorfServer=function(id){
  moduleServer(id,function(input,output,session){
    selected_smORF_ID=reactive({input$smORF_ID})
    return(selected_smORF_ID)
  })
}

############################## Home module
HomeTabUI=function(id){
  ns=NS(id)
  tabPanel(h5("Home"),
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
                        " and are yet to be characterised."),
                 style="text-align: justify;"),
               br(),
               p(paste0("In order to facilitate their function characterisation, we develop smORF-pro. ",
                        "smORF-pro is in an in silico platform that gather deep analysis of ",
                        "transcription profiles, coexpression patterns, conservations and protein motifs",
                        "; with the purpose of helping our community to generate hypothesis on smORF functions."),
                 style="text-align: justify;"),
               )
           )
  )
}

HometabServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    })
}


############################## Summary module
SmorfSummaryUI=function(id){
  ns=NS(id)
  tabPanel(h5("Summary"),
           titlePanel("smORF Summary"),
           textOutput(ns("smORF_summary"))
  )
}

SmorfSummaryServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$smORF_summary=renderText({
      load_smORF(selected_smORF())
      smORF_object$Name
    })
  })
}


############################## Conservation module
ConservationUI=function(id){
  ns=NS(id)
  tabPanel(h5("Conservation"),
           titlePanel("Alignment plot"),
           plotOutput(ns("ConservationAlignment")),
           titlePanel("Table of percentage conservation per species"),
           tableOutput(ns("ConservationTable")),
           titlePanel("Conservation Tree"),
           plotOutput(ns("ConservationTree"))
  )
}

ConservationServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$ConservationAlignment=renderPlot({
      load_smORF(selected_smORF())
      smORF_object[["Conservation"]]
      })
  })
}

############################## Tissue specificity Module
TissueSpecificityUI=function(id){
  ns=NS(id)
  tabPanel(h5(HTML(paste("Tissue", "Specificity",sep = "<br/>")),align='center'),
           titlePanel("Anatogram"),
           plotOutput(ns("TissueAnatogram")),
           titlePanel("Tissue TPM expression"),
           plotOutput(ns("TissueTPM"))
  )
}

TissueSpecificityServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$TissueAnatogram=renderPlot({
      load_smORF(selected_smORF())
      ggplot()
    })
    output$TissueTPM=renderPlot({
      load_smORF(selected_smORF())
      plot(smORF_object[["TissueTPM"]])},
      height = 900
    )
  })
}

################################ Cell type specificity
CellTypeUI=function(id){
  ns=NS(id)
  tabPanel(h5(HTML(paste("Cell Type", "Specificity",sep = "<br/>")),align='center'),
           titlePanel("Single-cell data"),
           selectInput(ns("SelectedTissue"),label = "Select a tissue",choices = c("Heart","Brain","Blood")),
           plotOutput("SingleCellUMAP"),
           titlePanel("Expression of selected gene"),
           htmlOutput(ns("GeneDescription")),
           selectInput(ns("GeneID"),label="",choices = Annotation$Gene_id),#selected = smORF_object$geneID),
           plotOutput(ns("CellTypePlot"))
  )
}

CellTypeServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$SingleCellUMAP=renderPlot({
    ggplot()
    })
    output$GeneDescription=renderUI({
      load_smORF(selected_smORF())
      HTML(CellTypeText(smORF_object))
    })
    output$CellTypePlot=renderPlot({
      load_smORF(selected_smORF())
      ggplot()
    })
  })
}

################################ Gene ontology tab. Bubble_GO for RNA and RIBO
GeneOntologyUI=function(id){
  ns=NS(id)
  tabPanel(h5(HTML(paste("Pathway", "Signatures",sep = "<br/>")),align='center'),
           selectInput(ns("gene_ID"),label = "Select an Ensembl gene ID",selected = "ENSG00000226688",choices=All_genes_files),
           textOutput(ns("geneID")),
           titlePanel("RNA-seq Gene Ontology"),
           #uiOutput(ns("geneID_ui")),
           plotOutput(ns("RNA_GO"),height = "100%"),
           titlePanel("RIBO-seq Gene Ontology"),
           plotOutput(ns("RIBO_GO"))
  )
}

GeneOntologyServer=function(id,selected_smORF){
   moduleServer(id,function(input,output,session){
  #   output$geneID_ui=renderUI({
  #     load_smORF(selected_smORF())
  #     selectInput("gene_ID",label = "Select an Ensembl gene ID",selected = smORF_object$geneID,choices=Annotation$Gene_id)
  #   })
    output$geneID=renderText({
      load_smORF(selected_smORF())
      paste0(smORF_object$Name," is located on gene ",smORF_object$geneID)
    })
    output$RNA_GO=renderPlot({
      geneID_i=input$gene_ID
      bubble_GO(geneID_i)
    },height=1000)
    output$RIBO_GO=renderPlot({
      load_smORF(selected_smORF())
      ggplot()
    })
  })
}

################################ Tissue Signature clustering
SignatureClusterUI=function(id){
  ns=NS(id)
  tabPanel(h5(HTML(paste("Signature", "Clustering",sep = "<br/>")),align='center'),
           textOutput(ns("SmorfGene")),
           titlePanel("Tissue Clusters"),
           selectInput(ns("Tissue"),label = "Select a GTex Tissue",choices = c("Heart","Brain","Blood")),
           plotOutput(ns("SignatureUMAP")),
           selectInput(ns("geneID"),label = "Select an Ensembl gene ID",choices = Annotation$Gene_id),
           textOutput(ns("GeneCluster")),
           selectInput(ns("ClusterNumber"),label="Select a cluster of interest",choices = 0:50),
           fluidRow(
             column(6,"Pathway Signatures enriched in this cluster:",
                    tableOutput(ns("ClusterSignatureTable"))),
             column(6,"Genes present in this cluster:",
                    tableOutput(ns("GenesInClusterTable")))
           )
  )
}

SignatureClusterServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$SmorfGene=renderText({
      load_smORF(selected_smORF())
      paste0(smORF_object$Name," is located on gene ",smORF_object$geneID)
    })
    output$SignatureUMAP=renderPlot({
      load_smORF(selected_smORF()) # will have to change to load_UMAP_tissue(input$Tissue)
      ggplot()
    })
    output$GeneCluster=renderText({
      paste0(input$geneID," is part of cluster X")
    })
    output$ClusterSignatureTable=renderTable({

    })
    output$GenesInClusterTable=renderTable({

    })
  })
}

################################ AA based analysis, Deeploc, TMHMM, InterProScan
AABasedAnalysisUI=function(id){
  ns=NS(id)
  tabPanel(h5(HTML(paste("AA Sequence", "Analysis",sep = "<br/>")),align='center'),
           titlePanel("Subcellular Localisation"),
           fluidRow(
             column(width=3,h3("Deeploc Table"),
                    htmlOutput(ns("DeeplocText"))),
             column(7,h3("Cell Anatogram"),#style = "background-color:#eeeeee;",
                    plotOutput(ns("CellAnatogram"))),
             column(width=1,h3(""),
                    plotOutput(ns("CellAnatogramScale"))),
           ),
           titlePanel("TMHMM"),
           titlePanel("SignalP"),
           titlePanel("InterProScan")
  )
}

AABasedAnalysisServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$DeeplocText=renderUI({
      load_smORF(selected_smORF())
      HTML(Deeploc_text(smORF_object))
    })
    output$CellAnatogram=renderPlot({
      load_smORF(selected_smORF())
      plot_cell_sublocation(deeploc_i = smORF_object$Deeploc)
    },height = 500)
    output$CellAnatogramScale=renderPlot({
      load_smORF(selected_smORF())
      plot_cell_sublocation_scale()
    },height = 200)
  })
}

################################ Mutation tab
MutationUI=function(id){
  ns=NS(id)
  tabPanel(h5(HTML(paste("Known", "Mutations",sep = "<br/>")),align='center'),
           tableOutput(ns("MutationTable"))
  )
}

MutationServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$MutationTable=renderTable({
      load_smORF(selected_smORF())
      data.frame()
    })
  })
}


################################ UTR Correlation
UTRCorrelationUI=function(id){
  ns=NS(id)
  tabPanel(h5(HTML(paste("UTR Correlation", "With Host",sep = "<br/>")),align='center')
           #titlePanel("Template_section_name"),
           #plotOutput(ns("TemplatePlot"))
  )
}

UTRCorrelationServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
#    output$TemplatePlot=renderPlot({
#      load_smORF(selected_smORF())
#      ggplot()
#    })
  })
}



################################ Template Module, 1 module = 1 tab
TemplateUI=function(id){
  ns=NS(id)
  tabPanel("Template_tab_name",
           titlePanel("Template_section_name"),
           plotOutput(ns("TemplatePlot"))
  )
}

TemplateServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$TemplatePlot=renderPlot({
      load_smORF(selected_smORF())
      ggplot()
    })
  })
}









