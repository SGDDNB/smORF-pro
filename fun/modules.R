##### smORF Selection #####

SelectSmorfUI=function(id){
  ns=NS(id)
  selectInput(ns("smORF_ID"),label = "Select a smORF ID",choices = Annotation$iORF_ID,selected = smORF_object$Annotation$iORF_ID)
}

SelectSmorfServer=function(id){
  moduleServer(id,function(input,output,session){
    selected_smORF_ID=reactive({input$smORF_ID})
    return(selected_smORF_ID)
  })
}

# short smORF description

smORF_descrUI=function(id){
  ns=NS(id)
  textOutput(ns("smORF_descr"))
}

smORF_descrServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$smORF_descr=renderText({
      load_smORF(selected_smORF())
      paste0(smORF_object$Annotation$iORF_ID," is a ", smORF_object$Annotation$ORF_type, " protein located on gene ",
                   smORF_object$Annotation$Gene_id, " - ", smORF_object$Annotation$Gene_name,
             " and was found in ", Annotation$Source[which(Annotation$iORF_ID==smORF_object$Annotation$iORF_ID)])
    })
  })
}


##### smORF Page Modules #####

##### smORF Summary tab #####

SmorfSummaryUI=function(id){
  ns=NS(id)
  fluidPage(fluidRow(h3("Protein predicted structure"),
    column(
      4, NGLVieweROutput(ns("ngl"), height = "500px")
    ),
    column(
      4,
      selectInput(ns("style"), "Representation",
                  c("cartoon", "ribbon", "trace", "line (licorice)" = "licorice", "spacefill"),
                  selected = "cartoon"),
      selectInput(ns("color"), "Color by",
                  c("B-factor (pLDDT)" = "bfactor",
                    "Residue index" = "residueindex",
                    "Uniform" = "uniform",
                    "Hydrophobicity"="hydrophobicity"),
                  selected = "bfactor"),
      checkboxInput(ns("spin"),"Spin",F)
    )
  ),fluidRow(
    column(width=4,h3("Summary"),
           htmlOutput(ns("smORF_summary"))),
    column(width=4,h3("Conservation"),
           plotOutput(ns("RefPlot")))
  ),fluidRow(
    column(width=4,h3("Expression"),
           plotOutput(ns("RiboTPM"))),
    column(width=4,h3("Ontology signature"),
           htmlOutput(ns("GSEA_pathways")))
  ))
}

SmorfSummaryServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    ns=session$ns

    output$ngl <- renderNGLVieweR({
      load_smORF(selected_smORF())
      code_i <- code_names$Code[code_names$Sequence == smORF_object$Annotation$Peptide_Seq]
      cif_path <- paste0("data_preparation/Structure/", code_i, "/", code_i, "_model_0.cif")

      NGLVieweR(cif_path) %>%
        stageParameters(list(showBackground = TRUE,
                             backgroundColor = "white",
                             backgroundOpacity = 0.8)) %>%
        addRepresentation(input$style, list(colorScheme = input$color)) %>%
        setFocus(0) %>%
        setSpin(FALSE)
    })

    proxy <- NGLVieweR_proxy("ngl", session = session)

    observeEvent(list(input$style, input$color, selected_smORF()), {
      output$ngl <- renderNGLVieweR({
        load_smORF(selected_smORF())
        code_i <- code_names$Code[code_names$Sequence == smORF_object$Annotation$Peptide_Seq]
        cif_path <- paste0("data_preparation/Structure/", code_i, "/", code_i, "_model_0.cif")

        NGLVieweR(cif_path) %>%
          stageParameters(list(showBackground = TRUE,
                               backgroundColor = "white",
                               backgroundOpacity = 0.8)) %>%
          addRepresentation(input$style, list(colorScheme = input$color)) %>%
          setFocus(0) %>%
          setSpin(FALSE)
      })
    }, ignoreInit = TRUE)

    observeEvent(input$spin, ignoreInit = TRUE, {
      proxy %>% updateSpin(input$spin)
    })


    output$smORF_summary=renderUI({
      load_smORF(selected_smORF())
      HTML(smORF_summary_text(smORF_object))
    })

    output$RefPlot=renderPlot({
      load_smORF(selected_smORF())
      if (length(smORF_object$Conservation)>1) {
        plot_ref(smORF=smORF_object)
      } else {
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("MSA failed to align this\n","smORF in other species."),
             cex = 2.2, col = "red")
      }
    },height = 300,width = 350)

    output$RiboTPM=renderPlot({
      load_smORF(selected_smORF())
      ggplot(smORF_object$TissueTPM,aes(x=Tissue,y=RIBO))+geom_boxplot(fill="darkorange2")+
        theme_bw()+xlab("")+ylab("TPM")+ggtitle("RIBO expression")+
        theme(plot.title =element_text(hjust = 0.5,size=15),
              axis.text.x = element_text(angle = 45,hjust=1,size = 10))},
      height = 400)

    output$GSEA_pathways=renderUI({
      load_smORF(selected_smORF())
      HTML(GO_summary_text(smORF_object))
      })
  })
}

##### AA Analysis #####

DomainUI=function(id){
  ns=NS(id)
  tableOutput(ns("DomainTable"))
}

DomainServer=function(id,selected_smORF,parent_session){
  moduleServer(id,function(input,output,session){
    output$DomainTable=renderTable({
      load_smORF(selected_smORF())
      smORF_object$Annotation
    })
  })
}

SubcellularUI=function(id){
  ns=NS(id)
  fluidRow(
    column(width=3,h3("Score per organelle"),
           htmlOutput(ns("DeeplocText"))),
    column(width=7,h3("Cell Anatogram"),
           uiOutput(ns("CellAnatogram"))),
    column(width=1,
           plotOutput(ns("CellAnatogramScale"))),
  )
}

SubcellularServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$DeeplocText=renderUI({
      load_smORF(selected_smORF())
      if (nchar(smORF_object$Annotation$Peptide_Seq)>9) {
        HTML(Deeploc_text(smORF_object))
      }
    })
    output$CellAnatogram=renderUI({
      load_smORF(selected_smORF())
      if (nchar(smORF_object$Annotation$Peptide_Seq)>9) {
        tags$iframe(style="height:500px; width:700px",src=paste0("Deeploc/",smORF_object$Annotation$iORF_ID,".png"))
      } else {
        br()
        h2("This smORF is too short to be analyzed by Deeploc")
      }
    })
    output$CellAnatogramScale=renderPlot({
      load_smORF(selected_smORF())
      if (nchar(smORF_object$Annotation$Peptide_Seq)>9) {
        plot_cell_sublocation_scale()
      }
    },height = 200,width = 80)
  })
}

#### Genetic Information #####

ConservationUI=function(id){
  ns=NS(id)
  page_fluid(
    titlePanel("Reference species"),
    plotOutput(ns("RefPlot")),
    titlePanel("Alignment plot"),
    uiOutput(ns("ConservationAlignment"))
  )
}


ConservationServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    ns=session$ns
    output$RefPlot=renderPlot({
      load_smORF(selected_smORF())
      if (length(smORF_object$Conservation)>1) {
        plot_ref(smORF=smORF_object)
      } else {
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("MSA failed to align this\n","smORF in other species."),
             cex = 2.2, col = "red")
      }
    },height = 300,width = 350)

    output$ConservationAlignment=renderUI({
      load_smORF(selected_smORF())
      if(paste0(smORF_object$Annotation$iORF_ID,".png")%in%list.files("www/Conservation/")){
        heighti=20+30*smORF_object$Conservation$nb_species
        widthi = 50+30*smORF_object$Conservation$max_L
        tags$iframe(style=paste0("height:",heighti,"px; width:",widthi,"px"),
                    src=paste0("Conservation/",smORF_object$Annotation$iORF_ID,".png"))
      }
    })
  })
}

GeneticVariantUI=function(id){
  ns=NS(id)
  page_fluid(dataTableOutput(ns("variants")),
             textOutput(ns("NoMutation")))
}

GeneticVariantServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$variants=renderDataTable({
      load_smORF(selected_smORF())
      if (smORF_object$Annotation$iORF_ID%in%GWAS_results$iORF_id) {
        FindMutation(smORF_object)
      }
    })
      if (!smORF_object$Annotation$iORF_ID%in%GWAS_results$iORF_id) {
        output$NoMutation=renderText({
           paste0("There is no GWAS hit overlapping with the CDS of this smORF from GWAS catalog.")
        })
      }
  })
}


##### Expression tab #####

TissueSpecificityUI=function(id){
  ns=NS(id)
  page_fluid(
    fluidRow(
      column(4,
             titlePanel("RNAseq TPM"),
             plotOutput(ns("RNAseqPlot"))),
      column(4,
             titlePanel("RIBOseq TPM"),
             plotOutput(ns("RIBOseqPlot")))),
    titlePanel("GTex database RNA Expression"),
    uiOutput(ns("RNAseqGTexPlot"),height = 700),
    fluidRow(
      column(4,
             titlePanel("Male RNA TPM Anatogram"),
             plotOutput(ns("MaleAnatogram"),height = "100%")),
      column(4,
             titlePanel("Female RNA TPM Anatogram"),
             plotOutput(ns("FemaleAnatogram"),height = "100%")),
      )
  )
}

TissueSpecificityServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$RNAseqPlot=renderPlot({
      load_smORF(selected_smORF())
      ggplot(smORF_object$TissueTPM,aes(x=Tissue,y=RNA))+geom_boxplot(fill="deepskyblue1")+
          theme_bw()+xlab("")+ylab("TPM")+ggtitle("RNA expression")+
          theme(plot.title =element_text(hjust = 0.5,size=15),
                axis.text.x = element_text(angle = 45,hjust=1))},
      height = 400)
    output$RIBOseqPlot=renderPlot({
      load_smORF(selected_smORF())
      ggplot(smORF_object$TissueTPM,aes(x=Tissue,y=RIBO))+geom_boxplot(fill="darkorange2")+
          theme_bw()+xlab("")+ylab("TPM")+ggtitle("RIBO expression")+
          theme(plot.title =element_text(hjust = 0.5,size=15),
                axis.text.x = element_text(angle = 45,hjust=1))},
      height = 400)
    output$RNAseqGTexPlot=renderUI({
      load_smORF(selected_smORF())
      if (paste0(smORF_object$Annotation$Gene_id,".png")%in%list.files("www/GTex_Expression")) {
        tags$iframe(style="height:700px; width:1800px",src=paste0("GTex_Expression/",smORF_object$Annotation$Gene_id,".png"))
      } else {
        h2("This smORF is on a gene that is not annotated in GTex database")
      }
    })
    output$MaleAnatogram=renderPlot({
      load_smORF(selected_smORF())
      if (!is.na(smORF_object$GTex_mean[1,2])) {
        Male_anatogram(smORF_object$GTex_mean)
      }
    },height = 600,width = 400)
    output$FemaleAnatogram=renderPlot({
      load_smORF(selected_smORF())
      if (!is.na(smORF_object$GTex_mean[1,2])) {
        Female_anatogram(smORF_object$GTex_mean)
      }
    },height = 600,width = 400)
  })
}

CellUI=function(id){
  ns=NS(id)
  page_fluid(
    selectInput(ns("scTissue"),"Select a Tissue",choices = c("Heart")),
    fluidRow(
      column(7,
             h3("UMAP clustering"),
             uiOutput(ns("dimplot"))
             ),
      column(5,
             uiOutput(ns("Feature"))
      )
    ),
    uiOutput(ns("Violin"))
  )
}

CellServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$dimplot=renderUI({
      tags$iframe(style="height:500px; width:1100px",src="heart_sc.png")
    })
    output$Feature=renderUI({
      load_smORF(selected_smORF())
      tags$iframe(style="height:500px; width:600px",src=paste0("Heart/",smORF_object$Annotation$Gene_name,"_F.png"))
    })
    output$Violin=renderUI({
      load_smORF(selected_smORF())
      tags$iframe(style="height:500px; width:100%",src=paste0("Heart/",smORF_object$Annotation$Gene_name,"_V.png"))
    })

  })
}

##### Correlation & Ontology #####

CorrelationUI=function(id){
  ns=NS(id)
  page_fluid(
    titlePanel("smORF Ribo Coexpression"),
    dataTableOutput(ns("CoexpressionTable"))
  )
}

CorrelationServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$CoexpressionTable=renderDataTable({
      load_smORF(selected_smORF())
      if (nrow(smORF_object$Coexpr)>0) {
        smORF_object$Coexpr$Gene_Name=gene_names$IDENTIFIER[match(smORF_object$Coexpr$ID,gene_names$geneID)]
        smORF_object$Coexpr=smORF_object$Coexpr[,c(1,4,2,3)]
        datatable(rownames = F,
                  smORF_object$Coexpr,
                  filter = list(position = 'top', clear = FALSE),
                  extensions = "Buttons",
                  options = list(paging = TRUE,
                                 scrollX=TRUE,
                                 searching = TRUE,
                                 ordering = TRUE,
                                 dom = 'lBtip',
                                 buttons = c('copy', 'csv', 'excel', 'pdf'),
                                 pageLength=10,
                                 lengthMenu=c(10,20,50,100) )
        )
      }
    })
  })
}

gseaUI=function(id){
  ns=NS(id)
  page_fluid(
    fluidRow(
             column(3,selectInput(ns("nb_pathway"),"Number of top pathways per tissue",choices = 1:10,
                                  selected = 5))),
    titlePanel("GTex RNAseq GSEA"),
    plotOutput(ns("GTex_gsea"),height = 800),
    fluidRow(
      column(6,
             titlePanel("RNAseq GSEA"),
             plotOutput(ns("RNAgsea"))
             ),
      column(6,
             titlePanel("RIBOseq GSEA"),
             plotOutput(ns("RIBOgsea"))
             )
      ),
    uiOutput(ns("UTRtitle")),
    uiOutput(ns("UTR"))
  )
}

gseaServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$GTex_gsea=renderPlot({
      load_smORF(selected_smORF())
      if (length(smORF_object$GSEA$GTex)>0) {
        Bubble_simple(smORF_object$GSEA$GTex,input$nb_pathway)
      }
    },height = 800)

    output$RNAgsea=renderPlot({
      load_smORF(selected_smORF())
      if (length(smORF_object$GSEA$RNA)>0) {
        Bubble_simple(smORF_object$GSEA$RNA,input$nb_pathway)
      }
    })
    output$RIBOgsea=renderPlot({
      load_smORF(selected_smORF())
      if (length(smORF_object$GSEA$Ribo)>0) {
        Bubble_simple(smORF_object$GSEA$Ribo,input$nb_pathway)
      }
    })

    output$UTRtitle=renderUI({
      load_smORF(selected_smORF())
      if (length(smORF_object$GSEA$UTR)>0) {
        titlePanel("Comparison smORF vs mORF")
      }
    })

    output$UTR=renderUI({
      load_smORF(selected_smORF())
      if (length(smORF_object$GSEA$UTR)>0) {
        ns=NS(id)
        plotOutput(ns("UTRgsea"))
        output$UTRgsea=renderPlot({
          Bubble_simple(smORF_object$GSEA$UTR,input$nb_pathway)
        },width = 800)
      }
    })
  })
}


GOClusterUI=function(id){
  ns=NS(id)
  page_fluid(
    titlePanel("Ribo-seq GSEA All Tissue clustering"),
    #selectInput(ns("ID_go"),"Select a gene or smORF",choices = colnames(GO_cluster)),
    textOutput(ns("cluster_ID_go")),
    selectInput(ns("cluster_nb"),"Select a cluster", choices = 0:20),
    uiOutput(ns("umap")),
    fluidRow(
      column(4,
             titlePanel("Proteins in selected cluster"),
             dataTableOutput(ns("cluster_gene_list"))),
      column(8,
             titlePanel("Pathways enriched in this cluster"),
             dataTableOutput(ns("cluster_markers")))
    )
  )
}

GOClusterServer=function(id,selected_smORF){
  moduleServer(id,function(input,output,session){
    output$cluster_ID_go=renderText({
      load_smORF(selected_smORF())
      cluster_ID_i=as.numeric(GO_cluster$Cluster[which(GO_cluster$geneID==smORF_object$Annotation$Gene_id)])
      if (smORF_object$Annotation$Gene_id%in%GO_cluster$geneID) {
        paste0(smORF_object$Annotation$iORF_ID," is part of cluster ",cluster_ID_i)
      }
      else {
        paste0(smORF_object$Annotation$iORF_ID,"'s gene was not detected in this scRNAseq.")
      }
    })
    output$umap=renderUI({
      tags$iframe(style="height:500px; width:700px",src="GO_cluster.png")
      })
    output$cluster_gene_list=renderDataTable({
      datatable(rownames = F,
                gene_names[which(gene_names$geneID%in%GO_cluster$geneID[which(GO_cluster$Cluster==input$cluster_nb)]),],
                filter = list(position = 'top', clear = FALSE),
                extensions = "Buttons",
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'lBtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf'),
                               pageLength=10,
                               lengthMenu=c(10,20,50,100)
                               )
                )
      })
    output$cluster_markers=renderDataTable({
      markers_i=markers[which(markers$group==input$cluster_nb),]
      markers_i=markers_i[order(markers_i$padj),]
      datatable(rownames = F,
                markers_i,
                filter = list(position = 'top', clear = FALSE),
                extensions = "Buttons",
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'lBtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf'),
                               pageLength=10,
                               lengthMenu=c(10,20,50,100)
                )
      )
    })
  })
}

##### Explore Page modules #####

ExploreTableUI=function(id){
  ns=NS(id)
  page_fluid(
    h3("Explore smORFs from Annotation"),
    p(""),
    p("Click on a smORF to jump to its individual page."),
    dataTableOutput(ns("ExploreTable")),
    br(),
    downloadButton(ns("Download"),"Download Original Table")
  )
}


ExploreTableServer=function(id,parent_session){
  moduleServer(id,function(input,output,session){

    output$ExploreTable=renderDataTable({
      datatable(rownames = F,
                df_DT,
                filter = list(position = 'top', clear = FALSE),
                extensions = "Buttons",
                escape = F,
                selection = "none",
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'lBtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf'),
                               pageLength=5,
                               lengthMenu=c(5,10,20,50,100) )
      )
    })

    output$Download=downloadHandler(
      filename = function(){"smORFs.csv"},
      content = function(fname){
        write.csv(Annotation, fname)
      }
    )

    proxy <- DT::dataTableProxy("ExploreTable", session = session)

    observeEvent(input$ExploreTable_cell_clicked,{
      if (length(input$ExploreTable_cell_clicked)>0) {
        clicked_row <- input$ExploreTable_cell_clicked$row
        smorf_id <- df_DT[clicked_row, 1]
        load_smORF(smorf_id)
        updateSelectInput(parent_session,"smORF-smORF_ID",selected = smorf_id)
        nav_select("navbar",selected = "smORF",session = parent_session)
      }
    })

    observeEvent(parent_session$input$navbar, {
      if (identical(parent_session$input$navbar, "Explore")) {
        DT::selectRows(proxy, NULL)
      }
    })
  })
}

ExplorePathwaysUI=function(id){
  ns=NS(id)
  fluidPage(
    h3("Explore smORFs from ontology enrichment"),
    p("Click on a smORF to jump to its individual page."),
    selectInput(ns("Pathway"),"Select a pathway",selected = "KEGG_ABC_TRANSPORTERS",
                choices = names(P_to_run2),width = "100%"),
    dataTableOutput(ns("ExplorePathways"))
  )
}


ExplorePathwaysServer=function(id,parent_session){
  moduleServer(id,function(input,output,session){

    output$ExplorePathways=renderDataTable({
      datatable(rownames = F,
                                  data.frame(cbind(Explore_smORFs_ID,readRDS(paste0("data/Explore_pathways/",
                                                                                    input$Pathway,".rds")))),
                                  filter = list(position = 'top', clear = FALSE),
                                  extensions = "Buttons",
                                  options = list(paging = TRUE,
                                                 scrollX=TRUE,
                                                 ordering = TRUE,
                                                 dom = 'lBtip',
                                                 buttons = c('copy', 'csv', 'excel', 'pdf'),
                                                 pageLength=10,
                                                 lengthMenu=c(10,20,50,100)
                                  ))
    })

    proxy <- DT::dataTableProxy("ExplorePathways", session = session)

    observeEvent(input$ExplorePathways_cell_clicked,{
      if (length(input$ExplorePathways_cell_clicked)>0) {
        clicked_row <- input$ExplorePathways_cell_clicked$row
        smorf_id <- Explore_smORFs_ID[clicked_row, 1]
        load_smORF(smorf_id)
        updateSelectInput(parent_session,"smORF-smORF_ID",selected = smorf_id)
        nav_select("navbar",selected = "smORF",session = parent_session)
      }
    })

    observeEvent(parent_session$input$navbar, {
      if (identical(parent_session$input$navbar, "Explore")) {
        DT::selectRows(proxy, NULL)
      }
    })
  })
}


ExploreMutationsUI=function(id){
  ns=NS(id)
  page_fluid(
    h3("Explore smORFs overlapping with GWAS hits"),
    p("Click on a smORF to jump to its individual page."),
    dataTableOutput(ns("ExploreMutations"))
  )
}

ExploreMutationsServer=function(id,parent_session){
  moduleServer(id,function(input,output,session){
    output$ExploreMutations=renderDataTable({
      datatable(rownames = F,
                GWAS_results,
                filter = list(position = 'top', clear = FALSE),
                extensions = "Buttons",
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               ordering = TRUE,
                               dom = 'lBtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf'),
                               pageLength=10,
                               lengthMenu=c(10,20,50,100) )
      )
    })

    proxy <- DT::dataTableProxy("ExploreMutations", session = session)

    observeEvent(input$ExploreMutations_cell_clicked,{
      if (length(input$ExploreMutations_cell_clicked)>0) {
        clicked_row <- input$ExploreMutations_cell_clicked$row
        smorf_id <- GWAS_results$iORF_id[clicked_row]
        load_smORF(smorf_id)
        updateSelectInput(parent_session,"smORF-smORF_ID",selected = smorf_id)
        nav_select("navbar",selected = "smORF",session = parent_session)
      }
    })

    observeEvent(parent_session$input$navbar, {
      if (identical(parent_session$input$navbar, "Explore")) {
        DT::selectRows(proxy, NULL)
      }
    })
  })
}


#### Find Page ####


FindTableUI=function(id){
  ns=NS(id)
  fluidPage(
    h3("Find your smORF"),
    p("Click on a smORF to jump to its individual page."),
    dataTableOutput(ns("FindTable")),
  )
}

FindTableServer=function(id,parent_session){
  moduleServer(id,function(input,output,session){
    output$FindTable=DT::renderDataTable({
      datatable(df_DT[,c(1,2,3,7)],
                filter = list(position = 'top', clear = FALSE),
                extensions = "Buttons",
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'lBtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf'),
                               pageLength=10,
                               lengthMenu=c(5,10,20,50,100)),
                rownames = F
    )},selection = list(mode="single",target="row"))

    proxy <- DT::dataTableProxy("FindTable", session = session)

    observeEvent(input$FindTable_cell_clicked,{
      if (length(input$FindTable_cell_clicked)>0) {
        clicked_row <- input$FindTable_cell_clicked$row
        smorf_id <- df_DT[clicked_row, 1]
        load_smORF(smorf_id)
        updateSelectInput(parent_session,"smORF-smORF_ID",selected = smorf_id)
        nav_select("navbar",selected = "smORF",session = parent_session)
      }
    })

    observeEvent(parent_session$input$navbar, {
      if (identical(parent_session$input$navbar, "Find")) {
        DT::selectRows(proxy, NULL)
      }
    })

  })
}






















