#### HomePage ####

HomeUI=function(id){
  ns=NS(id)
  page_fluid(sidebarLayout(sidebarPanel(
    h2("Documentation"),
    p("smORF-pro is a shiny app developped on Rstudio. The code can be found on github at:",
      style = "text-align: justify;"),
    br(),
    br(),
    p("(Link to paper)"),
    br(),
    br(),
    div(img(src = "github.png",height = 100,width = 100), style = "text-align: center;"),),
    mainPanel(h1("Introducing smORF-pro"),
              p(paste0("small Open Reading Frames (smORFs), are proteins that are shorter than 100 amino acids. ",
                       "Due to their small size, they have remained understudied until the development of technologies",
                       " like ribosome-profiling data. We believe that several thousands of them are expressed in human",
                       " and are yet to be characterised."),style = "text-align: justify;"),
              br(),
              p(paste0("In order to facilitate their function characterisation, we develop smORF-pro. ",
                       "smORF-pro is in an in silico platform that gather deep analysis of ",
                       "transcription profiles, coexpression patterns, conservations and protein motifs",
                       "; with the purpose of helping our community to generate hypothesis on smORF functions.",
                       " The smORFs  used in smORF-pro have been identified from one of those 3 studies: Chothani",
                       " et al. Mol Cell. 2022, van Heesch et al. Cell 2019 and Zhang et al. Nat Comms 2020."),
                style = "text-align: justify;"),
              br(),
              br()
    )
  ),
  br(),
  br(),
  fluidRow(
    column(4,
           titlePanel(h3("Individual smORF analysis")),
           p(paste0("Select a smORF ID to visualize all the analysis done for that smORF and to help you generate",
                    " hypothesis on its biological function."),
             style = "text-align: justify;"),
           align="center"
    ),
    column(4,
           titlePanel(h3("Explore smORFs database")),
           p(paste0("Shortlist which smORFs are of interest for you by filtering through the different features:",
                    " Subcellular localization, Conservation, PepScore, Ontology, ..."),
             style = "text-align: justify;"),
           align="center"
    ),
    column(4,
           titlePanel(h3("Find your smORF")),
           p(paste0("Find your smORF based of gene ID, gene IDENTIFIER or protein sequence."),
             style = "text-align: justify;"),
           align="center"
    ))
  )
}




HomeServer = function(id) {
  moduleServer(id, function(input, output, session) {
  })
}

#### smORF Page modules ####

smorfUI <- function(id) {
  ns <- NS(id)
  tagList(selectInput(ns("smORF_ID"),label = "Select a smORF ID",choices = Annotation$iORF_ID,selected = smORF_object$Annotation$iORF_ID),
          smORF_descrUI(ns("smORF_descr")),
  navset_card_tab(id=ns("tabs"),
                  nav_panel("Summary", SummarySmORFText(),
                            SmorfSummaryUI(ns("Summary"))),
                  nav_menu("Amino Acid\n Analysis",
                           nav_panel("Protein Domain", ProteinDomainText(),
                                     DomainUI(ns("Domain"))),
                           nav_panel("Subcellular Localisation", ProteinDeeplocText(),
                                     SubcellularUI(ns("Deeploc")))),
                  nav_menu("Genetic\n Information",
                           nav_panel("Conservation", ConservationText(),
                                     ConservationUI(ns("Conservation"))),
                           nav_panel("Genetic Variant", GeneticVariantText(),
                                     GeneticVariantUI(ns("GWAS")))),
                  nav_menu("Tissue and Cell\n Expression",
                           nav_panel("RNAseq & Riboseq\n Specificity", ExpressionText(),
                                     TissueSpecificityUI(ns("Tissue"))),
                           nav_panel("scRNAseq Specificity",SingleCellText(),
                                     CellUI(ns("scRNAseq")))),
                  nav_menu("Correlation\n Analysis",
                           nav_panel("Correlation", CoexpressionText(),
                                     CorrelationUI(ns("corr"))),
                           nav_panel("GSEA",gseaText(),
                                     gseaUI(ns("GSEA"))),
                           nav_panel("Ontology Clustering", clusterText(),
                                     GOClusterUI(ns("clustering"))))
  ))
}

smorfServer=function(id,parent_session){
  moduleServer(id, function(input,output,session){
    selected_smORF_ID=reactive({input$smORF_ID})
    smORF_descrServer("smORF_descr",selected_smORF_ID)
    SmorfSummaryServer("Summary",selected_smORF_ID)
    DomainServer("Domain",selected_smORF_ID,parent_session)
    SubcellularServer("Deeploc",selected_smORF_ID)
    ConservationServer("Conservation",selected_smORF_ID)
    GeneticVariantServer("GWAS",selected_smORF_ID)
    TissueSpecificityServer("Tissue",selected_smORF_ID)
    CellServer("scRNAseq",selected_smORF_ID)
    CorrelationServer("corr",selected_smORF_ID)
    gseaServer("GSEA",selected_smORF_ID)
    GOClusterServer("clustering",selected_smORF_ID)
  })
}


#### Explore Page ####

ExploreUI=function(id){
  page_fluid(
    navset_card_tab(id="Explore_tabs",
      nav_panel("Annotation",
                ExploreTableUI("ExploreTable")),
      nav_panel("Mutations",
                ExploreMutationsUI("ExploreMutations")),
      nav_panel("Pathways",
                ExplorePathwaysUI("ExplorePathway")))
  )
}


ExploreServer=function(id,parent_session){
    ExploreTableServer("ExploreTable",parent_session)
    ExplorePathwaysServer("ExplorePathway",parent_session)
    ExploreMutationsServer("ExploreMutations",parent_session)
}


#### Find Page ####


FindUI=function(id){
  page_fluid(
    FindTableUI("FindTable")
  )
}


FindServer=function(id,parent_session){
  FindTableServer("FindTable",parent_session)
}


