# Function to load the smORF of interest
load_smORF=function(smORF_ID="smORF_ID"){
  if (!exists("smORF_object")) {
    file_name=paste0(smORF_ID,".rds")
    smORF_selected=readRDS(file_name)
    assign("smORF_object",smORF_selected,envir = .GlobalEnv)
  }
  if (exists("smORF_object")&smORF_object[["Name"]]!=smORF_ID) {
    file_name=paste0(smORF_ID,".rds")
    smORF_selected=readRDS(file_name)
    assign("smORF_object",smORF_selected,envir = .GlobalEnv)
  }
}

# Function to generate text in the cell type specificity tab
CellTypeText=function(smORF_object){
  Text=paste0(
    smORF_object$Name,
    " is a smorf located on gene ",
    smORF_object$geneID,
    "<br/>","<br/>",
    "Enter any Ensembl gene ID to display in the UMAP and violin plot"
  )
  return(Text)
}

# Function to make the text from deeploc output

Deeploc_text=function(smORF_object){
  deeploc_i=smORF_object$Deeploc
  Text=paste0(
    "<br/>",
    "Deeploc gives a prediction score from 0 (low confidence) to 1 (high confidence). For ",
    smORF_object$Name,
    " the predicted score per organelle are as follow:",
    "<br/>","<br/>",
    "<br/>",
    paste0(colnames(deeploc_i)[c(1:5,7:10,6)],": ",round(deeploc_i[1,c(1:5,7:10,6)],3),collapse = "<br/>"),
    "<br/>",
    "<br/>","<br/>",
    "The most confident subcellular localization is ",
    colnames(deeploc_i)[which(deeploc_i[1,]==max(deeploc_i[1,]))],
    ".",
    "<br/>","<br/>","<br/>"
  )
  return(Text)
}


# Function to plot cell image with organelles according to deeploc output
plot_cell_sublocation=function(deeploc_i){
  binned_colors=colorRampPalette(c("white","blue"))(11)
  deeploc_i=round(deeploc_i,1)
  for (i in 1:ncol(deeploc_i)) {
    cell_img[organelles[[i]]]=binned_colors[10*deeploc_i[1,i]+1]
  }
  return(plot(cell_img))
}


# Plot the color scale of the cell sublocation plot
plot_cell_sublocation_scale=function(){
  binned_colors=colorRampPalette(c("white","blue"))(11)
  df_scale=data.frame(x=rep(1,11),y=rep(0.1,11),colors=binned_colors[11:1])
  ggplot(df_scale,aes(x=x,y=y,fill=colors))+
    scale_fill_manual(values = df_scale$colors)+geom_bar(stat="identity")+
    ylab("Deeploc Score")+
    theme(legend.position = "none",
          axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_rect(fill = "white",
                                          colour = "white"))+
    scale_y_continuous(limits = c(0,1), expand = c(0, 0),breaks = seq(0,1,0.2))+
    scale_x_continuous(limits = c(0.5,1.5), expand = c(0, 0))
}

# Read tree





