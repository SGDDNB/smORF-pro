library(gganatogram)

# Function to load the smORF of interest
load_smORF=function(smORF_ID="smORF_ID"){
  if (!exists("smORF_object")) {
    file_name=paste0(smORF_ID,".rds")
    smORF_selected=readRDS(paste0("~/smORF-pro/smORF_objects/",file_name))
    assign("smORF_object",smORF_selected,envir = .GlobalEnv)
  }
  if (exists("smORF_object")&smORF_object$Annotation$iORF_ID!=smORF_ID) {
    file_name=paste0(smORF_ID,".rds")
    smORF_selected=readRDS(paste0("~/smORF-pro/smORF_objects/",file_name))
    assign("smORF_object",smORF_selected,envir = .GlobalEnv)
  }
}

Deeploc_text=function(smORF_object){
  deeploc_i=smORF_object$Deeploc
  Text=paste0(
    "<br/>",
    "Deeploc gives a prediction score from 0 (low confidence) to 1 (high confidence). For ",
    smORF_object$Annotation$iORF_ID,
    " the predicted scores per organelle are as follow:",
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

plot_ref=function(smORF){
  df=data.frame(Species=factor(c("Human","Chimp","Rat","Mouse","Zebrafish"),
                               levels = c("Human","Chimp","Rat","Mouse","Zebrafish")),
                Conservation=unlist(smORF$Conservation$df[1,]))
  ggplot(df,aes(x=Species,y=Conservation,label=Conservation,fill=T))+geom_bar(stat = "identity")+
    geom_text(vjust=-0.25)+
    theme_bw()+scale_fill_manual(values="navy")+
    guides(fill="none")+
    theme(text = element_text(size = 15))
}


Male_key=hgMale_key[c(20,58,47,6,31,52,34,43,30,37,44,8,24,65,45,38,48,68,41,26,25,29,67,1,13),]
Male_key$colour="white"

Male_anatogram=function(GTex_mean){
  pal=colorRamp(c("white","deepskyblue"))
  Male_key$value[which(Male_key$organ=="spinal_cord")]=GTex_mean$RNA_TPM[which(GTex_mean$Tissue=="Brain - Spinal cord (cervical c-1)")]
  Male_key$value[which(Male_key$organ=="esophagus")]=mean(GTex_mean$RNA_TPM[grep("Esophagus",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="lung")]=GTex_mean$RNA_TPM[which(GTex_mean$Tissue=="Lung")]
  Male_key$value[which(Male_key$organ=="heart")]=mean(GTex_mean$RNA_TPM[grep("Heart",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="adipose_tissue")]=mean(GTex_mean$RNA_TPM[grep("Adipose",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="aorta")]=mean(GTex_mean$RNA_TPM[grep("Aorta",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="small_intestine")]=mean(GTex_mean$RNA_TPM[grep("Small Intestine",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="coronary_artery")]=mean(GTex_mean$RNA_TPM[grep("Coronary",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="urinary_bladder")]=mean(GTex_mean$RNA_TPM[grep("Bladder",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="brain")]=mean(GTex_mean$RNA_TPM[grep("Brain",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="colon")]=mean(GTex_mean$RNA_TPM[grep("Colon",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="left_ventricle")]=mean(GTex_mean$RNA_TPM[grep("Left Ventricle",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="liver")]=mean(GTex_mean$RNA_TPM[grep("Liver",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="kidney")]=mean(GTex_mean$RNA_TPM[grep("Kidney",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="adrenal_gland")]=mean(GTex_mean$RNA_TPM[grep("Adrenal Gland",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="salivary_gland")]=mean(GTex_mean$RNA_TPM[grep("Salivary",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="skeletal_muscle")]=mean(GTex_mean$RNA_TPM[grep("Muscle",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="pancreas")]=mean(GTex_mean$RNA_TPM[grep("Pancreas",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="prostate")]=mean(GTex_mean$RNA_TPM[grep("Prostate",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="skin")]=mean(GTex_mean$RNA_TPM[grep("Skin",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="spleen")]=mean(GTex_mean$RNA_TPM[grep("Spleen",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="stomach")]=mean(GTex_mean$RNA_TPM[grep("Stomach",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="testis")]=mean(GTex_mean$RNA_TPM[grep("Testis",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="thyroid_gland")]=mean(GTex_mean$RNA_TPM[grep("Thyroid",GTex_mean$Tissue)])
  Male_key$value[which(Male_key$organ=="breast")]=mean(GTex_mean$RNA_TPM[grep("Breast",GTex_mean$Tissue)])
  Male_key$colour=rgb(pal((Male_key$value - min(Male_key$value) )/ diff(range(Male_key$value))), max=255)
  gganatogram(data=Male_key, fillOutline='white', organism='human', sex='male', fill="colour") + theme_void()
}

rownames(hgFemale_key)=1:nrow(hgFemale_key)
Female_key=hgFemale_key[c(11,7,57,26,22,52,16,60,14,19,58,33,13,50,61,66,63,12,55,21,17,54,10,2,24,67,64,68,25),]
Female_key$colour="white"

Female_anatogram=function(GTex_mean){
  pal=colorRamp(c("white","deepskyblue"))
  Female_key$value[which(Female_key$organ=="spinal_cord")]=GTex_mean$RNA_TPM[which(GTex_mean$Tissue=="Brain - Spinal cord (cervical c-1)")]
  Female_key$value[which(Female_key$organ=="esophagus")]=mean(GTex_mean$RNA_TPM[grep("Esophagus",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="lung")]=GTex_mean$RNA_TPM[which(GTex_mean$Tissue=="Lung")]
  Female_key$value[which(Female_key$organ=="heart")]=mean(GTex_mean$RNA_TPM[grep("Heart",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="adipose_tissue")]=mean(GTex_mean$RNA_TPM[grep("Adipose",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="aorta")]=mean(GTex_mean$RNA_TPM[grep("Aorta",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="small_intestine")]=mean(GTex_mean$RNA_TPM[grep("Small Intestine",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="coronary_artery")]=mean(GTex_mean$RNA_TPM[grep("Coronary",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="urinary_bladder")]=mean(GTex_mean$RNA_TPM[grep("Bladder",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="brain")]=mean(GTex_mean$RNA_TPM[grep("Brain",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="colon")]=mean(GTex_mean$RNA_TPM[grep("Colon",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="left_ventricle")]=mean(GTex_mean$RNA_TPM[grep("Left Ventricle",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="liver")]=mean(GTex_mean$RNA_TPM[grep("Liver",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="kidney")]=mean(GTex_mean$RNA_TPM[grep("Kidney",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="adrenal_gland")]=mean(GTex_mean$RNA_TPM[grep("Adrenal Gland",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="salivary_gland")]=mean(GTex_mean$RNA_TPM[grep("Salivary",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="skeletal_muscle")]=mean(GTex_mean$RNA_TPM[grep("Muscle",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="pancreas")]=mean(GTex_mean$RNA_TPM[grep("Pancreas",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="skin")]=mean(GTex_mean$RNA_TPM[grep("Skin",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="spleen")]=mean(GTex_mean$RNA_TPM[grep("Spleen",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="stomach")]=mean(GTex_mean$RNA_TPM[grep("Stomach",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="thyroid_gland")]=mean(GTex_mean$RNA_TPM[grep("Thyroid",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="breast")]=mean(GTex_mean$RNA_TPM[grep("Breast",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="ectocervix")]=mean(GTex_mean$RNA_TPM[grep("Cervix",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="endometrium")]=mean(GTex_mean$RNA_TPM[grep("Cervix",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="fallopian_tube")]=mean(GTex_mean$RNA_TPM[grep("Fallopian",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="ovary")]=mean(GTex_mean$RNA_TPM[grep("Ovary",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="uterus")]=mean(GTex_mean$RNA_TPM[grep("Uterus",GTex_mean$Tissue)])
  Female_key$value[which(Female_key$organ=="vagina")]=mean(GTex_mean$RNA_TPM[grep("Vagina",GTex_mean$Tissue)])
  Female_key$colour=rgb(pal((Female_key$value - min(Female_key$value) )/ diff(range(Female_key$value))), max=255)
  gganatogram(data=Female_key, fillOutline='white', organism='human', sex='female', fill="colour") + theme_void()
}



Bubble_simple=function(GSEA,n){
  padj=GSEA[[1]]
  NES=GSEA[[2]]
  
  top_n=as.data.frame(apply(padj,2,order))
  top_n=unlist(top_n[1:n,])
  top_n=top_n[!duplicated(top_n)]
  top_n=rownames(NES)[top_n]
  
  NES=NES[top_n,]
  padj=padj[top_n,]
  
  score_df=as.data.frame(-log10(padj)*NES)
  score_df=as.data.frame(score_df)
  score_df[which(score_df=="NaN",arr.ind = T)]=0
  score_df[which(score_df=="Inf",arr.ind = T)]=sort(unique(unlist(c(score_df))),decreasing = T)[2]
  
  hclust_score <- hclust(dist(score_df), method = "complete")
  ordre=hclust_score$order
  
  toplot=data.frame(Tissue=rep(colnames(NES),each=nrow(NES)),Pathways=rownames(NES),NES=unlist(NES),padj=as.numeric(unlist(padj)))
  toplot$Pathways=factor(toplot$Pathways,levels=rownames(score_df)[ordre])
  
  toplot=toplot[toplot$Pathways%in%top_n,]
  
  ggplot(toplot,aes(x=Tissue,y=Pathways,size=-log10(as.numeric(padj)),
                    col=as.numeric(NES)))+
    geom_point()+theme_classic()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    scale_color_gradient2()+
    labs(color="NES",size="-log10(Padj)")
}



FindMutation=function(smORF_object){
  nt_seq=smORF_object$nt_seq
  nt_seq=gsub("-","",nt_seq)
  original_AA_seq=c2s(translate(s2c(nt_seq)))
  smORF_i_gtf=smORF_object$gtf

  GWAS_i=GWAS_results[which(GWAS_results$iORF_id==smORF_object$Annotation$iORF_ID),]
  GWAS_i$original_seq=original_AA_seq
  GWAS_i$New_AA_seq=NA
  GWAS_i$Seq_change="No"
  for (i in 1:nrow(GWAS_i)) {
    snp=str_sub(GWAS_i$STRONGEST.SNP.RISK.ALLELE[i],
                start=nchar(GWAS_i$STRONGEST.SNP.RISK.ALLELE[i]),
                end=nchar(GWAS_i$STRONGEST.SNP.RISK.ALLELE[i]))
    if (snp%in%c("A","C","G","T")) {
      pos_snp=GWAS_i$CHR_POS[i]
      if (nrow(smORF_i_gtf)==1) {
        if (smORF_i_gtf$Strand=="+") {
          pos_relative_snp=pos_snp-smORF_i_gtf$Start+1
          substr(nt_seq,pos_relative_snp,pos_relative_snp)=snp
        }
        if (smORF_i_gtf$Strand=="-") {
          pos_relative_snp=smORF_i_gtf$End-pos_snp+1
          substr(nt_seq,pos_relative_snp,pos_relative_snp)=chartr("ATGC","TACG",snp)
        }
      }
      new_seq=c2s(translate(s2c(nt_seq)))
    }
    GWAS_i$New_AA_seq[i]=new_seq
    if (GWAS_i$New_AA_seq[i]!=GWAS_i$original_seq[i]) {
      GWAS_i$Seq_change[i]="Yes"
    }
  }
  return(datatable(GWAS_i[,c(10,13,16,19,23,25,26,27,39,38,52,53,54,45:49)]))
}



ExploreAnnotationText=function(){
  p(paste0("Find annotation of all smORFs into one table. You can filter here to select specific smORFs of interest",
           " or directly download the whole table."),style = "text-align: justify;")
}

ExplorePathwayText=function(){
  p(paste0("Gene Set Enrichment Analysis (GSEA) was run from ranked coexpression analysis for all Ribo-seq tissues. ",
           "Select a pathway to find the NES and padj of smORFs for that pathway."),style = "text-align: justify;")
}

ExploreMutationText=function(){
  p(paste0("smORFs coding sequence (CDS) region coordinates were overlapped with GWAS catalog mutation coordinates.",
           " Here is the table summary of smORFs that had an overlap."),style = "text-align: justify;")
}

SummarySmORFText=function(){
  p(paste0("Summary information of selected smORF. Jump into other tabs for more",
           " details and comprehensive analysis."),style = "text-align: justify;")
}

ProteinDomainText=function(){
  p(paste0("DeepTMHMM, SignalP and InterProScan 2.0 were run to find protein domain within the smORF AA sequence.",
           " This is the one liner output summary of the domains found for the selected smORF."),style = "text-align: justify;")
}

ProteinDeeplocText=function(){
  p(paste0("Deeploc2.0 was run to predict the subcellular localization of the smORF just based on the amino acid sequence.",
           " DeepLoc 2.0 provides a score from 0 (not confident) to 1 (very confident) for every organelle of the cell.",
           " The blue scale represent the confidence score, the deeper the blue is, the highest the score is.")
    ,style = "text-align: justify;")
}

ConservationText=function(){
  p(paste0("Multiple sequence alignment (msa) was run across 100 vertebrates to identify in which species the smORF is conserved",
  " Nucleotide sequences per species were retrieved by using multiple alignment files (MAF). The protein alignment and distance ",
  " calculations were done using Clustal Omega."),style = "text-align: justify;")
}

GeneticVariantText=function(){
  p(paste0("The coding sequence region of this smORF was overlapped with GWAS catalog database to find whether any known mutation",
           "are potentially changing the sequence of the AA sequence of this smORF."),style = "text-align: justify;")
}

ExpressionText=function(){
  p(paste0("smORFs were counted at the gene level in GTex public RNA database to get information on the transcript ",
           "expression in different tissues. smORFs were also counted in in-house and public paired RNA- and Ribo-seq",
           " data to find expression specificity across available tissues."),style = "text-align: justify;")
}

SingleCellText=function(){
  p(paste0("RNA transcript of smORFswere counted in single cell data. scRNAseq data tables are coming from DISCO immune website. ",
           "You can look for more tissues or more specific sample details there."),style = "text-align: justify;")
}

CoexpressionText=function(){
  p(paste0("Spearman protein-protein pairwise correlation was run on Ribo-seq TPM counts of all tissues. ",
           "Here you can find the detailed correlation and pvalue for that specific smORF against all other ",
           "genes and smORFs expressed."),style = "text-align: justify;")
}

gseaText=function(){
  p(paste0("For every tissue in GTex, Gene Set Enrichment Analysis (GSEA) was run for pathways from C2.Kegg, C5 and h.all ",
           "from MSigDb. The same GSEA was also run in our in-house and publically available paired RNA-seq and Ribo-seq.",
           " While RNA may inform you on the transcript correlation signature, Ribo-seq is able to defferentiate reads ",
           "coming from specific parts of the transcript, allowing to count separately uORF/dORFs and their main ORF. ",
           ""),style = "text-align: justify;")
}

clusterText=function(){
  p(paste0("Gene Set Enrichment Analysis (GSEA) was run for every smORFs and genes in the Ribo-seq tissues. The resulting ",
           "NES and padj for every pathway of every smORFs were gathered into 1 table, PCA based clustering was then done to find ",
           "relative distances between smORFs based on their pathway enrichment. Finally, clusters are displayed in UMAP. For a ",
           "given smORF, you can verify which genes are clustering together with it and which pathways were driving those genes ",
           "to cluster together."),style = "text-align: justify;")
}


















