# Functions called in the interface for 1.shiny.R

Plot_DE_smorfs_per_exp=function(smORFs_padj){
  ggplot(smORFs_padj,aes(x=Experiment,y=Nb_smORFs_DEG))+
  geom_bar(stat = "identity",fill="dodgerblue")+
  theme_bw()+ylab("Number of DEG containing smORFs")+ggtitle("DEG containing smORFs per experiment")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
}


Plot_DE_smorfs_many_exp=function(smORFs_DEG_per_exp){
  ggplot(smORFs_DEG_per_exp,aes(x=Var1,y=Freq))+
    geom_bar(stat = "identity",fill="dodgerblue")+
    theme_bw()+ylab("Number of smORFs DEG")+xlab("Number of experiments")+
    ggtitle("DEG containing smORFs in multiple experiments")
}


Plot_volcano=function(Experiment_name){
  df=DEG[[Experiment_name]]
  df=df[!is.na(df$padj),]
  ggplot(df,aes(x=df$log2FC,y=-log10(df$padj),color=rownames(df)%in%Annotation$Gene_id))+ geom_point()+
    theme_classic()+xlab("Log2FoldChange")+ylab("-log(padj)")+
    ggtitle(paste0("smORFs fold change in ",Experiment_name))+
    guides(color=guide_legend(title = "smORFs"))+
    scale_color_manual(values = c("black","red"))
}

Plot_heatmap=function(Threshold=9){
  DEG_df=bind_cols(DEG)
  DEG_LFC=DEG_df[,1:(ncol(DEG_df)/2)*2-1]
  colnames(DEG_LFC)=names(DEG)
  DEG_LFC=DEG_LFC[rownames(DEG_LFC)%in%Annotation$Gene_id,]
  DEG_LFC=DEG_LFC[which(rowSums(abs(DEG_LFC)>=Threshold)>0),]
  pheatmap(DEG_LFC,color= brewer.pal(9,"RdBu"),
           main = "Heatmap of Log2FoldChange of smORFs per experiment",
           legend_breaks = c(-10,-5,0,5,10,15),
           legend_labels = c("-10","-5","0","5","10","Log2FC"),
           border_color = NA,angle_col = 45)
}

Plot_ORF_domain=function(SEP=ORF){
  if (ORF%in%TargetP$iORF_ID) {
    i=TargetP$`# ID`[which(TargetP$iORF_ID==ORF)]
    Target_P_i=fread(paste0(folder_targetP,"output_",i,"_pred.txt"))
  }
}

Plot_conservation=function(SEP=ORF){
  if(which(Annotation$iORF_id==SEP)<7768){
  fasta_name=paste0("data/Conservation_fasta/",SEP,".fa")
  ggmsa(fasta_name, char_width = 0.4, seq_name = T,) + geom_seqlogo() + geom_msaBar()
  }
}

# Functions used for smORF card tab

smORF_text=function(ORF="ENCT00000011417_ncRNA_3730"){
  df_i=Annotation[which(Annotation$iORF_id==ORF),]
  txt=paste(ORF,"is a" ,nchar(df_i$Peptide.seq)-1,"amino acid long smORF. It is localized on the gene",df_i$Gene_name,df_i$Gene_id,
            "which is annotated as",str_replace(df_i$Gene.type,"_"," "),"gene")
  return(txt)
}

smORF_DE_plot=function(ORF="ENCT00000011417_ncRNA_3730"){ # Plot the LFC of the smORF of interest accross the Fantom datasets
  df=data.frame(row.names = names(DEG),Log2FC=rep(0,length(DEG)),padj=rep(1,length(DEG)))
  ORF_gene=Annotation$Gene_id[which(Annotation$iORF_id==ORF)]
  if (ORF_gene%in%rownames(DEG[[1]])) {
    for (i in 1:length(DEG)) {
      df$Log2FC[i]=DEG[[i]]$log2FC[which(rownames(DEG[[i]])==ORF_gene)]
      df$padj[i]=DEG[[i]]$padj[which(rownames(DEG[[i]])==ORF_gene)]
    }
    ggplot(df,aes(x=rownames(df),y=Log2FC,color=-log10(padj),size=-log10(padj)))+geom_point()+
      theme_classic()+geom_hline(yintercept = c(-1,0,1),linetype=c("dashed","dotted","dashed"))+
      scale_color_binned(breaks=c(0,1,2,5,10),type="viridis")+xlab("")+
      scale_size_binned(breaks=c(0,1,2,5,10))+
      theme(axis.text.x = element_text(angle = 45,hjust=1))
  }
}



smORF_tpm_plot=function(ORF="ENCT00000011417_ncRNA_3730"){ # plot the tpm RNA and RIBO for a given smORF
  ribo_tpm=tpm_data[[3]]
  rna_tpm=tpm_data[[4]]
  coldata_ribo=tpm_data[[1]]
  coldata_rna=tpm_data[[2]]
  ribo_i=ribo_tpm[rownames(ribo_tpm)==ORF,]
  rna_i=rna_tpm[rownames(rna_tpm)==ORF,]
  df_i=data.frame(Tissue=coldata_ribo$V2,RNA=0,RIBO=0)
  for (i in unique(df_i$Tissue)) {
    df_i_tissue_i=which(df_i$Tissue==i)

    sample_ribo_tissue_i=coldata_ribo$V1[which(coldata_ribo$V2==i)]
    df_i$RIBO[df_i_tissue_i]=as.numeric(ribo_i[sample_ribo_tissue_i])

    sample_rna_tissue_i=coldata_rna$V1[which(coldata_rna$V2==i)]
    df_i$RNA[df_i_tissue_i]=as.numeric(rna_i[sample_rna_tissue_i])
  }
  df_i$Tissue=as.factor(df_i$Tissue)


  # GTex
  IDs_to_keep=GTex_Annotation$SAMPID[GTex_Annotation$SAMPID%in%colnames(GTex)]
  geneID=Annotation$Gene_id[Annotation$iORF_id==ORF]
  df_gtex=data.frame(samp_ID=IDs_to_keep,Tissue=GTex_Annotation$SMTSD[GTex_Annotation$SAMPID%in%IDs_to_keep],RNA=0)

  df_gtex$RNA=as.numeric(GTex[geneID,df_gtex$samp_ID])
  df_i$Tissue=as.factor(df_i$Tissue)

  plot1=ggplot(df_i,aes(x=Tissue,y=RNA))+geom_boxplot(fill="deepskyblue1")+
    theme_bw()+xlab("")+ylab("TPM")+ggtitle("RNA expression")+
    theme(plot.title =element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45,hjust=1))
  plot2=ggplot(df_i,aes(x=Tissue,y=RIBO))+geom_boxplot(fill="darkorange2")+
    theme_bw()+xlab("")+ylab("TPM")+ggtitle("RIBO expression")+
    theme(plot.title =element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45,hjust=1))
  plot3=ggplot(df_gtex,aes(x=Tissue,y=RNA))+geom_boxplot(fill="deepskyblue1")+
    theme_bw()+xlab("")+ylab("TPM")+ggtitle("RNA expression")+
    theme(plot.title =element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45,hjust=1))
  grid.arrange(plot1,plot2,ncol=2,
               nrow=2,plot3,
               layout_matrix=rbind(c(1,2),c(3,3)))
}

min_padj=function(vec){
  order(vec)[1:n]
}

bubble_GO=function(Gene="ENSG00000204954"){
  n=5
  NES=fread(paste0("Coexpression/NES_per_gene/",Gene,".txt"))
  padj=fread(paste0("Coexpression/padj_per_gene/",Gene,".txt"))

  NES$V1=substr(NES$V1,14,nchar(NES$V1)-nchar(Gene)-7)
  colnames(NES)=c("Tissue",names(P_to_run))
  padj$V1=substr(padj$V1,16,nchar(padj$V1)-nchar(Gene)-7)
  colnames(padj)=c("Tissue",names(P_to_run))
  NES=as.data.frame(NES)
  padj=as.data.frame(padj)
  NES_melted=melt(NES)
  padj_melted=melt(padj)
  df=data.frame(Tissue=NES_melted$Tissue,Pathways=NES_melted$variable,NES=NES_melted$value,padj=padj_melted$value)
  path_selected=apply(padj[,-1],1,min_padj) # Select top n padj pathways
  for (i in 1:nrow(padj)) {
    path_selected[which(padj[i,path_selected[,i]+1]>0.05),i]=0
  }
  # unique pathways
  path_selected=unique(c(path_selected))
  if (0%in%path_selected) {
    path_selected=path_selected[-which(path_selected==0)]
  }
  path_selected=colnames(padj[,-1])[path_selected]

  toplot=df[df$Pathways%in%path_selected,]

  score_df=as.data.frame(-log10(padj[,-1]*NES[,-1]))
  score_df=score_df[,colnames(score_df)%in%path_selected]
  rownames(score_df)=NES$Tissue
  score_df[which(score_df=="NaN",arr.ind = T)]=0
  score_df[which(score_df=="Inf",arr.ind = T)]=sort(unique(unlist(c(score_df))),decreasing = T)[2]

  hclust_score <- hclust(dist(t(score_df)), method = "complete")
  ordre=hclust_score$order

  toplot$Pathways=factor(toplot$Pathways,levels=colnames(score_df)[ordre])

  ggplot(toplot,aes(x=Tissue,y=Pathways,size=-log10(padj),col=NES))+
    geom_point()+theme_classic()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    scale_color_gradient2()
}

