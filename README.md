smORF-pro
================

<!-- # smORF-pro -->

## Structure

This readme is for running all the prerequisites data analysis to be
done before being able to start the shiny app. Once you have completed
this readme, you can integrate all the results into one shinyapp by
following the steps in the README_ShinyApp.Rmd .

  
  
There are 4 categories in the shiny app with all their different
analysis:  
1. Amino Acid analysis: Deeploc, DeepTMHMM, Target P, Interproscan  
2. Genetic Information: Conservation, GWAS  
3. Expression Specificity: GTex, RNA and RIBO, scRNAseq  
4. Ontology and clustering: Correlation, GSEA, GSEA clustering  
  
Each of these categories uses different tools to run analyses and
require a specific format. This page describes the data needed and
format required for each. Categories 1 2 and 3 are independent so it is
possible to integrate only one tab for new smORFs if you don’t have the
information required for that tab. Category 4 depends on category 3 in
order to run.  
  
Once all the analysis from point I to point IV are done, if you want to
build your own ShinyApp, you’ll need to create one rds object per smORF
and pre-generate plots for every smORF so that the ShinyApp can run
smoothly. Finally, you can head to the 2nd readme to visualize your own
ShinyApp.

## I - Amino Acid Analysis

In this part, we describe the tools that are being used which rely
exclusively on the amino acid sequences of the smORFs. Most of these
tools need a fasta file as input so the first step is to create a fasta
file from the smORFs name and amino acid sequence. The results will be
added as columns to the data.frame “Annotation”.

``` r
# Load example files
df=read.csv("data/input_example.csv")

# Annotation format
Annotation=data.frame(iORFID=df$iORF_ID,
                      Gene_id=df$Gene_id,
                      Gene_name=df$Gene_name,
                      ORF_type=df$ORF_type,
                      Gene_type=df$Gene_type,
                      Source=NA,
                      peptide_seq=df$Peptide_Seq,
                      Start=df$Start_codon)```
                      
# Make fasta
library(seqinr)
library(stringr)

write.fasta(as.list(Annotation$peptide_seq),names=Annotation$iORFID,
            as.string = F,"data_preparation/Amino_Acid_Analysis/Example_smORFs.fa")
```

### DeepTMHMM

To run DeepTMHMM, go to <https://dtu.biolib.com/DeepTMHMM> , load your
fasta file and run the analysis. Download results in 3line format and
save it in this folder data_preparation/Amino_Acid_Analysis/ . Then
process the results to integrate the results in your Annotation
data.frame:

``` r
File_3line=read.delim("data_preparation/Amino_Acid_Analysis/predicted_topologies.3line",
                      header = F)
File_3line=File_3line[substr(File_3line$V1,0,1)==">",]
File_3line=as.data.frame(sub(">","",File_3line))
File_3line[c("iORFID","TMHMM")]=str_split_fixed(File_3line[,1]," \\| ",2)

Annotation$deepTMHMM=File_3line$TMHMM[match(Annotation$iORFID,File_3line$iORFID)]
```

### TargetP 2.0

To run Target P 2.0 , go to
<https://services.healthtech.dtu.dk/services/TargetP-2.0/> ,load your
fasta file. Select “Non-plant” and “Short output (no figures)”. Submit
your job.Once completed, in the downloads list, select “Prediction
summary” and save your file to this folder
data_preparation/Amino_Acid_Analysis/ . Then process the results to
integrate the results in your Annotation data.frame:

``` r
targetP=fread("/home/baptiste/smORF-pro/data_preparation/Amino_Acid_Analysis/output_protein_type.txt")
Annotation$TargetP=targetP$Prediction[match(Annotation$iORFID,targetP$`# ID`)]
```

### Deeploc 2.0

Go to <https://services.healthtech.dtu.dk/services/DeepLoc-2.0/> . Load
your fasta file and select “High-quality (Slow)” and “Short output (no
figures)”. Submit your job. Download the “CSV summary” file to your
folder data_preparation/Amino_Acid_Analysis/ . If you have too many
proteins to run deeploc through the web interface, follow the
instructions to download it to run it locally, you can install it
through Anaconda and run it with the command:  
deeploc2 -f fasta_file.fa -m Accurate -p -o Output_folder  
It will create one figure per protein in that folder and 1 summary
table.  
Move the summary table to your folder
data_preparation/Amino_Acid_Analysis/

``` r
deeploc=read.csv("data_preparation/Amino_Acid_Analysis/deeploc_results.csv") # change name of file accordingly
Annotation$Deeploc=deeploc$Localizations[match(Annotation$iORFID,deeploc$Protein_ID)]
```

### Interproscan

Go to <https://www.ebi.ac.uk/interpro/search/sequence/> . Load your
fasta file and click on search to submit your job. Once completed you
can click on “Group Actions” and select download TSV output and save it
in the folder data_preparation/Amino_Acid_Analysis/ .  
If you have too many proteins (\>100 proteins) to run interproscan
through their web interface, follow the instructions to download it to
run it locally. <https://interproscan-docs.readthedocs.io/en/latest/> .
You can install it by following the website’s instructions, then, run
this command line in the folder where you installed it:  
./interproscan.sh -i new_smORFs.fa -f tsv -o
New_smORFs_interproscan.tsv  
It will create a tsv file with the summary from interproscan. Move this
file to data_preparation/Amino_Acid_Analysis/

``` r
library(data.table)
interpro=fread("data_preparation/Amino_Acid_Analysis/interpro_example_output.tsv") # change name accordingly
colnames(interpro)=c("smORF_ID","md5","length","Analysis","Signature_accession",
                     "Signature_description","start","stop","score","Status","Date",
                     "interpro_accession","interpro_description", "GO_Annotation", "Pathway_Annotation")
smORFs_with_signals=unique(interpro$smORF_ID)

Annotation$interproscan="No domain found"
for (smORF_i in smORFs_with_signals) {
  interpro_i=interpro[which(interpro$smORF_ID==smORF_i),]
  signals_i=unique(c(unlist(interpro_i[,c(6,13)])))
  signals_i=signals_i[-which(signals_i=="-")]
  Annotation$interproscan[which(Annotation$iORFID==smORF_i)]=paste(signals_i,collapse = "; ")
}

saveRDS(Annotation,"data/Annotation_Example.rds")


library(DT)
df_DT=Annotation
df_DT[,c(4,5,8,9,10,12)]=as.factor(df_DT[,c(4,5,8,9,10,12)])

saveRDS(df_DT,"data/Annotation_DT.rds")
```

## II - Genetic information

In this part, we describe how the methods used for analysis based on
genetic information such as conservation and finding existing mutations
from GWAS catalog. For both analysis, a gtf file of the smORFs is
required. If you don’t have a gtf format yet but you have all the CDS
coordinates, you can process it to turn it into a gtf.

``` r

Example_genetic_info=read.csv("Example_smORFs.csv")

# Create 1 line per exon
Example_genetic_info=data.frame(rbindlist(list(Example_genetic_info[,1:19],
                                     Example_genetic_info[,c(1:17,20,21)],
                                     Example_genetic_info[,c(1:17,22,23)],
                                     Example_genetic_info[,c(1:17,24,25)])))
Example_genetic_info$iORF_id=Example_genetic_info$ORF_id

# Remove lines of smORFs with no exon information and reorder the annotation
Example_genetic_info=Example_genetic_info[which(!is.na(Example_genetic_info$S_exon1)),]
Example_genetic_info=Example_genetic_info[order(Example_genetic_info$iORF_id),]

# make dataframe in format for the perl ExonToGTF.pl script

Exon_file=data.frame(chr=Example_genetic_info$Chr, # 0 line number for later perl script
                     start=Example_genetic_info$S_exon1, # 1
                     end=Example_genetic_info$E_exon1, # 2
                     strand=Example_genetic_info$strand, # 3
                     gene=Example_genetic_info$gene_id, # 4
                     ORF=Example_genetic_info$iORF_id, # 5
                     exon_id=paste0(Example_genetic_info$Chr,Example_genetic_info$S_exon1,Example_genetic_info$E_exon1), # 6
                     gene_name=Example_genetic_info$gene_name, # 7
                     orf_type=Example_genetic_info$iORF_type, # 8
                     gene_type=Example_genetic_info$gene_biotype, # 9
                     pept=Example_genetic_info$Peptide.seq, # 10
                     Source="Ho et al. 2020") # 11

write.table(Exon_file,"Example_exon_file.txt",row.names = F,col.names = F,
            quote = F,sep = "\t")
```

Once the CDS file is ready, you can convert it into a gtf by running the
2 perl scripts in the perl_script directory:

    perl exonToGTF.pl Example_smORFs_format.gtf Example_exon_file.txt
    perl collapseGTF Example_smORFs.gtf Example_smORFs_format.gtf 

Now this GTF is usable for calculating the conservation and finding if
the CDS parts of the smORFs are overlaping with GWAS hits.  
NOTE: This GTF only contains CDS regions and can only be used for
riboseq data analysis.

### Conservation

Make sure your gtf is in the same format as the example gtf given on
github, or that it has been created with the previous section. This
portion of the code is going to use the gtf coordinates of the smORFs
CDS to align them with the 100 vertebrates maf to reveal which proteins
are conserved in which species.  
  
Need to detail how to obtain the maf files and how to process them.  
  
Create a folder fasta, where the nucleotides alignments will be stored,
and create a folder protein within the fasta folder, this is where the
protein alignment will be stored. Modify the python script
bio_parse_orf_gtf_no_length_restriction_maf.py to specify your input
files gtf and maf accordingly as well as your paths to the output
directories fasta and protein. Run the python script:

    python bio_parse_orf_gtf_no_length_restriction.py

If quotes are present in the names of your smORFs in the output folders
fasta and proteins you can remove them in R

``` r
# nucleotide fasta files are not gonna be used further.
setwd("Conservation/fasta/")
to_rename=list.files(pattern = "*.fa")
new_names=gsub('"',"",to_rename)
file.rename(to_rename,new_names)

# Protein fasta are the ones used for plots
setwd("protein/")
to_rename=list.files(pattern = "*.txt")
new_names=gsub('"',"",to_rename)
file.rename(to_rename,new_names)
```

If you prefer to use the common species names instead of the genome
name, you can convert the names

``` r
# Make sure your setwd is in the protein folder before running this code
for (i in list.files(pattern = ".txt")) {
  file_i=fread(i,header = F)
  idx_to_rm=which(file_i[,1]=="")
  if (length(idx_to_rm)>0) {
    idx_to_rm=c(idx_to_rm,idx_to_rm-1)
    file_i=file_i[-idx_to_rm,]
  }
  order_to_move=c()
  for (sp in 1:100) {
    if (length(grep(species$GENOME[sp],unlist(file_i)))>0) {
      order_to_move=c(order_to_move,grep(species$GENOME[sp],unlist(file_i)))
      order_to_move=c(order_to_move,order_to_move[length(order_to_move)]+1)
    }
    file_i=as.data.frame(apply(file_i, 1, str_replace, pattern=species$GENOME[sp],
      replacement=species$COMMON[sp]))
    reordered=as.data.frame(file_i[order_to_move,])
  }
  write.table(reordered,i,row.names = F,col.names = F,quote = F)
}
```

You can then create the alignments using clustal omega. Download
clustalo from <http://www.clustal.org/omega/> . Create a msa and a
dist_mat folders in the protein folder. Run the run_clustalo.sh script.

### Overlaping smORFs CDS and known GWAS mutations

In this section we describe how to use the smORF gtf to find mutations
in their CDS region from GWAS catalog. You can download GWAS catalog
from XXXX . The goal is to find if there are any overlap, and in the
case of an overlap, if it results to a misense or non-sense mutation.

``` r
library(rtracklayer)
library(GenomicRanges)


GWAS=fread("data_preparation/Genetic_Information/GWAS/gwas_catalog_v1.0-associations_e111_r2024-01-19.tsv")
GWAS$CHR_POS=as.numeric(GWAS$CHR_POS)
GWAS=GWAS[!is.na(GWAS$CHR_POS),]

GTF=import("Example_smORFs.gtf")
GTF=GTF[GTF$type=="orfCDS",]


# Start the overlap
smORF_granges=GRanges(IRanges(start=GTF$start,end=GTF$end),seqnames=GTF$seqnames,iORF_ID=GTF$iORF_id)
GWAS_granges=GRanges(IRanges(start=GWAS$CHR_POS,end=GWAS$CHR_POS),seqnames=GWAS$CHR_ID,SNP_ID=GWAS$SNPS)

smORF_granges_list  = GNCList(smORF_granges)
overlap = findOverlaps(GWAS_granges, smORF_granges_list)

GWAS_results=cbind(GTF[overlap@to,],GWAS[overlap@from,])
write.csv(GWAS_results,"data_preparation/Genetic_Information/GWAS/GWAS_smORF_overlap.csv")


GWAS_results=read.csv("~/smORF-pro/data_preparation/Genetic_Information/GWAS/GWAS_smORF_overlap.csv")[,-1]
GWAS_results$peptide_seq=NA
GWAS_results$seq_after_mut=NA
GWAS_results$seq_change=NA

for (i in 1:nrow(GWAS_results)) {
  print(i)
  smORF_object=readRDS(paste0("~/smORF-pro/smORF_objects/",GWAS_results$iORF_id[i],".rds"))
  nt_seq=smORF_object$nt_seq
  nt_seq=gsub("-","",nt_seq)
  nt_seq=toupper(nt_seq)
  original_AA_seq=c2s(translate(s2c(nt_seq)))
  GWAS_results$peptide_seq[i]=original_AA_seq
  smORF_i_gtf=smORF_object$gtf
  
  GWAS_i=GWAS_results[i,]
  change="No"
  snp=str_sub(GWAS_i$STRONGEST.SNP.RISK.ALLELE[],
              start=nchar(GWAS_i$STRONGEST.SNP.RISK.ALLELE),
              end=nchar(GWAS_i$STRONGEST.SNP.RISK.ALLELE))
  if (snp%in%c("A","C","G","T")) {
    pos_snp=GWAS_i$CHR_POS
    print(nrow(smORF_i_gtf)==1)
    if (nrow(smORF_i_gtf)==1) {
      if (smORF_i_gtf$strand=="+") {
        pos_relative_snp=pos_snp-smORF_i_gtf$start+1
        substr(nt_seq,pos_relative_snp,pos_relative_snp)=snp
      }
      if (smORF_i_gtf$strand=="-") {
        pos_relative_snp=smORF_i_gtf$end-pos_snp+1
        substr(nt_seq,pos_relative_snp,pos_relative_snp)=chartr("ATGC","TACG",snp)
      }
    }
    new_seq=c2s(translate(s2c(nt_seq)))
    GWAS_results$seq_after_mut[i]=new_seq
    if (new_seq!=original_AA_seq) {
      change="Yes"
    }
  }
  GWAS_results$seq_change[i]=change
}
```

## III - Expression specificity

By using available public or your own RNA-seq and ribo-seq data, you can
identify which cell lines / tissues your smORFs transcripts are being
expressed and translated.

### GTex

For a public RNA-seq database, we are using GTex, which has over 10k
samples of different tissues. Download the TPM count file as well as the
annotation file from GTex website XXXX. There are no processing to be
done at this steps, the plots will directly be generated in the
“Pre-generate png” section of this readme.

### RNA specificity

If you have your own RNA-seq samples or public dataset that you want to
use for your smORFs, recount the dataset using a regular nsembl gtf or
your prefered version and your smORF GTF if you have smORFs in
unnanotated regions.

All you have to do is to make TPM from your count table and create a
sample annotation file.

``` r
# Example of RNA annotation file
coldata_rna=read.table("coldata_rna.txt",header = F)
rnaseq=fread("count_rna.txt)

# Make TPM for RNA
L=rnaseq$Length
x=rnaseq[,7:ncol(rnaseq)]/L
TPM=t(t(x)*1e6/colSums(x))
TPM=as.data.frame(TPM)
rownames(TPM)=make.names(rnaseq$Geneid,unique = T)

# colnames can't start with a number and - to replace with .
colnames(TPM)=str_replace_all(colnames(TPM),"-",".")

write.csv(TPM,"RNA_TPM.csv",row.names = T)
```

### Ribo specificity

Similar to RNA-seq, you can count your ribo-seq or public ribo-seq data
using your gtf. The gtf made earlier can be used for that if you count
using the orfCDS feature. Then just make your TPM and annotation file
from your count table.

``` r
# Example of ribo annotation file
coldata_ribo=read.table("coldata_ribo.txt",header = F)
riboseq=fread("count_ribo.txt)

# Make TPM for ribo
L=riboseq$Length
x=riboseq[,7:ncol(riboseq)]/L
TPM=t(t(x)*1e6/colSums(x))
TPM=as.data.frame(TPM)
rownames(TPM)=make.names(riboseq$Geneid,unique = T)

# colnames can't start with a number and - to replace with .
colnames(TPM)=str_replace_all(colnames(TPM),"-",".")

write.csv(TPM,"ribo_TPM.csv",row.names = T)
```

## IV - Correlation analysis and Ontology

To functionalize the smORFs based on their expression pattern, we rely
on coexpression analysis and gene ontology. These parts are quite
computationally heavy and may require the use of high-performance
computer (HPC) depending on your computer specs.

### Gene pairwise correlation analysis

The first step is to run gene-gene and protein-protein pairwise
genome-wide Spearman correlation. This can be done over all tissues /
cell lines or per tissue / cell line if you have enough samples per
tissue. We recommend at least 20 samples per tissue / cell type. It can
be done the same way for GTex, RNA, and ribo data.

``` r
library(data.table)
library(WGCNA)
library(gridExtra)
library(stringr)
library(tidyverse)
library(dplyr)

# Load data
Annotation_smORFs=read.csv("~/smORF-pro/V2/data/Annotation_full.csv")
TPM=read.csv("Expression/Counts/Atlas_Ribo_All_TPM.csv",row.names = 1)

# Select tissues of interest


load("~/Desktop/PhD/smorfs functionalization/data/tpm_data.Rdata")
coldata_ribo=tpm_data[[1]]
table(coldata_ribo$V2)
#tissues=unique(coldata_ribo$V2)

tissues=c("Fibroblasts","Heart DCM","All")
counts=TPM

for (t in tissues) {
  print(t)
  if (t!="All") {
    counts_i=counts[,colnames(counts)%in%coldata_ribo$V1[coldata_ribo$V2==t]]
    counts_i=counts_i[rowSums(counts_i>0)>=10,]

    # Run correlation
    coexpr=corAndPvalue(t(counts_i), nThreads = 20) # change your nThreads number depending on your computer specs
    corr=coexpr$cor
    pval=coexpr$p

    # Set no significant correlations to a correlation of 0
    #corr[which(pval>0.05,arr.ind = T)]=0
    saveRDS(corr,file=paste0("Expression/Coexpression/Ribo/Coexpression_tables/coexpr_smORFs_Ribo_",t,".RData"))
    saveRDS(pval,file=paste0("Expression/Coexpression/Ribo/Coexpression_tables/coexpr_smORFs_pval_Ribo_",t,".RData"))
  }
  if (t=="All") {
    counts_i=counts_i[rowSums(counts_i>0)>=10,]
    # Run correlation
    coexpr=corAndPvalue(t(counts_i), nThreads = 20)
    corr=coexpr$cor
    pval=coexpr$p

    # Set no significant correlations to a correlation of 0
    # corr[which(pval>0.05,arr.ind = T)]=0
    saveRDS(corr,file=paste0("Expression/Coexpression/Ribo/Coexpression_tables/coexpr_smORFs_Ribo_",t,".RData"))
    saveRDS(pval,file=paste0("Expression/Coexpression/Ribo/Coexpression_tables/coexpr_smORFs_pval_Ribo_",t,".RData"))
  }
}
```

### Ontology GTex, RNA, Ribo

For ontology

### Ribo ontology clustering
