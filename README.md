smORF-pro
================

<!-- # smORF-pro -->

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

## Amino Acid Analysis

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
your job. Once completed, in the downloads list, select “Prediction
summary” and save your file to this folder
data_preparation/Amino_Acid_Analysis/ . Then process the results to
integrate the results in your Annotation data.frame:

``` r
targetP=fread("/home/baptiste/smORF-pro/data_preparation/Amino_Acid_Analysis/output_protein_type.txt")
Annotation$TargetP=targetP$Prediction[match(Annotation$iORFID,targetP$`# ID`)]
```
