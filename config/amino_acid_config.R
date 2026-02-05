## Configuration for amino acid analysis workflow

amino_acid_config <- list(
  # Relative paths within the project
  paths = list(
    annotation_default_csv = "Example/Example_smORFs.csv",
    fasta_output           = "data_preparation/Amino_Acid_Analysis/Example_smORFs.fa",
    deeptmhmm_example      = "data_preparation/Amino_Acid_Analysis/predicted_topologies.3line",
    targetp_example        = "data_preparation/Amino_Acid_Analysis/output_protein_type.txt",
    deeploc_example        = "data_preparation/Amino_Acid_Analysis/Example_Deeploc_output.csv",
    interpro_example       = "data_preparation/Amino_Acid_Analysis/interpro_example_output.tsv",
    annotation_rds         = "data/Annotation.rds",
    annotation_dt_rds      = "data/Annotation_DT.rds"
  ),

  # External tool URLs
  urls = list(
    deeptmhmm = "https://dtu.biolib.com/DeepTMHMM",
    targetp   = "https://services.healthtech.dtu.dk/services/TargetP-2.0/",
    deeploc   = "https://services.healthtech.dtu.dk/services/DeepLoc-2.0/",
    interpro  = "https://www.ebi.ac.uk/interpro/search/sequence/"
  )
)

