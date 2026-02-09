## Configuration for genetic information analysis workflow

genetic_information_config <- list(
  # Relative paths within the project (from project root)
  paths = list(
    annotation_default_csv = "Example/Example_smORFs.csv",
    exon_file_output      = "data_preparation/Genetic_information/Example_exon_file.txt",
    intermediate_gtf       = "data_preparation/Genetic_information/Example_smORFs_format.gtf",
    gtf_output             = "data_preparation/Genetic_information/Example_smORFs.gtf",
    perl_script_dir        = "data_preparation/Genetic_information",
    python_script          = "data_preparation/Genetic_information/bio_parse_orf_gtf_Hubner_no_length_restriction_maf.py",
    species_names_repo     = "data_preparation/Genetic_information/species_names.txt"
  ),
  
  # Assets folder structure (relative to project root)
  assets = list(
    base_dir        = "assets",
    maf_files       = "assets/maf_files",
    gwas_catalog    = "assets/gwas_catalog",
    species_names   = "assets/species_names",
    clustalo        = "assets/clustalo"
  ),
  
  # External tool URLs
  urls = list(
    maf_download    = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/",
    gwas_catalog    = "https://www.ebi.ac.uk/gwas/docs/file-downloads",
    clustalo        = "http://www.clustal.org/omega/"
  ),
  
  # File patterns to look for
  file_patterns = list(
    maf_files       = "\\.maf$",
    gwas_catalog    = "gwas_catalog.*\\.tsv$",
    species_names   = "species_names\\.txt$"
  ),
  
  # Tool detection
  tools = list(
    clustalo_cmd    = "clustalo"
  ),
  
  # Default settings
  defaults = list(
    source_name     = "Ho_et_al_2020",
    expected_gwas_size_mb = 300
  )
)
