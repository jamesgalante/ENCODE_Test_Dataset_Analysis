This folder includes the raw data analysis / processing for the three perturb-seq datasets used in the held out dataset. All datasets were processed linearly through an R markdown file. Raw data can be found on GEO/ENCODE for each paper. The results of processing all dataset can be found in "results/" of each subdirectory (e.g. Xie/results/)


Morris/
- supplementary_tables/science.adh7699_table_s3.xlsx - supplementary table 3 from the Morris et al. 2023 publication. Used for both STINGv1 and STINGv2. Contains guide information
- Morris/STINGv1
  - blat_results.rds - the results from using UCSCs BLAT on all of the guides used in STINGv1. Note that this uses the BLAT API as there are so few guides.
  - data/ - the raw data downloaded from GEO. This is done automatically in the Rmd.
- Morris/STINGv2
  - out.psl - the results from using UCSCs BLAT on all of the guides used in STINGv2. Note that this output is from using BLAT on a local server. Instructions to run BLAT locally can be found in `using_BLAT.txt`
  - data/ - the raw data must be manually downloaded from GEO. Instructions are found in lines 92-94 in the STINGv2 Rmd file.
  
Xie/
- Global-analysis-K562-enhancers - this folder and it's subfolders are taken from the original Xie analysis github. Only a small part of this directory is used to create the guide UMI matrices. Detailed can be found in the Rmd.
- data/ - raw and filtered data from the original screen. Filtered data is generated in the Rmd file using python functions from Global-analysis-K562-enhancers. Instructions can be found in lines 33-79 of the Rmd
- out.psl - the results from using UCSCs BLAT on all of the guides used in STINGv2. Note that this output is from using BLAT on a local server. Instructions to run BLAT locally can be found in `using_BLAT.txt`

Klann/
- data/ - raw guide and gene count data downloaded from ENCODE. Instructions for downloading can be found in lines 27-55 of the Rmd.
- gr_object_final_df.csv - results from lifting the hg19 guide coordinates to hg38 using easylift. Instructions for doing so can be found in lines 44-115 in the `guide_info_creation.Rmd` file
- supplementary_table_16_scCERES_library_gRNAs.xlsx - supplementary table 16 from the original Klann et al. 2021 study used for information on the guides