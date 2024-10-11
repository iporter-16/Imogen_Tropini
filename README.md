# Imogen @ Tropini lab

This contains all (productive and/or successful) scripts written during my time at the lab, alongside the data files required to run said scripts.

### CC BP Figure 6F

The script **CC_Bowel_Prep_6F.Rmd** is used to produce the two components of Fig. 6F in Clayton et al. (*Bowel preparation-induced changes in the gut environment facilitate pathogen and pathobiont colonization*). 

**C3_basic_metadata.xlsx** contains the sample metadata (tissue type, timepoint, and BP). 

The following data were derived from the QIIME2 pipeline described in my lab notebook (*IP-CC1.1 C3 MiniSeq 16S analysis & translocation focus*):

**C3_L6_data_w_mapping.csv** contains the QIIME2-produced relative abundance table at the 6th taxonomic level (genus), with metadata attached.

**C3_MiSeq_pivot_translocation.csv** contains the same data, filtered for taxa that were found translocated to the MLN by 24h post-BP.

### CC BP Figure 6F

For _Salmonella_ Tn library translocation and infection graphs. Data is from September 19 2024. Produces violin plots of GI infection, organ translocation. Takes in a .csv file that requires columns:

- Day : contains day at which mouse was sacc'd 
- Day_0, Day_1, Day_2 : CFU count estimates at sampling day (three columns)
- Liver, MLN, Spleen : CFU count estimate from each tissue

### 16S_MB_Clostridia

Analysis of AS's 16S data for mice gavaged with either a MB slurry, Clostridia slurry, or both slurries. 

**Startup.Rmd** takes in the components required to generate a phyloseq object (constructed originally in qiime2). The phyloseq object and taxonomy matrix are also contained as RData files for easy access.

- Alpha and beta diversity analyses performed
- Plotting relative abundance of specific taxa of interest, as defined in lines 170 onwards.

**DESEq2_analyses.Rmd** takes in the saved phyloseq object and taxa table from Startup.Rmd and performs differential abundance analyses using the DESeq2 package. Produces a graph of taxa significantly different in abundance between treatments.

Contains subdirectory **DP Data** which is 16S data from four treatment groups : conventional control (C), MB-naive (N), MB-positive (M), and MB-introduced (I). Data analysis takes in the phyloseq object and taxa table, converting it to the same dataframe format as in the main directory. This is then used to calculate log-ratio values and produce box plots, heatmaps, for each treatment.

### GH Analysis

Designed to read through glycoside hydrolase annotations (created by dbCAN3) on multiple genomes and output heatmaps summarising substrate and source specificities. 

**GH_annotations.Rmd** script run to produce said heatmaps. Requires family information (fam-substrate-mapping-08012023.tsv), and a directory containing our annotated genomes (see below)

**fam-substrate-mapping-08012023.tsv** obtained from dbCAN3, contains information about each GH family (substrate, enzyme name, etc.)

**Muribaculaceae_genomes/Annotations/** a directory containing the annotated genomes on which analysis is performed. Each sub-directory is named according to its genome, and contains the many outputs of dbCAN3 annotation. For this pipeline, specifically the **hmmer.out** file is required for each.

The data (substrate and source counts) produced by this code are also included for ease of access. The 'detailed' file adds the coverage and e values if this is of interest, as well as the more detailed GH profiles (for those with sub-profiles).

### AS_PEG_F0_Salmonella

Translocation and infection analysis performed on AS's F0 Salmonella data (mice treated with PEG in past or present generation). Produces bar + box plots (with error bars for SE) , as well as graph of weight change over time. 



