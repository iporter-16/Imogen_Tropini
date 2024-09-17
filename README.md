# Imogen @ Tropini lab

This contains all (productive and/or successful) scripts written during my time at the lab, alongside the data files required to run said scripts.

### CC BP Figure 6F

The script **CC_Bowel_Prep_6F.Rmd** is used to produce the two components of Fig. 6F in Clayton et al. (*Bowel preparation-induced changes in the gut environment facilitate pathogen and pathobiont colonization*). 

**C3_basic_metadata.xlsx** contains the sample metadata (tissue type, timepoint, and BP). 

The following data were derived from the QIIME2 pipeline described in my lab notebook (*IP-CC1.1 C3 MiniSeq 16S analysis & translocation focus*):

**C3_L6_data_w_mapping.csv** contains the QIIME2-produced relative abundance table at the 6th taxonomic level (genus), with metadata attached.

**C3_MiSeq_pivot_translocation.csv** contains the same data, filtered for taxa that were found translocated to the MLN by 24h post-BP.

### GH Analysis

Designed to read through glycoside hydrolase annotations (created by dbCAN3) on multiple genomes and output heatmaps summarising substrate and source specificities. 

**GH_annotations.Rmd** script run to produce said heatmaps. Requires family information (fam-substrate-mapping-08012023.tsv), and a directory containing our annotated genomes (see below)

**fam-substrate-mapping-08012023.tsv** obtained from dbCAN3, contains information about each GH family (substrate, enzyme name, etc.)

**Muribaculaceae_genomes/Annotations/** a directory containing the annotated genomes on which analysis is performed. Each sub-directory is named according to its genome, and contains the many outputs of dbCAN3 annotation. For this pipeline, specifically the **hmmer.out** file is required for each.