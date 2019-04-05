# C239SPeronosporaleanEvoution
This is a project that studies the evolution of a C239S mutation in Peronosporalean oomycetes that coincides with resistance to ethaboxam

This project is an extension of a project to develop a high-throughput fungicide sensitivity assay for oomycetes and and characterize sensitivity to mefenoxam and ethaboxam, and can be found here: [Soybean Associated Oomycete Fungicide Sensitivity](https://github.com/noelzach/Community_Fungicide_Sensitivity)

Analysis of oomycete fungicide sensitivity can be found in that repository

 
## Data associated with this project are either phylogeny files or data analysis files listed below 

* Figure 1 was generated based on data presented based on the ethaboxam phenotypes presented in [Soybean Associated Oomycete Fungicide Sensitivity](https://github.com/noelzach/Community_Fungicide_Sensitivity). [Supplemental table 1](SupplementalTable1/SupplementalTable1_isolatelist.xlsx) contains the ethaboxam phenotypes used in Figure 1. 
* Figure 2 was generated based on a protein alignment by cloning and sequencing beta tubulin from oomycete cultures (Accession Numbers MK752959-MK753004). [Supplemental Table 2](SupplementalTable2/SupplementalTable2_primers.xlsx) includes the primers and vectors used to clone oomycete beta tubulin. 
* Figure 3 phylogeny along with the ancestral sequence reconstruction information can be found in the Fig3 folder.The MCMC runs were too large to include in this repository. If you would like them please email me (noelzach@msu.edu). 
* Figure 4 was generated with the code contained within the Fig4 foler. [Treated Seed Virulence Assay](Fig4/seedrot.md) uses the seedinfectionassay_clean.csv file.[Supplemental Figure 1](SupplementalFigure1/SupplementalFigure1.pdf) was used as a guid to score disease severity and convert to DSI score to then do statistics. 
* Figure 5 was generated based on information in the Fig. 5 folder. Commands in PeronosporaleanPhylogeny.txt were used to make Figure 5 and expanded version in Supplemental figure 2. These commands use Peronosporalean_betatubulin.fasta to make the alignment with MAFFT. Then use FastTree to make the .newick file. Then upload the .newick file to iTOL. Then load in labels.txt and StarBanches.txt to fully reproduce Supplemental Figure 2. 

* Additionally all sequences associated with this project can be found under Genbank accession numbers MK752959-MK753004. Or in the [NCBI uploaded folder](NCBI_uploaded)



