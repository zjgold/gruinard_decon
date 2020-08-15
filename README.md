# Gruinard decon
An eDNA metabarcoding decontamination pipeline

This R script conducts a 6 step decontamination protocol on output Anacapa community tables. The objective of this script is to fit in between the Anacapa Toolkit and ranacapa, providing a series of cleaning and pre-processing steps to remove contaminant ASVs and poorly sequenced samples prior to data analysis and exploration.

**NOTE:** This script is currently underdevelopment and in version 0.0 . However, I decided to share this code incase any of the pieces maybe of interest to the larger eDNA metabarcoding community. The code will be updated during the summer and fall of 2020 to improve functionality.

I want to acknowledge that this script is built off of the work of many other great metabarcoding scientists including heavily relying on code from Ramon Gallego, Ryan Kelly, and the microDecon package.

#Cleaning Process
**Cleaning Process 0: Remove all Forward, Reverse, and Unmerged reads & Remove Singletons	Used to subset ASVs by read type.** Currently options are a) merged or b) all reads. Will add merged & forward selection
Cleaning Process 1: Estimation of *Tag-jumping* or sample *Cross-talk*	 With this procedure, we substract the proportion of reads that jumped into control samples from each environmental sample. First we determine which ASVs came from controls vs environmental samples. Next we create a vector of ASVs in positive controls. We then calculate what proportion of the reads found in the positive controls are found in environment samples.  Next, we substract the composition of the positive controls from the environment samples. The idea behind this procedure is that we now know how many reads from each ASV appeared in the controls. These come from 2 processes: sequences we know should appear in the positive controls, and sequences that have *hopped* from the environment to the positive controls.

**Cleaning Process 2: Discarding PCR replicates with low number of reads**	We fit the sample read depth to a beta distribution. We then discard samples below the 5% probability level.

**Cleaning Process 3: microDecon - Clearance of Negative Control Contaminants**	We use the microDecon package to identify contaminant ASVs that occurred in field, extraction, and PCR negative controls. microDecon compares the prevalence of ASVs in blanks to the prevalence in samples and removes contaminant ASVs and subtracts contaminant sequences.

**Cleaning Process 4: Site Occupancy Modeling**	Site occupancy modeling is used to asses whether an ASV is a true rare sequence or likely contaminant. We first identify the pattern of presence of each ASV across biological and technical replicates. We then estimate the occupancy probability of each ASV pattern across replicates based on the pattern of presesence across all samples. We then filter out ASVs with low occupancy probabilities.

**Cleaning Process 5: Dissimilarity between PCR (biological) replicates**	This step removes samples for which the dissimilarity between PCR replicates exceeds the normal distribution of dissimilarities observed in samples. The objective of this step is to remove any technical replicates that look like they do not belong.


#Background on Gruinard Decon
We decided to name this decontamination script after the infamous decontamination of Gruinard Island in Scotland. This island has a long and dark history of contamination arising from anthrax biological warfare testing during WW2 in which local peoples were displaced from the island and hundreds of domestic animals were sacrificed. Gruinard Island was abandoned as left uninhabitable for decades, but was eventually decontaminated through the dumping of formaldehyde across the entire island to remove anthrax spores. This was both one of the largest biological warfare testing as well as one of the largest decontamination efforts undertaken. We hope this pipeline will provide a useful tool for removing pesky contaminant ASVs and leave your data slightly less contaminated than Gruinard Island.
