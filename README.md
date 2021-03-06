# Gruinard decon
An eDNA metabarcoding decontamination pipeline.  

This R script conducts a 6 step decontamination protocol on output [Anacapa Toolkit](https://github.com/limey-bean/Anacapa) community tables. The objective of this script is to fit in between the Anacapa Toolkit and [ranacapa](https://github.com/gauravsk/ranacapa), providing a series of cleaning and pre-processing steps to remove contaminant ASVs and poorly sequenced samples prior to data analysis and exploration.  

**NOTE:** This script is currently underdevelopment and in version 0.0 . However, I decided to share this code incase any of the pieces maybe of interest to the larger eDNA metabarcoding community even at this early stage. The code will be updated during the summer and fall of 2020 to improve functionality. Also please feel free to contact me and recommend any suggestions, utilities, functions, etc.

I want to acknowledge that this script is built off of the work of many other great metabarcoding scientists. I heavily relying on code from [Ramon Gallego](https://github.com/ramongallego?tab=repositories), [Ryan Kelly](https://github.com/invertdna), and the [microDecon package](https://github.com/donaldtmcknight/microDecon). Thank you for your dedication to open access software to provide coding resources to all.

## Cleaning Process
### **Cleaning Process 0: Remove all Forward, Reverse, and Unmerged reads & Remove Singletons	Used to subset ASVs by read type.**
Currently options are a) merged or b) all reads. Will add merged & forward selection.  
### **Cleaning Process 1: Estimation of Tag-jumping or Index Hopping**	 
With this procedure, we substract the proportion of reads that jumped into control samples from each environmental sample. First we determine which ASVs came from controls vs environmental samples. Next we create a vector of ASVs in positive controls. We then calculate what proportion of the reads found in the positive controls are found in environment samples.  Next, we substract the composition of the positive controls from the environment samples. The idea behind this procedure is that we now know how many reads from each ASV appeared in the controls. These come from 2 processes: sequences we know should appear in the positive controls, and sequences that have *hopped* from the environment to the positive controls.  

### **Cleaning Process 2: Discarding PCR replicates with low number of reads**
We fit the sample read depth to a beta distribution. We then discard samples below the 5% probability level.  

### **Cleaning Process 3: microDecon - Clearance of Negative Control Contaminants**
We use the microDecon package to identify contaminant ASVs that occurred in field, extraction, and PCR negative controls. microDecon compares the prevalence of ASVs in blanks to the prevalence in samples and removes contaminant ASVs and subtracts contaminant sequences.  

### **Cleaning Process 4: Site Occupancy Modeling**
Site occupancy modeling is used to asses whether an ASV is a true rare sequence or likely contaminant. We first identify the pattern of presence of each ASV across biological and technical replicates. We then estimate the occupancy probability of each ASV pattern across replicates based on the pattern of presesence across all samples. We then filter out ASVs with low occupancy probabilities.  

### **Cleaning Process 5: Dissimilarity between PCR (biological) replicates**
This step removes samples for which the dissimilarity between PCR replicates exceeds the normal distribution of dissimilarities observed in samples. The objective of this step is to remove any technical replicates that look like they do not belong.  

## Running Gruinard decon

### Example Command Line Script

```
Rscript /path/generalized_decontam_script.R -w /path/example_run -a path/example_anacapa_table.txt -y one -m path/example_metadata.txt -r merged_only -l 0.05 -u 10000 -n 2 -t 0.7 -p 5e-05 -g 0 -z double -c 10 -i 1000 -k 0.975 -s som_min_threshold_1
```
### Parameters
  **-w** --working directory , is the path to working directory (required)  
  **-a** --input_anacapa_path_1 , path to anacapa table 1 (required). This species (row) -site (column) community table needs to be in the Anacapa Toolkit output format. This requires the first column *seq_number* to contain the names of each unique ASV. The middle columns are the samples with read counts of each ASV within each sample. The final column *sum.taxonomy* contains the assigned taxonomy for each ASV.  
  **-b** --input_anacapa_path_2 , path to anacapa table 2  
  **-y** --number_anacapa_tables, Currently only able to run with *one* or  *two* anacapa tables. (required) Default is *two* tables. Other option is *one* table.  
  **-m** --input_metadata_path , path to metadata (required). The metadata table needs to be implemented in the exact same format as the example metadata table. The first column *Seq_number* should contain the names exactly matching the final output from the Anacapa Toolkit. *New_name* should contain updated names for each sample. Typically Seq_number names contain extraneous information and are not in an elegant format. This allows for the user to update the names generated by the Illumina sequencing machine. *Site* is the highest hierarchical level of replication with multiple *Bio_rep* or bottle replicates taken within a site and with multiple *Tech_rep* or PCR replicates taken from each bottle replicate. This is a typical hierarchical replication stratification employed in eDNA studies. However, alternative hierarchies can be employed for the site occupancy modeling. For example, *Bio_rep* could be replicate days in which a site was sampled with *Tech_rep* representing multiple bottles collected on the same day. However, the names of the columns must match these exactly. Finally the two final columns are *Sample_Control* which identifies which samples are a sample or a control as well as *Control_Type* which identifies whether a control is a *Pos* positive control or a *Blank* negative control.  
  **-r** --read_type , options are *merged_only*, *merged_and_forward*, or *all*. All includes all read types (forward,reverse,unmerged, and merged).  
  **-l** --low_dist_probability_cutoff , Threshold applied to the beta distribution of sample sequencing depth used to eliminate samples with low sequencing depth. Default is *0.05*. Set to 0.00000001 to ignore.  
  **-u** --minimum_read_cutoff , Minimum sample read depth threshold.   Default is 10000. Set to 1 to ignore.  
  **-n** --step3_runs , The run term used in step 3 for microDecon. See microDecon for detailed instructions. Default is 2.  
  **-t** --step3_thresh , The thresh term used in step 3 for microDecon. See microDecon for detailed instructions. Default is 0.7.  
  **-p** --step3_prop.thresh , The prop.thresh term used in step 3 for microDecon. See microDecon for detailed instructions. Default is 5e-5.  
  **-g** --step3_regression , The regression term used in step 3 for microDecon. See microDecon for detailed instructions. Default is 0.
  m**-z** --som_level , Options are double or single. Double hierarchy has 2 levels of replication (e.g. biological and technical replicates or biological replicates and multiple sampling events). A single hierarchy has one level of replication (e.g. only technical replicates, only biological replicates, or only multiple sampling events.  
  **-c** --som_chains , Number of chains used in site occupancy modeling. Default is 10.  
  **-i** --som_iterations , Number of iterations used in site occupancy modeling. Default is 1000.  
  **-k** --max_dist_probability_cutoff , Threshold applied to the normal distribution of pairwise replicate dissimilarity. This term is used to eliminate replicates that are highly dissimilar from eachother. Default is 0.975. Set to 0.00000001 to ignore.  
  **-s** --som_filter_choice , Site occupancy probability cutoff. Default is som_min_threshold_1 which is the max occupancy probability of a single technical replicate detection at a given site (e.g. removes species that occurred in only 1 technical replicate at any site). Other options are b) som_min_threshold_2 which is the max occupancy probability of a single technical replicate detection across all biological replicates at a given site (e.g. removes species that only occurred in only 1 technical replicate across all biological at any site even if it occurred once in multiple biological replicates), c) som_min_threshold_3 which is the minimum occupancy probability of 3 technical replicate detections at a given site (e.g. removes species that do not occur in at least 3 technical replicates at any site).  

## Background on Gruinard Decon
We decided to name this decontamination script after the infamous decontamination of [Gruinard Island](https://en.wikipedia.org/wiki/Gruinard_Island) in Scotland. This island has a long and dark history of contamination arising from anthrax biological warfare testing during WW2 in which local peoples were displaced from the island and hundreds of domestic animals were sacrificed. Gruinard Island was abandoned as left uninhabitable for decades, but was eventually decontaminated through the dumping of formaldehyde across the entire island to remove anthrax spores. This was both one of the largest biological warfare testing as well as one of the largest decontamination efforts undertaken. We hope this pipeline will provide a useful tool for removing pesky contaminant ASVs and leave your data slightly less contaminated than Gruinard Island.
