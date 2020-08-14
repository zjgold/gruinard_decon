
#' Function to summarize the ASVtable to see how many samples, hashes, and reads are preserved at each decontamination step
#' @param ASVtable Long-format ASV table generated from an Anacapa formated OTU table
#' @author Ramon Gallego
#' @export
#Function Name: How many function
#Use: To calculate number of hashes (ASVs), number of reads, and number of samples at ech stage of the decontamination pipeline
how.many <- function(ASVtable, round){
  ASVtable %>% ungroup() %>% 
    dplyr::summarise(nsamples = n_distinct(sample),
                     nHashes = n_distinct(seq_number),
                     nReads = sum(nReads), 
                     Stage = paste0("Step_", round)) %>% 
    gather(starts_with("n"), value = "number", key = "Stat")
}


#' Function to convert the cleaned.tibble generated in step 4 to a vegan matrix across all samples
#' @param tb Long-format table generated from step 4
#' @author Ramon Gallego
#' @export
#Function Name: tibble_to_vegdist_all
#Use: Caclulate vegan distance matrix across all samples
tibble_to_vegdist_all <- function (tb) {
  tb %>% 
    group_by(Miseq_run, sample, seq_number) %>% 
    dplyr::summarise(nReads = sum(Normalized.reads)) %>% 
    spread ( key = "seq_number", value = "nReads", fill = 0) %>% 
    ungroup() %>% unite(Miseq_run,sample, col="sample", sep=":")-> matrix_1
  samples <- pull (matrix_1, sample)
  matrix_1[,-1] -> matrix_1
  data.matrix(matrix_1) -> matrix_1
  dimnames(matrix_1)[[1]] <- samples
  vegdist(matrix_1) -> matrix_1
}


#' Function to convert the nested.cleaning generated to a vegan matrix for each biological replicate 
#' @param ASVtable Long-format ASV table generated from an Anacapa formated OTU table
#' @author Ramon Gallego
#' @export
#Function Name: How many function
#Use: Caclulate vegan distance matrix for each biological replicate
tibble_to_vegdist_bio_rep <- function (tb) {
  tb %>% 
    unite(col="location_id", Miseq_run,Site,Bio_rep,Tech_rep, remove = FALSE, sep=":") %>%
    group_by(location_id, sum.taxonomy) %>% 
    dplyr::summarise(nReads = sum(Normalized.reads)) %>% 
    spread ( key = "sum.taxonomy", value = "nReads", fill = 0) %>% 
    ungroup() -> matrix_1
  samples <- pull (matrix_1, location_id)
  matrix_1[,-1] -> matrix_1
  data.matrix(matrix_1) -> matrix_1
  dimnames(matrix_1)[[1]] <- samples
  vegdist(matrix_1) -> matrix_1
}

#' Function to convert the nested.cleaning generated to a vegan matrix for each biological replicate 
#' @param ASVtable Long-format ASV table generated from an Anacapa formated OTU table
#' @author Ramon Gallego
#' @export
#Function Name: How many function
#Use: Caclulate vegan distance matrix for each biological replicate
tibble_to_vegdist_bio_rep_single <- function (tb) {
  tb %>% 
    unite(col="location_id", Miseq_run,Site,Bio_rep, remove = FALSE, sep=":") %>%
    group_by(location_id, sum.taxonomy) %>% 
    dplyr::summarise(nReads = sum(Normalized.reads)) %>% 
    spread ( key = "sum.taxonomy", value = "nReads", fill = 0) %>% 
    ungroup() -> matrix_1
  samples <- pull (matrix_1, location_id)
  matrix_1[,-1] -> matrix_1
  data.matrix(matrix_1) -> matrix_1
  dimnames(matrix_1)[[1]] <- samples
  vegdist(matrix_1) -> matrix_1
}


#' Function for calculating distance to centroid from nested.cleaning
#' @param ASVtable Long-format ASV table generated from an Anacapa formated OTU table
#' @author Ramon Gallego
#' @examples
#' good_taxon_table <- data.frame(sum.taxonomy = c("a;b;c;d;f;u", "p;q;r;s;t;u"),
#' site_1 = c(0,1), site_2 = c(10, 20))
#' group_anacapa_by_taxonomy(good_taxon_table)
#' @export
#Function Name: How many function
#Use: Caclulate distance to the centroid for each biological replicate
dist_to_centroid <- function (x,y) {
  biol <- rep(y, dim(x)[[1]])
  
  if (length(biol) == 1) {
    output = rep(x[1]/2,2)
    names(output) <- attr(x, "Labels")
  }else{ 
    
    dispersion <- betadisper(x, group = biol)
    output = dispersion$distances
  }
  output
}