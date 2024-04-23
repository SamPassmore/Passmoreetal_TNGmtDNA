# make fasta files 

# install DNA pacakges
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# 
# BiocManager::install(c("Biostrings", "Rsamtools"))

suppressPackageStartupMessages({
  library(Biostrings)
  library(stringr)
  library(pbapply)
  library(tictoc)
  library(assertthat)
})

# rcrs reference genome
rcrs = Biostrings::readDNAStringSet("data/rcrs.fasta")

genetic_data = read.csv('processed_data/genetic_samples.csv')

x = assert_that(length(unique(genetic_data$ID)) == nrow(genetic_data))

mutations = genetic_data$mutations_fromRCRS[1]

## To make this faster and cleaner, the for loop within create_reference_genome could be vectorised 
create_reference_genome = function(mutations, reference, sep = "\\, |\\s"){
  require(stringr)
  
  # Splits the character into mutation sequences
  m_list = unlist(strsplit(mutations, sep))
  m_list = trimws(m_list) # remove any white space 
  m_list_length = length(m_list)
  m_list = unique(m_list) # FXR: If there is a duplication it is a typo error
  
  for(i in seq_along(m_list)){
    m = m_list[i]
    
    # FXR: We checked for that sample and in that position is a C. So this mutation is similar to the rCRS and can just be ignored
    if(m == "T15238C"){
      return(reference)
    }
    
    # identify index in reference genome
    index = unlist(str_extract_all(m, "[0-9]+"))
    index = as.numeric(index)
    
    # If the mutation is a substitution, it will have this format
    # Currently this ignore deletions and insertions
    if(str_detect(m, "^[A-Z]{1}[0-9]+[A-Z]+")){
      # capture the nucleotide at the beginning and the end
      nucleotides = str_match(m, "([A-Z]{1})[0-9]+([A-Z]+)")[1,2:3]
      
      locus = reference[[1]][index]
      
      if(nchar(nucleotides[2]) > 1){
        # mask slippage
        reference[[1]][index] = "N"
        next()
      }
      
      # check that you have found the right locus
      if(as.character(locus[1]) == nucleotides[1]){
        reference[[1]][index] = nucleotides[2]
        } else {
          stop("Leading nucleotide is not what is expected")
        }  
    }
    
    ## If the mutation is a deletion (e.g.) 185A mask it
    if(str_detect(m, "^[0-9]+[A-Z]{1}")){
      ## Mask deletions
      reference[[1]][index] = "N"
      
    }
    
    ## If the mutation is an insertion e.g. 309.1CT mask it
    if(str_detect(m, "^[0-9]+\\.[0-9][A-Z]+")){
      ## Mask insertions
      reference[[1]][index] = "N"
    }
  }
  ## Mask problem regions of the genome, 
  reference[[1]][16182] = "N"
  reference[[1]][16183] = "N"
  reference[[1]][16519] = "N"
  
  reference
}

# remove comma seperators
genetic_data$formatted_mutations = str_replace_all(genetic_data$mutations_fromRCRS, 
                                               pattern = ", ", 
                                               replacement = " ")

## There is a duplicate mutation in rows: 18

# for loop: 292.532 sec elapsed
tic("for loop")
reference_genomes = list()
pb = txtProgressBar(min = 0, max = nrow(genetic_data), initial = 0)
stepi = 0
for(i in seq_along(genetic_data$formatted_mutations)){

  stepi = stepi + 1
  setTxtProgressBar(pb, stepi)

  m = genetic_data$formatted_mutations[i]
  reference_genomes[[i]] = create_reference_genome(m, 
                                  reference = rcrs,
                                  sep = " ")
}
close(pb)
toc()

names(reference_genomes) = NULL
fasta_set = do.call(c, reference_genomes)
# rename fasta set
names(fasta_set) = genetic_data$ID
# # duplicated genomes 
# idx_dup = duplicated(fasta_set)
# fasta_set = fasta_set[idx_dup]

x = assert_that(length(unique(names(fasta_set))) == length(fasta_set))
x = assert_that(all(unique(names(fasta_set)) %in% names(fasta_set)))

## Split M and N Haplotypes
fasta_set_m = fasta_set[genetic_data$haplogroup_MN == "M"]
fasta_set_n = fasta_set[genetic_data$haplogroup_MN == "N"]

## For each Haplotype set, add one out group from the other haplotype.

# Outgroup for M is: HGDP00546
fasta_set_m = c(fasta_set_m, fasta_set_n[1])

# Outgroup for N is: HGDP00545
fasta_set_n = c(fasta_set_n, fasta_set_m[1])

## Save fasta file
# M Haplotype
Biostrings::writeXStringSet(fasta_set_m, filepath = 'processed_data/papua_m.fasta')
# ape::write.nexus.data(fasta_set_m, file = "processed_data/papua_m.nex")

# N Haplotype
Biostrings::writeXStringSet(fasta_set_n, filepath = 'processed_data/papua_n.fasta')
# ape::write.nexus.data(fasta_set_n, file = "processed_data/papua_n.nex")

