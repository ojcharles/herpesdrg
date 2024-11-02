context("calling resistance and file handling")

library(herpesdrg)


test_that("variant files return resistance", {
  #----- HCMV
  # read?
  df = call_resistance(infile = system.file("testdata",  "HCMV_A10.tab", package = "herpesdrg"),
                       all_mutations = TRUE, virus = "HCMV")
  expect_equal(nrow(df), 2981)
  expect_equal(nrow(df[! is.na(df$mutation_id),]) , 6)
  
  # calls resistance?
  df = call_resistance(infile = system.file("testdata",  "HCMV_A10.vcf", package = "herpesdrg"),
                       all_mutations = FALSE, virus = "HCMV")
  expect_equal(unique(df$change), c("UL54_709_frameshift", "UL54_883_frameshift", "UL54_D588N", "UL97_C480F", "UL97_C592G", "UL97_H411Y", "UL97_Q126L", "UL97_T409M" ))
  
  #----- HSV1
  # read?
  df = call_resistance(infile = system.file("testdata",  "HSV1_F716L.vcf", package = "herpesdrg"), 
                       all_mutations = TRUE, virus = "HSV1")
  expect_equal(nrow(df), 5)
  
  # calls resistance?
  df = call_resistance(infile = system.file("testdata",  "HSV1_F716L.vcf", package = "herpesdrg"), 
                       all_mutations = FALSE, virus = "HSV1")
  expect_equal(unique(df$change), "UL30_F716L")
  
  #----- HSV2
  df = call_resistance(infile = system.file("testdata",  "HSV2_Q34H.vcf", package = "herpesdrg"), 
                       all_mutations = TRUE, virus = "HSV2")
  expect_equal(nrow(df), 4)
  
  # calls resistance?
  df = call_resistance(infile = system.file("testdata",  "HSV2_Q34H.vcf", package = "herpesdrg"), 
                       all_mutations = FALSE, virus = "HSV2")
  expect_equal(unique(df$change), "UL30_Q34H")
  
  
  #----- VZV
  df = call_resistance(infile = system.file("testdata",  "VZV_K25R.vcf", package = "herpesdrg"), 
                       all_mutations = TRUE, virus = "VZV")
  expect_equal(nrow(df), 4)
  
  # calls resistance?
  df = call_resistance(infile = system.file("testdata",  "VZV_K25R.vcf", package = "herpesdrg"), 
                       all_mutations = FALSE, virus = "VZV")
  expect_equal(unique(df$change), "ORF36_K25R")
  
  #----- HHV6b
  # todo

})



test_that("insertions and deletions are handled with vcf", {
  
  # RC and insertion and deletion
  # the current method provides an alignment, which is great. But we want to identify stop codons also and return them
  virus = "HCMV"
  all_mutations = T
  infile = system.file("testdata",  "HCMV_frameshift_residueloss_gain.vcf", package = "herpesdrg")
  df = call_resistance(infile, virus, all_mutations)
  expect_equal(nrow(df), 5)
  expect_equal( sum(grepl("frameshift", df$consequence)), 1)
  expect_equal(sum(grepl("residue", df$change)), 2)
  
  
  
})



test_that("fasta PCR products that need to be reverse complented are", {
  
  # this input DNA sequence is the wrong orientation
  # mafft should align it both verbatim and as a RC, identify which is the optimal alignment -> return the optimal alignments data
  # if we have not implemented RC then this will provide essentially random mutations in random genes scattered across the genome
  infile = system.file("testdata",  "HSV2_tk_RC.fasta", package = "herpesdrg")
  virus = "HSV2"
  all_mutations = TRUE
  df = call_resistance(infile, virus, all_mutations)
  expect_equal(nrow(df), 3)
  expect_equal(unique(df$gene), "UL23")
  expect_equal(unique(df$change), "UL23_G39E") # there is only a single nt change
  
})






# test_that("insertions are handled", {
#   
#   # RC and insertion and deletion
#   # the current method provides an alignment, which is great. But we want to identify stop codons also and return them
#   virus = "HSV2"
#   all_mutations = TRUE
#   infile = system.file("testdata",  "HSV2_tk_RC_frameshift_insertion.fasta", package = "herpesdrg")
#   df = call_resistance(infile, virus, all_mutations)
#   expect_equal(nrow(df), 4)
#   expect_equal(grep("frameshift", df$change), 1)
#   expect_equal(grep("frameshift", df$consequence), 1)
# })

test_that("insertions and deletions are handled with fasta", {
  
  # RC and insertion and deletion
  # the current method provides an alignment, which is great. But we want to identify stop codons also and return them
  virus = "HSV2"
  all_mutations = T
  infile = system.file("testdata",  "HSV2_tk_RC_frameshift_indel.fasta", package = "herpesdrg")
  df = call_resistance(infile, virus, all_mutations)
  expect_equal(length(grep("frameshift", df$consequence)), 1)
  expect_equal(length(grep("_residue_loss", df$change)), 1)
  
  
  # example where indels are great in number, and occur at the end of the genome too.
  # previously logic broke when number of indels was not a multiple of 3.
  # this broke some logic looking ahead for if the del position was part of contiguous set.
  virus = "VZV"
  all_mutations = T
  infile = system.file("testdata",  "VZV_many_indel_and_at_end.fasta", package = "herpesdrg")
  df = call_resistance(infile, virus, all_mutations)
  expect_equal(length(grep("indel", df$consequence)), 31)
  expect_equal(length(grep("_residue_loss", df$change)), 29)
  
})




test_that("vcf 4.2 with only allele freq not counts is handles", {
  
  virus = "HCMV"
  all_mutations = T
  infile = system.file("testdata",  "HCMV_vcf4.2_only_freq_not_counts.vcf", package = "herpesdrg")
  df = call_resistance(infile, virus, all_mutations)
  expect_equal(nrow(df), 146)
})


