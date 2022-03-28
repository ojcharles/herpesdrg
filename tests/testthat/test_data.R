context("calling resistance and file handling")

library(herpesdrg)


test_that("variant files return resistance", {
  #----- HCMV
  # read?
  df = call_resistance(infile = system.file("testdata",  "HCMV_A10.tab", package = "herpesdrg"),
                       all_mutations = TRUE, virus = "HCMV")
  expect_equal(nrow(df), 2926)
  
  # calls resistance?
  df = call_resistance(infile = system.file("testdata",  "HCMV_A10.vcf", package = "herpesdrg"),
                       all_mutations = FALSE, virus = "HCMV")
  expect_equal(unique(df$change), c("UL54_D588N", "UL97_C592G", "UL97_H411Y", "UL97_Q126L", "UL97_T409M"))
  
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

test_that("insertions and deletions are handled", {
  
  # RC and insertion and deletion
  # the current method provides an alignment, which is great. But we want to identify stop codons also and return them
  virus = "HSV2"
  all_mutations = TRUE
  infile = system.file("testdata",  "HSV2_tk_RC_frameshift_indel.fasta", package = "herpesdrg")
  df = call_resistance(infile, virus, all_mutations)
  expect_equal(nrow(df), 5)
  expect_equal(grep("frameshift", df$consequence), 1)
  expect_equal(grep("residue_loss_gain", df$consequence), 2)
})


