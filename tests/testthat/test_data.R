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
  expect_equal(unique(df$change), c("UL54_D588N", "UL97_C592G", "UL97_H411Y", "UL97_T409M"))
  
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

# travis seems a real faff to get working for snp-sites 2.3 as it's r builds are xenial not bionic ubuntu. so ignoring fasta tests for now.
# tested regularly locally. mafft and snp-sites don't need to be tested alone.



