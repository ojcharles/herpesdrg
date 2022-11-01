context("read and parse file formats")

library(herpesdrg)


test_that("VCFv4.1 is handled", {
  infile = system.file("testdata",  "HCMV_A10.vcf", package = "herpesdrg")
  df = herpesdrg::read_input(infile)
  expect_true(nrow(df) == 2501)
  expect_true(sum(grepl("+-",df$Var)) == 0)
  t = readRDS(system.file("test_references",  "readInput_vcfv4_1.rds", package = "herpesdrg"))
  expect_true( nrow(setdiff(df, t)) == 0 )
})


test_that("VCFv4.2 is handled", {
  infile = system.file("testdata",  "HCMV_vcf4.2.vcf", package = "herpesdrg")
  df = herpesdrg::read_input(infile)
  expect_true( nrow(df) == 2 )
  expect_true(sum(grepl("+-",df$Var)) == 0)
})


test_that("VCFv4.1 from snp=sites is handled", {
  infile = system.file("testdata",  "VZV_K25R.vcf", package = "herpesdrg")
  df = herpesdrg::read_input(infile)
  expect_true(nrow(df) == 4)
  expect_true(sum(grepl("+-",df$Var)) == 0)
})


test_that("varscan tab is handled", {
  infile = system.file("testdata",  "HCMV_A10.tab", package = "herpesdrg")
  df = herpesdrg::read_input(infile)
  expect_true(nrow(df) == 3600)
  expect_true(sum(grepl("+-",df$Var)) == 0)
})


