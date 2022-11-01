context("displaying visuals")

library(herpesdrg)
library(stringr)

cls <- c("datatables", "htmlwidget")

test_that("Data table is functioning", {
  dat = call_resistance(infile = system.file("testdata",  "HCMV_A10.vcf", package = "herpesdrg") ,all_mutations = FALSE, virus = "HCMV")
  clin_table = make_clin_table(dat)  
  expect_equal(sum(str_count(as.character(clin_table$x$data[1,]),pattern = "Low level")),0)
  expect_equal(sum(str_count(as.character(clin_table$x$data[1,]),pattern = "High level")),7)
  expect_equal(sum(str_count(as.character(clin_table$x$data[2,]),pattern = "In vitro")),2)
  expect_equal(length(attr(clin_table$x,"colnames")),12)
})


test_that("lollipops are correct", {
  #package variables
  dat = call_resistance(infile = system.file("testdata",  "HCMV_A10.vcf", package = "herpesdrg"),all_mutations = FALSE, virus = "HCMV")
  plot = plot_lollipop(dat, f.gene = "UL97" )
  
  expect_equal(length(plot$layers),5)
  expect_equal(class(plot$layers[[1]]$geom)[1], "GeomSegment")
  expect_equal(plot$labels$y, "Mutation Frequency")
  
})
