context("displaying visuals")

library(herpesdrg)
library(stringr)

cls <- c("datatables", "htmlwidget")

test_that("Data table is functioning", {
  dat = call_resistance(infile = system.file("testdata",  "HCMV_A10.tab", package = "herpesdrg") ,all_mutations = FALSE, virus = "HCMV")
  clin_table = make_clin_table(dat)  
  expect_equal(sum(str_count(as.character(clin_table$x$data[1,]),pattern = "Low level")),1)
  expect_equal(sum(str_count(as.character(clin_table$x$data[2,]),pattern = "good, in vitro")),3)
  expect_equal(length(attr(clin_table$x,"colnames")),12)
})


test_that("lollipops are correct", {
  #package variables
  global = list()
  global$res_table = system.file("herpesdrg-db", "herpesdrg-db.csv", package = "herpesdrg")
  #create unique session folder
  global$date <- format(Sys.time(), "%Y-%m-%d")
  global$dir = ""
  
  dat = call_resistance(infile = system.file("testdata",  "HCMV_A10.tab", package = "herpesdrg"),all_mutations = FALSE, virus = "HCMV")
  plot = plot_lollipop(dat, f.gene = "UL54", global = global)
  
  expect_equal(length(plot$layers),5)
  expect_equal(class(plot$layers[[1]]$geom)[1], "GeomSegment")
  expect_equal(plot$labels$y, "Mutation Frequency")
  
})
