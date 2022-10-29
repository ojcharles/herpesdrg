context("read and parse file formats")

library(herpesdrg)


test_that("VCFv4.1 is handled", {
  
  infile = system.file("testdata",  "HCMV_A10.vcf", package = "herpesdrg")
  df = herpesdrg::read_input(infile)
  df2 = annotate_variants(toannotate = df, global = global)
  
  t = readRDS(system.file("test_references",  "readInput_vcfv4_1.rds", package = "herpesdrg"))
  expect_true( nrow(setdiff(df, t)) == 0 )
})





#package variables
infile = "cmvdrg_test.vcf"
infile = "cmvdrg_test_fs.vcf"
virus = "HCMV"
global = list()
global$res_table = system.file("herpesdrg-db", "herpesdrg-db.tsv", package = "herpesdrg")
global$date <- format(Sys.time(), "%Y-%m-%d")
global$dir = tempdir()
global$virus_genome = utils::read.csv(system.file("", "virus-genome.csv", package = "herpesdrg"),stringsAsFactors = F)
global$genome = global$virus_genome[global$virus_genome$virus == virus,2]
global$path_gff3_file=system.file("ref", paste0(global$genome,".gff3"), package = "herpesdrg")
global$path_fasta_file=system.file("ref", paste0(global$genome,".fasta"), package = "herpesdrg")
global$path_txdb=system.file("ref", paste0(global$genome,".sqlite"), package = "herpesdrg")

#dat1 = read_input(infile, global = global)

#dat2 = annotate_variants(toannotate = dat1, global = global)




#a = VariantAnnotation::readVcfAsVRanges(infile)
#a = as.data.frame(a)


txdb = AnnotationDbi::loadDb(global$path_txdb)
fa = Rsamtools::FaFile(global$path_fasta_file)


vcf = VariantAnnotation::readVcf(infile)
#rd = rowRanges(vcf)
#loc = locateVariants(rd, txdb, VariantAnnotation::CodingVariants())
coding = suppressWarnings(VariantAnnotation::predictCoding(vcf, txdb, seqSource=fa , ignore.strand=FALSE))

coding[mcols(coding)$CONSEQUENCE == "frameshift"]
coding[mcols(coding)$CONSEQUENCE == "synonymous"]





