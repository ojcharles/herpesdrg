#' HSV Resistance Genotyping from command line
#'
#' Command line version of the website application
#' Takes as input a VCF, varscan tab or fasta file.
#' The program assumes variant files are generated relative to Merlin strain.
#' Fasta files if not Whole Genome, or not aligned / assembled relative to Merlin 
#' are Processed using MAFFT & snp-sites.
#' In this case the output files are returned to your working directory.
#'
#' @param infile the input fasta, vcf or varscan.tab file
#' @param virus HCMV HSV1 HSV2 VZV
#' @param all_mutations when FALSE only recognised resistant variants present are returned.
#' @param outdir for fasta input files intermediate alignment fasta & vcf files are generated, this defines the directory they are saved to. "out.fasta" "out.vcf"
#' @return A data.frame containing resistance information for variants identified
#' @export

call_resistance = function(infile = system.file("testdata",  "HSV1_F716L.vcf", package = "herpesdrg"), virus = "HSV1", all_mutations = TRUE, outdir = ""){
  
  #package variables
  global = list()
  global$res_table = system.file("herpesdrg-db", "herpesdrg-db.csv", package = "herpesdrg")
  global$date <- format(Sys.time(), "%Y-%m-%d")
  global$dir = outdir
  global$virus_genome = utils::read.csv(system.file("", "virus-genome.csv", package = "herpesdrg"),stringsAsFactors = F)
  global$genome = genome=global$virus_genome[global$virus_genome$virus == virus,2]
  global$path_gff3_file=system.file("ref", paste0(global$genome,".gff3"), package = "herpesdrg")
  global$path_fasta_file=system.file("ref", paste0(global$genome,".fasta"), package = "herpesdrg")
  global$path_txdb=system.file("ref", paste0(global$genome,".sqlite"), package = "herpesdrg")
  
  
  dat1 = read_input(infile, global = global)
  
  ### annotate variants
  dat2 <- annotate_variants(f.dat = dat1, global = global)
  
  ### add res info
  dat3 <- add_resistance_info(f.dat = dat2, resistance_table=global$res_table, all_muts = all_mutations, virus = virus)
  
  dat3=dat3
  
  return(dat3)
}
