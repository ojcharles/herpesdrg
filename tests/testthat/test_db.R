context("Database Integrity")

# this will fail if:
# gene names are not perfect
# virus names are not perfect
# wt AA's are not the same as ref strains
# theres any kink in the reference data

library(herpesdrg)
library(stringr)
library(ape)
library(AnnotationDbi)


dat = utils::read.delim(system.file("herpesdrg-db", "herpesdrg-db.tsv", package = "herpesdrg"),sep = "\t",stringsAsFactors = F,quote = "")
dat$loc = as.integer(str_extract_all(dat$aa_change, pattern = "[0-9]{1,6}", simplify = T )[,1])
dat$from = stringr::str_extract_all(dat$aa_change, pattern = "[A-Z]", simplify = T )[,1]
dat$to = stringr::str_extract_all(dat$aa_change, pattern = "[A-Z]", simplify = T )[,2]
dat = dat[!base::grepl("del",dat$aa_change),] # ignore deletions

genomes = utils::read.csv(system.file("", "virus-genome.csv", package = "herpesdrg"),stringsAsFactors = F)

test_that("Databse Wild Type Integrity", {
  
  out = ""
  
  for(virus in unique(dat$virus)){
    if(!virus %in% genomes$virus){next} # skip if not 
    t1 = dat[dat$virus == virus & dat$status == "A",] 
    for(gene in unique(t1$gene)){
      t2 = t1[t1$gene == gene,]
      # extract sequence of relevant gene
      genome = genomes[genomes$virus == virus,2]
      txdb=system.file("ref", paste0(genome,".sqlite"), package = "herpesdrg")
      txdb <- AnnotationDbi::loadDb(txdb)
      cds = data.frame(GenomicFeatures::cdsBy(txdb, "gene")) # extract coding regions - gene includes introns
      cds = cds[cds$group_name == gene,]
      inseq = ape::read.dna(system.file("ref",paste0(genome,".fasta"), package = "herpesdrg"),format = "fasta",as.character = T,as.matrix = T)
      seq = ""
      for(r in 1:nrow(cds)){ # for each cds 
      start = cds$start[r]
        end = cds$end[r]
        seq = c(seq,inseq[,start:end])
      }
      seq = ape::as.DNAbin(seq[-1])
      if(unique(cds$strand) %in% "-"){ # will be equal for both
        seq = ape::complement(seq)
      }
      outseq = suppressWarnings(ape::trans(seq))
      outseq = as.matrix(as.character(outseq))

      # check mutations lign up with fasta
      for(i in 1:length(outseq)){
        if(nrow(t2[t2$loc ==i,]) > 0){
          ref = outseq[i]
          db = unique(t2[t2$loc == i,]$from)
          db = db[grepl("[A-Z]",db)] # remove any indelstops
          if(length(db) > 1){
            out = paste(out,virus, gene,i,ref,db, "\n")
            #stop()
          }
          if(db != ref){
            out = paste(out,virus, gene,i,ref,db , "\n")
            #stop()
          }
        }
      }
    }
  }
  
  expect_equal(nchar(out), 0)
  print(out)
  
})
