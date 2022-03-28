#' Annotates variant table
#'
#' @param toannotate an intermediate reduced vcf style dataframe
#' @param global internal list with runtime vars
#' @return intermediate data.frame with genome level annotation
#' @keywords internal
#' @export
#' 

annotate_variants = function(toannotate,global){
  check = IRanges::IRanges(start=toannotate$Position, end=toannotate$Position, width=1)
  gr = GenomicRanges::GRanges(seqnames = global$genome ,ranges=check)
  S4Vectors::values(gr) = S4Vectors::DataFrame(id = toannotate$Sample, freq = toannotate$VarFreq, RefCount= toannotate$Ref.count, VarCount= toannotate$Var.count, VarAllele=toannotate$Var)
  varallele = Biostrings::DNAStringSet(toannotate$Var)
  #txdb = makeTxDbFromGFF(file=global$path_gff3_file, format="gff3") # takes 1 sec, save and load.
  #saveDb(txdb, file="inst/ref/NC_006273.2.sqlite")
  txdb = AnnotationDbi::loadDb(global$path_txdb)
  suppressMessages(gn <- GenomicFeatures::genes(txdb))
  vcf = gr
  GenomeInfoDb::seqlevels(txdb) = global$genome
  GenomeInfoDb::seqlengths(vcf) = GenomeInfoDb::seqlengths(txdb)[names(GenomeInfoDb::seqlengths(vcf))]
  GenomeInfoDb::isCircular(vcf) = GenomeInfoDb::isCircular(txdb)[names(GenomeInfoDb::seqlengths(vcf))]
  fa = Rsamtools::FaFile(global$path_fasta_file)
  codvar = suppressMessages(suppressWarnings(VariantAnnotation::locateVariants(vcf, txdb, VariantAnnotation::CodingVariants())))
  coding = suppressWarnings(VariantAnnotation::predictCoding(vcf, varAllele = varallele, txdb, seqSource=fa , ignore.strand=FALSE))
  txids = S4Vectors::values(GenomicFeatures::transcripts(txdb))$tx_name
  names(txids) = S4Vectors::values(GenomicFeatures::transcripts(txdb))$tx_id
  coding_df = as.data.frame(coding)
  coding_df$aachange = paste(coding_df$REFAA,coding_df$PROTEINLOC,coding_df$VARAA,sep="")
  
  # frameshifts need to be labelled - expecting these for stop codons
  good = grep("^[A-Z]{1}$",coding_df$VarAllele)
  which.drop.gain.aa = grep("[-+][A-Z]{3}",coding_df$VarAllele)
  good = c(good,which.drop.gain.aa)
  # anything else is a frameshift
  which.fs = setdiff(1:nrow(coding_df), good)
  coding_df$aachange[which.fs] = paste0(coding_df$REFAA[which.fs],
         coding_df$aachange[which.fs],
         "frameshift")
  coding_df$change = paste(coding_df$GENEID,coding_df$aachange,sep="_")
  coding_df$CONSEQUENCE = as.character(coding_df$CONSEQUENCE)
  coding_df$CONSEQUENCE[which.fs] = "frameshift"
  coding_df$CONSEQUENCE[which.drop.gain.aa] = "residue_loss_gain"
  
  # reduce fluff data
  coding_df = coding_df[,c("change", "freq","start", "strand", "REFCODON", "VARCODON", "RefCount", "VarCount", "CONSEQUENCE")]
  names(coding_df) = c("change", "freq","genome_pos", "strand", "ref_codon", "var_codon", "ref_count", "var_count", "consequence")
  
  return(coding_df)
}


