#' Annotates variant table
#'
#' @param toannotate an intermediate reduced vcf style dataframe
#' @param global internal list with runtime vars
#' @return intermediate data.frame with genome level annotation
#' @keywords internal
#' @export
#' 

annotate_variants = function(toannotate,global){
  width = nchar(toannotate$Ref) - 1
  check = IRanges::IRanges(start=toannotate$Position, end=toannotate$Position+width)
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
  

  
  coding_df$change = ""
  
  
  
  #### initial cleaning
  # handle subenes, take first CDSID
  which_2_cdsid = which(lengths(coding_df$CDSID) >1 )
  if( length(which_2_cdsid) > 0 ){ for(i in which_2_cdsid){ coding_df$CDSID[i] = unlist(coding_df$CDSID[i])[1] } }
  
  
  ### handle frameshifts

  which_frameshift = which(coding_df$CONSEQUENCE == "frameshift")
  if( length(which_frameshift) > 0 ){
    
    # many pos -> last pos
    for(i in which_frameshift ){
      if( length( unlist(coding_df$PROTEINLOC[i]) ) > 1 ) {
        pos_s = length(unlist(coding_df$PROTEINLOC[i]))
        coding_df$PROTEINLOC[i] = as.integer(unlist(coding_df$PROTEINLOC[i])[pos_s])
      }else{
        coding_df$PROTEINLOC[i] = as.integer(unlist(coding_df$PROTEINLOC[i])[1])
      }
    }
    
    #now generate change string
    coding_df$change[which_frameshift] = paste0(coding_df$GENEID[which_frameshift], "_",
                                                coding_df$REFAA[which_frameshift],
                                                coding_df$PROTEINLOC[which_frameshift],
                                                "_frameshift"
                                                )
  }
  

  
  
  
  # handle multi-residue and deletions
  which_many_changes = which( coding_df$REFAA != coding_df$VARAA & coding_df$change == "")   # will remove these rows at end
  if( length(which_many_changes) > 0 ){
    for( i in which_many_changes){
      refs = unlist(strsplit(as.character(coding$REFAA[i]),""))
      vars = unlist(strsplit(as.character(coding$VARAA[i]),""))
      refcodons = sapply(seq(from=1, to=nchar(coding_df$REFCODON[i]), by=3), function(x) substr(coding_df$REFCODON[i], x, x+2))
      varcodons = sapply(seq(from=1, to=nchar(coding_df$VARCODON[i]), by=3), function(x) substr(coding_df$VARCODON[i], x, x+2))
      pos_s = as.integer(unlist(coding_df$PROTEINLOC[i])) ; pos_s = pos_s[1] : pos_s[length(pos_s)]
      
      # create template row on current mixed SAV data, and input updates per SAV and rbind
      t =  coding_df[i,]
      for( new_entry in 1:length(refs) ){
        t1 = t
        
        new_pos = pos_s[new_entry]
        new_refcodon = refcodons[new_entry]
        new_varcodon = varcodons[new_entry]
        new_refAA = refs[new_entry]
        new_varAA = vars[new_entry]
        if( is.na(new_refAA) ){ new_refAA = "" }
        if( is.na(new_varAA) ){ new_varAA = ""}
        
        # only certain columns need updating
        t1$PROTEINLOC = as.integer(new_pos)
        t1$REFCODON = new_refcodon
        t1$VARCODON = new_varcodon
        t1$REFAA = new_refAA
        t1$VARAA = new_varAA
        t1
        
        
        if( new_entry == length(refs) & length(vars) > length(refs) ){ # insertion
          t1$VARAA = paste0(vars[ length(refs):length(vars)],collapse = "")
          t1$VARCODON = paste0(varcodons[ length(refcodons):length(varcodons)],collapse = "")
          t1$change = paste0(t1$GENEID, "_", t1$REFAA, t1$PROTEINLOC, "_residue_gain")
          t1$CONSEQUENCE = "indel"
          
        }else if( nchar(new_refAA) == 1 & nchar(new_varAA) == 1 ){ # SAV
          t1$change = paste0(t1$GENEID, "_", t1$REFAA, t1$PROTEINLOC, t1$VARAA)
          if(t1$REFCODON == t1$VARCODON){t1$CONSEQUENCE = "none"}
        
        }else if( nchar(new_refAA) == 1 & new_varAA == "" ){ # deletion
          t1$VARAA = ""
          t1$VARCODON = ""
          t1$change = paste0(t1$GENEID, "_", t1$REFAA, t1$PROTEINLOC, "_residue_loss")
          t1$CONSEQUENCE = "indel"
        }else{ stop ("wha?")}
        
        t1
        
        coding_df = rbind(coding_df , t1)
        
      }
    }
    # now we have appended the well formatted SAV rows, we remove the multipleSAV rows
    coding_df = coding_df[-which_many_changes,]
    
  }
  
  
  
  # # handle residue insertions
  # which_aa_gain = which(nchar(coding_df$REFAA) != nchar(coding_df$VARAA) & coding_df$change == "" )
  # if( length(which_aa_gain) > 0){
  #   g = coding_df$GENEID[which_aa_gain]
  #   r = coding_df$REFAA[which_aa_gain]
  #   v = coding_df$VARAA[which_aa_gain]
  #   pos_s = c()
  #   for(i in which_aa_gain){ pos_s = c(pos_s, unlist(coding_df$PROTEINLOC[i])[1] ) }
  #   coding_df$change[which_aa_gain] = paste0(g, "_",r, pos_s, v, "_residue_gain")
  # }
  

  
  
  

  
  
  
  

  
  
  # clean
  coding_df = coding_df[coding_df$CONSEQUENCE != "none",]
  coding_df$PROTEINLOC = as.integer(unlist(coding_df$PROTEINLOC))
  coding_df$CDSID = as.integer(unlist(coding_df$CDSID))
  coding_df = coding_df[order(coding_df$start),]
  which_SAV = which(coding_df$change == "")
  if( length(which_SAV) > 0 ){
    coding_df$change[which_SAV] = paste(coding_df$GENEID[which_SAV],
                                        "_", 
                                        coding_df$REFAA[which_SAV],
                                        coding_df$PROTEINLOC[which_SAV],
                                        coding_df$VARAA[which_SAV],sep="")
    
  }
  
  
  

  # reduce fluff data
  coding_df = coding_df[,c("change", "freq","start", "strand", "REFCODON", "VARCODON", "RefCount", "VarCount", "CONSEQUENCE")]
  names(coding_df) = c("change", "freq","genome_pos", "strand", "ref_codon", "var_codon", "ref_count", "var_count", "consequence")
  
  rownames(coding_df) = NULL
  
  return(coding_df)
}


