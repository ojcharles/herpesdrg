#' Handles fasta files
#'
#' Decides whether the provided fasta file is, whole genome Merlin aligned/ assembled,
#' Whole genome other, or only covers a fraction of the genome.
#' If not WG Merlin type, then it aligns using MAFFT --add --keeplength.
#' 
#' Then converts to VCF using snp-sites
#'
#' @param fasta_in Path to the input file
#' @param fasta_out filename of fasta alignment produced
#' @param fasta_ref the Merlin reference fasta file
#' @return vcf file and fasta alignment
#' @keywords internal
#' @export
#' 
handle_fasta = function(fasta_in, fasta_out, fasta_ref){
  ### fasta alignment
  # in - fasta file
  # out - fasta alignment
  
  print("fasta found")
  fa.text <- readLines(fasta_in)
  fa.header = fa.text[1]
  fa.body <- fa.text[2:length(fa.text)]
  fa.genome_len = as.numeric(sum(nchar(fa.body)))
  ref.text = readLines(fasta_ref)
  
  
  
  ### test for binary dependency
  test.software <- Sys.which("snp-sites")
  if(test.software == ""){
    stop("snp-sites is either not installed or not added to PATH. Required for generating VCF from alignments.")
  }
  ####
  
  
  
  ### mafft notes
  ## fragment add, i.e. for a gene fasta
  #   % mafft --addfragments fragments --reorder --6merpair --thread -1 existing_alignment > output
  ## full length alignment
  # mafft --thread 8 --reorder --keeplength --mapout --ep 0.0 --add new_sequences input > output
  
  # tested and --add works just fine for the use case, no need for fragment add.
  
  ### test for binary dependency
  test.software <- Sys.which("mafft")
  if(test.software == ""){
    stop("mafft is either not installed or not added to PATH. Required for alignments to reference genome.")
  }
  ###
  
  # if is fairly long i.e. not a pcr product
  if(fa.genome_len > 10000){
    command = paste("mafft --keeplength --mapout  --add  ",
                    fasta_in,
                    fasta_ref,
                    ">", fasta_out)
    
    
    system(command)
    
  }else{
    # we need to do the forward and reverse alignment, then decide which is best without any prior info
    rc.fa = ape::read.dna(fasta_in,format = "fasta")
    rc.fa = ape::complement(rc.fa) # reverse complement
    rc.fa = ape::as.alignment(rc.fa)$seq
    
    # append this reverse complement sequence to the original input sequence
    fa.text[3] = ">RC"
    fa.text[4] = rc.fa
    writeLines(fa.text , fasta_in)
    
    # now run mafft add for both
    command = paste("mafft --keeplength --mapout  --thread 4 --add  ",
                    fasta_in,
                    fasta_ref,
                    ">", fasta_out)
    system(command)
    
    # decide whether the original or RC is the better match
    t.fa = ape::read.dna(fasta_out, "fasta",as.matrix = T,as.character = T)
    
    # which alignment has fewer gap regions? - lower is better
    score_original = stringr::str_count(paste0(t.fa[2,],collapse = ""),
                       "[acgt]{1,1000}")
    score_rc = stringr::str_count(paste0(t.fa[3,],collapse = ""),
                                        "[acgt]{1,1000}")
    
    # update fasta_out, to either be ref+ provided seq, or ref + revcomp seq
    if(score_original <= score_rc){
      writeLines(c(paste0(">",names(t.fa[1,1])),
                   paste0(t.fa[1,], collapse = ""),
                   paste0(">",names(t.fa[2,1])),
                   paste0(t.fa[2,], collapse = "")), fasta_out)
    }else{
      writeLines(c(paste0(">",names(t.fa[1,1])),
                   paste0(t.fa[1,], collapse = ""),
                   paste0(">",names(t.fa[3,1])),
                   paste0(t.fa[3,], collapse = "")), fasta_out)
    }
    
    
  }


  
  
  ## now we have a fasta file called out.fasta in the session directory
  
  
  ## snp-sites fasta -> vcf. This is great!
  # https://github.com/sanger-pathogens/snp-sites/blob/master/snp-sites.txt
  vcf_out = stringr::str_replace(fasta_out,pattern = ".fasta", replacement = ".vcf")
  command = paste("snp-sites -vc -o", vcf_out, fasta_out)
  system(command)
  
  
  #return vcf file location to pass onto 
  return(vcf_out)
    
     
}
