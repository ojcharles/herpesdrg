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
  
  
  
  if(fa.genome_len == 152222){#if is whole genome, aligned to merlin
    print("fasta whole genome")
    # cat fasta together
    fa.fileConn<-file(fasta_out)
    writeLines(c(ref.text, fa.text), fa.fileConn)
    close(fa.fileConn)
    
    
    #todo do i need add and addfragment? just try the
  }else{
    
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
    
      command = paste("mafft --keeplength --mapout --ep 0.0 --add",
                      fasta_in,
                      fasta_ref,
                      ">", fasta_out)
      
      
      system(command)
      
      
      # }else{
      #   command = paste("mafft --addfragments",
      #                   fasta_in,
      #                   "--6merpair",
      #                   fasta_ref,
      #                   ">", fasta_out)
      # }
      # source(command)
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
