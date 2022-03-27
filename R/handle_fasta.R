#' Handles fasta files
#'
#' Decides whether the provided fasta file is, whole genome Merlin aligned/ assembled,
#' Whole genome other, or only covers a fraction of the genome.
#' If not WG Merlin type, then it aligns using MAFFT --add --keeplength.
#' 
#' Then converts to VCF using snp-sites
#'
#' @param dir the directory containing 
#' @return vcf file and fasta alignment
#' @keywords internal
#' @export
#' 
handle_fasta = function(dir){
  ### fasta alignment
  # define files
  in_fasta = paste0(dir, "/in_query.fasta")
  in_rc_fasta = paste0(dir, "/in_query_rc.fasta")
  ref_fasta = paste0(dir, "/in_ref.fasta")
  in_msa_fasta = paste0(dir, "/in_query_msa.fasta")
  in_rc_msa_fasta = paste0(dir, "/in_query_rc_msa.fasta")
  out_msa = paste0(dir, "/out_msa.fasta")
  out_vcf = stringr::str_replace(out_msa,pattern = ".fasta", replacement = ".vcf")
  out_vcf_deletions = stringr::str_replace(out_msa,pattern = ".fasta", replacement = "_deletions.vcf")
  
  
  print("fasta found")
  fa.text <- readLines(in_fasta)
  fa.header = fa.text[1]
  fa.body <- fa.text[2:length(fa.text)]
  fa.genome_len = as.numeric(sum(nchar(fa.body)))
  ref.text = ape::read.dna(ref_fasta,format = "fasta",as.matrix = T,as.character = T)
  
  
  
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
                    in_fasta,
                    ref_fasta,
                    ">", fasta_out)
    
    
    system(command)
    
  }else{
    # we need to do the forward and reverse alignment, then decide which is best without any prior info
    orig.fa = ape::read.dna(in_fasta,format = "fasta")
    
    # write a RC sequence to file
    rc.fa = ape::complement(orig.fa)
    rc.fa = ape::as.alignment(rc.fa)$seq
    rc.text = c(">RC",rc.fa) 
    writeLines(rc.text , in_rc_fasta)
    

    
    # now run mafft add for both
    command = paste("mafft --keeplength --mapout --add  ",
                    in_fasta,
                    ref_fasta,
                    ">", in_msa_fasta)
    system(command)
    
    command = paste("mafft --keeplength --mapout --add  ",
                    in_rc_fasta,
                    ref_fasta,
                    ">", in_rc_msa_fasta)
    system(command)

    # decide whether the original or RC is the better match
    
    
    # which alignment has fewer gap regions? - lower is better
    t.fa = ape::read.dna(in_msa_fasta, "fasta",as.matrix = T,as.character = T)
    score_original = stringr::str_count(paste0(t.fa[2,],collapse = ""),
                       "[acgt]{1,1000}")
    t.fa = ape::read.dna(in_rc_msa_fasta, "fasta",as.matrix = T,as.character = T)
    score_rc = stringr::str_count(paste0(t.fa[2,],collapse = ""),
                                        "[acgt]{1,1000}")
    optimal_msa = ifelse(score_original <= score_rc, 1,2)
    
    # identify which msa is best, and save this as out.fa. things point to this now
    if(score_original <= score_rc){
      file.copy(in_msa_fasta , out_msa)
      file.copy(stringr::str_replace(in_fasta,pattern = ".fasta", replacement = ".fasta.map") ,
                stringr::str_replace(out_msa,pattern = ".fasta", replacement = ".fasta.map"))
    }else{
      file.copy(in_rc_msa_fasta , out_msa)
      file.copy(stringr::str_replace(in_rc_fasta,pattern = ".fasta", replacement = ".fasta.map") ,
                stringr::str_replace(out_msa,pattern = ".fasta", replacement = ".fasta.map"))
    }
    
    
    # there will be no insertions as mafft --add --keeplength
    # there my be deletions
    # any deletions should either be dodgy sequencing i..e "n" or are frameshifts (we don't expect residue drop in res genes)
    command = paste("snp-sites -vc -o", out_vcf, out_msa)
    system(command)
    
    
    # the mafft .map file tells us which pos has any indels.
    t = readLines(paste0(out_msa, ".map"))
    t = t[3 : length(t)]
    t2 = read.table(text = t, sep = ",")
    names(t2) = c("nt", "query_pos", "ref_pos")
    # now if this is a good assembly then ambiguous positions will be "n", and handled with the variant caller
    # here we simply want to find deletions, within the sequences segment
    deletions_relative_to_reference =  setdiff((min(t2$ref_pos):max(t2$ref_pos)),
            t2$ref_pos)
    
    # now append this to the vcf file
    text = readLines(out_vcf)
    last_vcf_entry = read.table(text = text[length(text)], sep = "\t",colClasses = c("V4" = "character", "V5" = "character"))
    for(del in 1:length(deletions_relative_to_reference)){
      # what is the reference base at this position?
      "1\t47924\t.\tC\tT\t.\t.\t.\tGT\t0\t1"
      text = c(text,
               paste(last_vcf_entry$V1 + del,
                     deletions_relative_to_reference[del],
                     last_vcf_entry[,3],
                     toupper(ref.text[deletions_relative_to_reference[del]]), # reference base
                     "-", # deletion
                     last_vcf_entry[,6],
                     last_vcf_entry[,7],
                     last_vcf_entry[,8],
                     last_vcf_entry[,9],
                     last_vcf_entry[,10],
                     last_vcf_entry[,11], sep = "\t")
      )
    }
    # write the new vcf
    writeLines(text, out_vcf_deletions)
    
    return(out_vcf_deletions)
   
    
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
