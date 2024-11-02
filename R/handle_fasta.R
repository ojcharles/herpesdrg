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
handle_fasta = function(dir) {
  ### fasta alignment
  # define files
  in_fasta = paste0(dir, "/in_query.fasta")
  in_rc_fasta = paste0(dir, "/in_query_rc.fasta")
  ref_fasta = paste0(dir, "/in_ref.fasta")
  in_msa_fasta = paste0(dir, "/in_query_msa.fasta")
  in_rc_msa_fasta = paste0(dir, "/in_query_rc_msa.fasta")
  out_msa = paste0(dir, "/out_msa.fasta")
  out_vcf = stringr::str_replace(out_msa, pattern = ".fasta", replacement = ".vcf")
  out_vcf_deletions = stringr::str_replace(out_msa, pattern = ".fasta", replacement = "_deletions.vcf")
  
  
  #print("fasta found")
  fa.text <- readLines(in_fasta)
  fa.header = fa.text[1]
  fa.body <- fa.text[2:length(fa.text)]
  # replace N, n , - with nothing in query fasta
  #fa.body = gsub("n|N|-", "", fa.body)
  #writeLines(c(fa.header, fa.body), in_fasta)
  # continue
  fa.genome_len = as.numeric(sum(nchar(fa.body)))
  ref.text = ape::read.dna(
    ref_fasta,
    format = "fasta",
    as.matrix = T,
    as.character = T
  )
  
  
  
  ### test for binary dependency
  test.software <- Sys.which("snp-sites")
  if (test.software == "") {
    stop(
      "snp-sites is either not installed or not added to PATH. Required for generating VCF from alignments."
    )
  }

  
  
  ### mafft notes
  ## fragment add, i.e. for a gene fasta
  #   % mafft --addfragments fragments --reorder --6merpair --thread -1 existing_alignment > output
  ## full length alignment
  # mafft --thread 8 --reorder --keeplength --mapout --ep 0.0 --add new_sequences input > output
  
  # tested and --add works just fine for the use case, no need for fragment add.
  
  ### test for binary dependency
  test.software <- Sys.which("mafft")
  if (test.software == "") {
    stop(
      "mafft is either not installed or not added to PATH. Required for alignments to reference genome."
    )
  }
  
  
  # if is fairly long i.e. not a pcr product
  if (fa.genome_len > 10000) {
    #print("mafft alignment started")
    command = paste("mafft --keeplength --mapout  --add  ",
                    in_fasta,
                    ref_fasta)
    t = system(command, intern = T, ignore.stderr = T) # by internalising, we remove the many mafft messages
    writeLines(t, out_msa)
    #print("mafft alignment finished")
    
    file.copy(
      stringr::str_replace(in_fasta, pattern = ".fasta", replacement = ".fasta.map") ,
      stringr::str_replace(out_msa, pattern = ".fasta", replacement = ".fasta.map"),
      overwrite = T
    )
    
    
    
  } else{
    # we need to do the forward and reverse alignment, then decide which is best without any prior info
    orig.fa = ape::read.dna(in_fasta, format = "fasta")
    
    # write a RC sequence to file
    rc.fa = ape::complement(orig.fa)
    rc.fa = ape::as.alignment(rc.fa)$seq
    rc.text = c(">RC", rc.fa)
    writeLines(rc.text , in_rc_fasta)
    
    
    
    # now run mafft add for both
    #print("mafft alignment started")
    command = paste("mafft --keeplength --mapout --add  ",
                    in_fasta,
                    ref_fasta)
    t = system(command, intern = T, ignore.stderr = T) # by internalising, we remove the many mafft messages
    writeLines(t, in_msa_fasta)
    
    command = paste("mafft --keeplength --mapout --add  ",
                    in_rc_fasta,
                    ref_fasta)
    t = system(command, intern = T, ignore.stderr = T)
    writeLines(t, in_rc_msa_fasta)
    #print("mafft alignment finished")
    
    # decide whether the original or RC is the better match
    
    # which alignment has fewer gap regions? - lower is better
    t.fa = ape::read.dna(in_msa_fasta,
                         "fasta",
                         as.matrix = T,
                         as.character = T)
    score_original = stringr::str_count(paste0(t.fa[2, ], collapse = ""),
                                        "[acgt]{1,1000}")
    t.fa = ape::read.dna(in_rc_msa_fasta,
                         "fasta",
                         as.matrix = T,
                         as.character = T)
    score_rc = stringr::str_count(paste0(t.fa[2, ], collapse = ""),
                                  "[acgt]{1,1000}")
    
    # identify which msa is best, and save this as out.fa. things point to this now
    if (score_original <= score_rc) {
      file.copy(in_msa_fasta , out_msa, overwrite = T)
      file.copy(
        stringr::str_replace(in_fasta, pattern = ".fasta", replacement = ".fasta.map") ,
        stringr::str_replace(out_msa, pattern = ".fasta", replacement = ".fasta.map"),
        overwrite = T
      )
    } else{
      file.copy(in_rc_msa_fasta , out_msa, overwrite = T)
      file.copy(
        stringr::str_replace(
          in_rc_fasta,
          pattern = ".fasta",
          replacement = ".fasta.map"
        ) ,
        stringr::str_replace(out_msa, pattern = ".fasta", replacement = ".fasta.map"),
        overwrite = T
      )
    }
    
    
  } # now whole genome and pcr products are handled together
  
  # there will be no insertions as mafft --add --keeplength
  # there my be deletions
  # any deletions should either be dodgy sequencing i..e "n" or are frameshifts (we don't expect residue drop in res genes)
  command = paste("snp-sites -v -o", out_vcf, out_msa)
  system(command)
  
  
  # the mafft .map file tells us which pos have indels
  t = readLines(paste0(out_msa, ".map"))
  t = t[3:length(t)]
  df_map_query_ref_pos = utils::read.table(text = t, sep = ",")
  names(df_map_query_ref_pos) = c("nt", "query_pos", "ref_pos")
  # now if this is a good assembly then ambiguous positions will be "n", and handled with the variant caller
  
  
  text = readLines(out_vcf)
  # as we allow non ACGT characters, we now need to wrangle the vcf output
  text = handle_ambiguous_bases_from_snp_sites ( text)
  
  last_vcf_entry = read.table(
    text = text[length(text)],
    sep = "\t",
    colClasses = c("V4" = "character", "V5" = "character")
  )
  
  
  ############ now we need to append to the vcf entry, indels, in the format we get from an actual vcf
  
  # insertions
  mafft_map_index_insertion =  grep("-", df_map_query_ref_pos$ref_pos)
  mafft_map_insertion_query_pos = df_map_query_ref_pos$query_pos[mafft_map_index_insertion]
  mafft_map_insertion_ref_pos = df_map_query_ref_pos$ref_pos[mafft_map_index_insertion - 1]
  mafft_map_num_insertions = length(mafft_map_index_insertion)
  num_insertions = length(mafft_map_index_insertion)
  
  if (num_insertions > 0) {
    # identity the ref position the insertion occurs at
    rp_num = 1
    for(rp in 1:length(mafft_map_insertion_ref_pos)){
      # first entry should always be numeric
      if( grepl("[0-9]",mafft_map_insertion_ref_pos[rp]) ){
        rp_num = mafft_map_insertion_ref_pos[rp]
      }else{
        # this entry is a "-" character, and should be set at the left most numeric value we have
        mafft_map_insertion_ref_pos[rp] = rp_num
      }
    }
    mafft_map_insertion_ref_pos = as.numeric(mafft_map_insertion_ref_pos)
  }
  
  
  # now per insertion location, add a vcf line entry
  if (num_insertions > 0) {
    
    for(i_site in unique(mafft_map_insertion_ref_pos)){
      
      # which var pos map to this ref site?
      i_site_varposs = which(mafft_map_insertion_ref_pos == i_site)
      # which positions in df_map_query_ref_pos are this?
      i_site_varposs_index = mafft_map_index_insertion[i_site_varposs]
      # what is the ref nt?
      i_ref_nt = df_map_query_ref_pos$nt[i_site_varposs_index[1] - 1]
      #what are the nt's inserted?'
      i_var_nts = paste0(i_ref_nt , paste0(df_map_query_ref_pos$nt[i_site_varposs_index],collapse = "") )
      
      text = c(
        text,
        paste(
          last_vcf_entry$V1, # + iter,
          i_site,
          last_vcf_entry[, 3],
          toupper(i_ref_nt),
          toupper(i_var_nts),
          # reference base concat with insertion char
          last_vcf_entry[, 6],
          last_vcf_entry[, 7],
          last_vcf_entry[, 8],
          last_vcf_entry[, 9],
          last_vcf_entry[, 10],
          last_vcf_entry[, 11],
          sep = "\t"
        )
      )
    }
  }
  
  # remove insertions
  last_vcf_entry = read.table(
    text = text[length(text)],
    sep = "\t",
    colClasses = c("V4" = "character", "V5" = "character")
  )
  df_map_query_ref_pos = df_map_query_ref_pos[- mafft_map_index_insertion, ]
  df_map_query_ref_pos$ref_pos = as.numeric(df_map_query_ref_pos$ref_pos)
  
  
  if( nrow(df_map_query_ref_pos) == 0){
    num_deletions = 0
  }else{
    deletions_relative_to_reference =  setdiff((min(df_map_query_ref_pos$ref_pos):max(df_map_query_ref_pos$ref_pos)),
                                               df_map_query_ref_pos$ref_pos)
    num_deletions = length(deletions_relative_to_reference)
  }

  # now append this to the vcf file
  del = 1
  while (del <= num_deletions) {
    iter = 1
    del_ref_pos = deletions_relative_to_reference[del]
    del_ref_prior_pos = del_ref_pos - 1
    del_ref_prior_nt = ref.text[del_ref_prior_pos]
    
    # what is the reference base at this position?
    # "1\t47924\t.\tC\tT\t.\t.\t.\tGT\t0\t1"
    # are there contiguous deletions?
    dels = del
    if( del + 1 <= num_deletions){
      if (deletions_relative_to_reference[del + 1] == del_ref_pos + 1) {
        dels = del:(del + 1)
      }
    }
    if(del + 2 <= num_deletions){
      if( deletions_relative_to_reference[del + 2] == del_ref_pos + 2) {
      dels = del:(del + 2)
      }
    }
    
    
    text = c(
      text,
      paste(
        last_vcf_entry$V1 + del,
        del_ref_prior_pos,
        # we care about the nt before deletion to anchor
        last_vcf_entry[, 3],
        paste(toupper(ref.text[c(del_ref_prior_pos, deletions_relative_to_reference[dels])]), collapse = ""),
        # reference base
        toupper(del_ref_prior_nt),
        # deletion
        last_vcf_entry[, 6],
        last_vcf_entry[, 7],
        last_vcf_entry[, 8],
        last_vcf_entry[, 9],
        last_vcf_entry[, 10],
        last_vcf_entry[, 11],
        sep = "\t"
      )
    )
    del = max(dels) + 1
  }
  
  
  # write the new vcf
  writeLines(text, out_vcf_deletions)
  
  return(out_vcf_deletions)
  
  
  
  
  
}



#' Converts ambiguous ( dimer) base calls in the vcf from snp-sites, to 50% alt frequency
#'
#' This is to better support sanger seq fasta files typically used in hospitals
#' for an ambiguous base (AB), identifies positions with AB, 
#' then replaces the vcf line with the alt choice. and sets it to 50% frequency
#' this assumes you do not have a mixture of two alternate alleles
#' Filters away all over AB
#'
#'
#' @param vcf_text a text list of the vcf file from snp-sites
#' @return the same vcf text but altered AB's
#' @keywords internal
#'
handle_ambiguous_bases_from_snp_sites = function(vcf_text){

  chrom_line = grep("#CHROM", vcf_text)
  
  t_df = utils::read.table(text = vcf_text,colClasses = c("V4"="character",
                                                   "V5"="character")) # will ignore the header automatically
  
  # remove n or N - these cannot be analysed
  # also for now removes trimer calls
  t_df = t_df[!toupper(t_df$V5) %in% c("N", "B", "D", "H", "V"),]
  
  # handle bases can be 1 of 2
  ambig_base_codes = c("R", "Y", "K", "M", "S", "W")
  ambig_base_opt1 = c("A", "C", "G", "A", "C", "A")
  ambig_base_opt2 = c("G", "T", "T", "C", "G", "T")
  
  for( i in 1:length(ambig_base_codes) ){
    ambig_code = ambig_base_codes[i]
    nuc_pos_of_ambig_code = t_df[toupper(t_df$V5) == ambig_code,]$V2
    for(nuc_pos in nuc_pos_of_ambig_code){
      if( t_df[t_df$V2 == nuc_pos,]$V4 == ambig_base_opt1[i] ){
        t_df[t_df$V2 == nuc_pos,]$V5 = ambig_base_opt2[i]
      }else{
        t_df[t_df$V2 == nuc_pos,]$V5 = ambig_base_opt1[i]
      }
      t_df[t_df$V2 == nuc_pos,]$V10 = 1
    }
  }

  t = readr::format_tsv(t_df)
  output = c(vcf_text[1:chrom_line],
             unlist(strsplit(t,"\n"))[-1])
  return(output)
}







