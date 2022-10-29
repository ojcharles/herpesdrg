#' Handles input file
#'
#' This function handles vcf, varscan2 tab variant & fasta file formats
#' onto handle_fasta().
#' Returns an intermediate data.frame, which contains variant information to be
#' annotated by annotate_variants().
#'
#' @param infile Path to the input file
#' @param global Package object for consistent runtime variables
#' @return An intermediate data.frame
#' @keywords internal
#' @export
#' 
read_input <- function(infile, global){
  #takes .vcf .tab .fasta inputs and returns a standard dataframe for resistance calling
  infile <- as.character(infile)
  ### tab ###
  if(tools::file_ext(infile) == "tab"){
    #if file appears as a varscan tab file
    # delimit
    tab.dat <- utils::read.table(file = infile, header = T, as.is = T, sep = "\t")
    out <- read_varscan_data(tab.dat)
  }
  ### vcf ###
  else if(tools::file_ext(infile) == "vcf"){
    
    #get file VCF version
    text = readLines(infile)
    vcf_format = ""
    if( grepl("VCFv4.1",text[1]) ){ vcf_format = "4.1"}
    if( grepl("VCFv4.2",text[1]) ){ vcf_format = "4.2"}
    #if( grepl("VCFv4.3",text[1]) ){ vcf_format = "4.3"}
    if(vcf_format == ""){ stop("only VCF formats 4.1 -> 4.2 are currently supported! please log a github issue")}
    
    
    if(vcf_format == "4.1"){
      out = parse_vcf_4_1(infile)
    } else if(vcf_format == "4.2"){
      out = parse_vcf_4_2(infile)
    }
    
    
    
    
    
  }
  ### fasta ###
  else if(tools::file_ext(infile) %in% c("fa", "fasta", "fas", "fna")){
    # in each case output is a vcf file, which then gets processed as above into the out data structure.
    query_file_loc = file.path(global$dir, "in_query.fasta")
    ref_file_loc = file.path(global$dir, "in_ref.fasta")
    file.copy(infile, query_file_loc,overwrite = T)
    file.copy(global$path_fasta_file, ref_file_loc,overwrite = T)
    fasta_out = file.path(global$dir, "out_msa.fasta")
    vcf_file = handle_fasta(dir =  global$dir) 
    text <- readLines(vcf_file)
    start <- base::grep('chrom',ignore.case = T, text)
    vcf = utils::read.delim(vcf_file, sep = "\t", as.is = T, skip = start - 1,colClasses = c("character"))
    vcf[,2] = as.numeric(vcf[,2])
    vcf[,10] = as.numeric(vcf[,10])
    vcf[,11] = as.numeric(vcf[,11])
    
    # if has a format column & genotype column, split to extract ref.count, var.count per position
    for(i in 1:nrow(vcf)){#clean up vcf indel format to be as in varscan tab
      ref = vcf$REF[i]
      var = vcf$ALT[i]
      if(nchar(ref) > 1){#if deletion
        out.ref = var
        out.var = ref
        base::substr(out.var, 1, 1) <- "-"
        vcf$REF[i] = out.ref
        vcf$ALT[i] = out.var
      }
      if(nchar(var) > 1){#if insertion
        out.ref = ref
        out.var = var
        base::substr(out.var, 1, 1) <- "+"
        vcf$REF[i] = out.ref
        vcf$ALT[i] = out.var
      }
    }
    
    t.vcf <- data.frame(Position = vcf$POS,
                        Ref = vcf$REF,
                        Var = vcf$ALT,
                        Ref.count = vcf[,10], #diff from vcf proc
                        Var.count = vcf[,11], # diff from vcf proc
                        VarFreq = "100%",
                        Sample = "single run",
                        stringsAsFactors = F)
    out <- t.vcf
    
    # remove all temprary files
    
    

    
  }else{
    stop("Check your variant call file is in .tab, .vcf format \n or check your fasta file has the .fa, .fas or .fasta extension")
    
  }
  
  return(out)
  #remember to update read_Varscan input functions to read a dataframe not a file location
}













parse_vcf_4_1 = function(infile){
  # read vcf
  r = VariantAnnotation::readVcf(infile)
  t = data.frame(r@assays@data@listData)
  
  # positions
  pos_s = r@rowRanges@ranges@start
  
  # adjust ref alt so it matches varscan tab formatting
  ref = as.character(r@fixed$REF)
  var = as.character(as.data.frame(r@fixed$ALT)$value)
  for( i in 1:length(ref) ){#clean up vcf indel format to be as in varscan tab
    r = ref[i]
    v = var[i]
    if(nchar(r) > 1){#if deletion
      out.ref = v
      out.var = r
      base::substr(out.var, 1, 1) <- "-"
      ref[i] = out.ref
      var[i] = out.var
    }
    if(nchar(v) > 1){#if insertion
      out.ref = r
      out.var = v
      base::substr(out.var, 1, 1) <- "+"
      ref[i] = out.ref
      var[i] = out.var
    }
  }
  
  # format for internal handling
  out = data.frame(Position = pos_s,
                   Ref = ref,
                   Var = var,
                   Ref.count = t$RD, #diff from vcf proc
                   Var.count = t$AD, # diff from vcf proc
                   VarFreq = t$FREQ,
                   Sample = "single run",
                   stringsAsFactors = F)
  return(out)
}

parse_vcf_4_2 = function(infile){
  # this is really terrible code
  
  # read vcf
  r = VariantAnnotation::readVcf(infile)
  t = data.frame(r@assays@data@listData)
  
  # positions
  pos_s = r@rowRanges@ranges@start #
  refs = as.character(r@fixed$REF)
  alts = r@fixed$ALT
  var_counts = data.frame(t$AO)
  ref_counts = t$RO
  
  a = data.frame(alts)
  colnames(a) = c("rownum", "NA", "var")
  
  a$pos = 0 ; a$ref = "" ;a$ref_count = 0; a$var_counts = 0
  na = nrow(a)
  for(i in 1:na){
    index = a$rownum[i]
    a[i,4] = pos_s[index]
    a[i,5] = refs[index]
    a[i,6] = ref_counts[index]
    
    # which genotype index am i at this position?
    t1 = a[1:i,1]
    t1 = sum(t1 == index,na.rm = T)
    a[i,7] = var_counts[t1,index]
  }
  
  
  
  # adjust ref alt so it matches varscan tab formatting
  ref = a$ref
  var = a$var
  for( i in 1:length(ref) ){#clean up vcf indel format to be as in varscan tab
    r = ref[i]
    v = var[i]
    if(nchar(r) > 1){#if deletion
      out.ref = v
      out.var = r
      base::substr(out.var, 1, 1) <- "-"
      ref[i] = out.ref
      var[i] = out.var
    }
    if(nchar(v) > 1){#if insertion
      out.ref = r
      out.var = v
      base::substr(out.var, 1, 1) <- "+"
      ref[i] = out.ref
      var[i] = out.var
    }
  }
  
  paste0( round( 100 * a$var_counts / (a$var_counts + a$ref_count), 2) , "%" )
  
  # format for internal handling
  out = data.frame(Position = a$pos,
                   Ref = ref,
                   Var = var,
                   Ref.count = a$ref_count, #diff from vcf proc
                   Var.count = a$var_counts, # diff from vcf proc
                   VarFreq = paste0( round( a$var_counts / (a$var_counts + a$ref_count) , 2) , "%" ),
                   Sample = "single run",
                   stringsAsFactors = F)
  return(out)
}

