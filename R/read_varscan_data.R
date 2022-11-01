#' Handles varscan data
#'
#' reads varscan files & returns a formatted data.frame
#'
#' @param f.df data.frame of read varscan file
#' @return formatted data.frame
#' @keywords internal
#' @export
#' 

read_varscan_data <- function(infile){
  dat <- utils::read.table(file = infile, header = T, as.is = T, sep = "\t")
  dat$id = "single run"
  df.het <- as.data.frame(matrix(unlist(strsplit(dat[,5], split=":")), ncol=6, byrow="T"), stringsAsFactors=F)
  all <- cbind(dat[,1:4], df.het[,1:5], dat[,6],dat[,12]) # shifted from in-house 
  colnames(all)[5]<-"Var_touse"
  colnames(all)[6]<-"something"
  colnames(all)[7]<-"Ref.count"
  colnames(all)[8]<-"Var.count"
  colnames(all)[9]<-"VarFreq"
  colnames(all)[10]<-"StrandFilter"
  colnames(all)[11]<-"Sample"
  all$Var_touse <- NULL
  all$StrandFilter <- NULL
  all$something <- NULL
  all$Chrom <- NULL
  
  # make the ref and var encoding standard
  which_loss = grep("-", all$Var)
  if( length(which_loss) > 0 ){
    r = all$Ref[which_loss]
    v = all$Var[which_loss]
    new_ref = paste0(r,v) ; new_ref = gsub("-","", new_ref)
    new_var = r
  }
  all$Ref[which_loss] = new_ref
  all$Var[which_loss] = new_var
  
  which_ins = grep("\\+", all$Var)
  if( length(which_ins) > 0 ){
    r = all$Ref[which_ins]
    v = all$Var[which_ins]
    new_var = paste0(r,v) ; new_var = gsub("\\+","", new_var)
  }
  all$Var[which_ins] = new_var
  
  return(all)
  }