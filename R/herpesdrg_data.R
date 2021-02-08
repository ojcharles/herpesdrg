#' Returns a dataframe of the herpesdrg database
#' @return dataframe of the herpesdrg database
#' @export

herpesdrg_data <- function(){
  
  file = system.file("herpesdrg-db", "herpesdrg-db.tsv", package = "herpesdrg")
  dat = utils::read.delim(file, stringsAsFactors = F,sep = "\t")
  colnames(dat)[[1]] = "mutation_id"
  return(dat)
}
