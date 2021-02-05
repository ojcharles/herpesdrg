#' Returns a dataframe of the herpesdrg database
#' @return dataframe of the herpesdrg database
#' @export

herpesdrg_data <- function(){
  
  file = system.file("herpesdrg-db", "herpesdrg-db.csv", package = "herpesdrg")
  dat = utils::read.csv(file, stringsAsFactors = F)
  colnames(dat)[[1]] = "mutation_id"
  return(dat)
}
