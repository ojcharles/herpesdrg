#' Runs the HSV Resistance Phenotyping Shiny Application
#' 
#' @import shiny
#' 
#' @export

ShinyHerpes = function() {
  appDir <- system.file("shiny-app", package = "herpesdrg")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `herpesdrg`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}