#' Runs the HSV Resistance Phenotyping Shiny Application
#' 
#' @import shiny
#' 
#' @export

runShinyHSV = function() {
  appDir <- system.file("shiny-app", package = "hsvdrg")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `hsvdrg`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}