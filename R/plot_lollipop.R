#' Make gene lollipop plot
#'
#' Produces figures for the web application.
#' Explores spatial location of mutations in resistance genes.
#'
#' @param f.dat resistance data frame from cmvdrg, where all_muts == F
#' @param f.gene Which gene to plot
#' @param resistance_table location of the csv resistance database
#' @return intermediate data.frame with genome level annotation
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' 
#' @export


plot_lollipop <- function(f.dat, f.gene = "UL54", resistance_table = system.file("herpesdrg-db", "herpesdrg-db.tsv", package = "herpesdrg")){
  ## plot mutations in specific genes with resistance mutation.
  # this could be cleaned up a lot.

  # load data
  mut <- f.dat
  
  
  # if sample has any mutations in f.gene
  if(length(base::grep(unique(mut$GENE), pattern = f.gene)) > 0){
    
    drugs = unlist(utils::read.csv(system.file("", "drugs.csv", package = "herpesdrg"),stringsAsFactors = F))
    #manually force value types, should be oneline or sorted in data
    mut$RefCount <- as.numeric(mut$RefCount)
    mut$VarCount <- as.numeric(mut$VarCount) 
    mut$PROTEINLOC <- as.numeric(mut$PROTEINLOC)
    
    # call res mutations from mutations
    mut_res <- mut %>% dplyr::filter(.data$GENEID==f.gene & .data$CONSEQUENCE=="nonsynonymous")
    mut_res$depth <- mut_res$RefCount + mut_res$VarCount 
    resistance <- utils::read.delim(resistance_table, header = TRUE,as.is = TRUE, sep = "\t")
    resistance$change <- paste(resistance$gene,resistance$aa_change,sep="_")
    resistance$aapos <- readr::parse_number(resistance$aa_change)
    resistance <- resistance %>% dplyr::filter(.data$gene== f.gene)
    resistance = reshape2::melt(resistance, measure.vars = drugs)
    resistance = resistance[resistance$value != "",]
    resistance = resistance[!is.na(resistance$value),]
    resistance = resistance %>% dplyr::group_by(.data$change) %>% dplyr::arrange(.data$value) %>% dplyr::top_n(1, .data$value) # reduces data to 1 row per mutation.
    
    #f.gene_colour = "Ganciclovir"
    #d.resall <- resistance[resistance$variable == f.gene_colour,]
    d.resall = resistance
    #d.resall$resistance <-as.character(cut(as.numeric(d.resall$value), c(0,1,1.5,99999999), right=FALSE, labels=c("none", "low", "high")))
    # no need to alter "sus.." or "res.."
    d.resall$resistance = d.resall$value
    d.resall$resistance[d.resall$resistance > 2 & d.resall$resistance != "Resistant" & d.resall$resistance != "Polymorphism"] = "Resistant"
    d.resall$resistance[d.resall$resistance <= 2 & d.resall$resistance != "Resistant" & d.resall$resistance != "Polymorphism"] = "Polymorphism"
    
    
    d.resmuts <- data.frame(x = mut_res$PROTEINLOC,
                            y = readr::parse_number(mut_res$freq),
                            label = mut_res$aachange)
    d.resmuts = d.resmuts[!duplicated(d.resmuts),]
    d.resmuts$resistance = d.resall$resistance[d.resall$aa_change %in% d.resmuts$label]
    t.y <- 1 # updated to hard as now fixed y axis
    g <- ggplot2::ggplot() +
      #must add resistance mut along bottom colour by fold change etc?
      #all res muts
      ggplot2::geom_segment(data = d.resall, ggplot2::aes(x = .data$aapos, xend = .data$aapos, y = -t.y, yend = t.y, colour = .data$resistance)) +
      ggplot2::geom_hline(yintercept = 0) +
      #lollipop called res muts
      ggplot2::geom_segment( data = d.resmuts, ggplot2::aes(x = .data$x, xend = .data$x, y = 0, yend = .data$y, colour = .data$resistance)) +
      ggplot2::geom_point(data = d.resmuts, ggplot2::aes(x = .data$x, y = .data$y), show.legend=FALSE) +
      #ggplot2::geom_text(data = d.resmuts, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), angle = 0, nudge_y = 1)  + 
      ggplot2::scale_colour_manual(values = c("#00BA38", "#F8766D", "#619CFF"),
                                   labels = c("Polymorphism", "Resistant", "Sample Mutation"), 
                                   drop = FALSE) +
      ggrepel::geom_text_repel(data = d.resmuts, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                               point.padding = 0.2,
                               nudge_x = 0,
                               nudge_y = .5,
                               segment.size = 0.2,
                               ylim = c(-Inf, Inf),
                               show.legend = F
      ) +
      ggplot2::theme_light() +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position="bottom"
      ) +
      ggplot2::ylim(-2,100) + # force y axis. alow for db x axis lines
      ggplot2::xlab(paste(f.gene, "AA location")) +
      ggplot2::ylab("Mutation Frequency") +
      ggplot2::guides(colour = ggplot2::guide_legend(title="Mutation Association"))
    
    
    
  }else{
    g <- NULL
  }
  return(g)
  #references
  # https://bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html
}
