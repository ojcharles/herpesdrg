#' Make gene lollipop plot
#'
#' Produces figures for the web application.
#' Explores spatial location of mutations in resistance genes.
#'
#' @param f.dat resistance data frame from cmvdrg, where all_muts == F
#' @param f.gene Which gene to plot
#' @param global Package object for consistent runtime variables
#' @return intermediate data.frame with genome level annotation
#' 
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' 
#' @export


plot_lollipop <- function(f.dat, f.gene = "UL54", global){
  ## plot mutations in specific genes with resistance mutation.
  # this could be cleaned up a lot.
  
  # load data
  mut <- f.dat

  
  # if sample has any mutations in f.gene
  if(length(base::grep(unique(mut$GENE), pattern = f.gene)) > 0){
    
    #manually force value types, should be oneline or sorted in data
    mut$RefCount <- as.numeric(mut$RefCount)
    mut$VarCount <- as.numeric(mut$VarCount) 
    mut$PROTEINLOC <- as.numeric(mut$PROTEINLOC)

    # call res mutations from mutations
    mut_res <- mut %>% dplyr::filter(.data$GENEID==f.gene & .data$CONSEQUENCE=="nonsynonymous")
    mut_res$depth <- mut_res$RefCount + mut_res$VarCount 
    resistance_table=global$res_table
    resistance <- utils::read.csv(resistance_table, header = TRUE,as.is = TRUE)
    resistance$change <- paste(resistance$gene,resistance$aa_change,sep="_")
    resistance$aapos <- readr::parse_number(resistance$aa_change)
    resistance <- resistance %>% dplyr::filter(.data$gene== f.gene)
    resistance = reshape2::melt(resistance, measure.vars = colnames(resistance[,6:14]))
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
      ggplot2::geom_point(data = d.resmuts, ggplot2::aes(x = .data$x, y = .data$y, colour = "Sample Mutations" , size = 8), show.legend=FALSE) +
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



 #old method
# # pass sample gene info with res
# sample.gr <- GRanges("NC_006273.2", IRanges(mut_res$PROTEINLOC,width = 1, names = paste0(mut_res$aachange,sep=" N=",mut_res$depth))) 
# sample.gr$score <- readr::parse_number(mut_res$freq)
# 
# # all resistance muts - to add to plot as smlal black lines
# features_res <- GRanges("NC_006273.2", IRanges(resistance$aapos,width = 1))
# 
# # plot res gene lollipop
# features_res$fill <- "orange"
# sample.gr$color <- sample(c("orange"), length(mut_res$PROTEINLOC), replace=TRUE)
# sample.gr$border <- sample(c("gray80"), length(mut_res$PROTEINLOC), replace=TRUE)  
# 
#
#
# this function takes 10 seconds!!!!!!!!!
#g <- lolliplot(sample.gr, features_res,ylab="freq",xlab= paste(f.gene, "AA location"))
