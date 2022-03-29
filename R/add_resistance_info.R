#' Resistance Genotyping
#'
#' Calls resistance for variants provided, joins the variant and resistance tables by mutational change columns. <GENE_A123B>
#'
#' @param f.dat intermediate-annotated data.frame
#' @param resistance_table the current version of the resistance db in csv format
#' @param all_muts when TRUE all variants passed are returned even if they conferred no resistance
#' @return data.frame of resistance variants
#' @keywords internal
#' @export
#'
add_resistance_info <-
  function(f.dat,
           resistance_table,
           all_muts = FALSE,
           virus) {
    coding_df <- f.dat
    resistance = utils::read.delim(resistance_table, header = TRUE,sep = "\t")
    
    # filter status - records on-revision may be below the data quality we expect, and are flagged.
    resistance = resistance[resistance$status == "A", ]
    resistance = resistance[resistance$virus == virus, ]
    resistance$change <- paste(resistance$gene, resistance$aa, sep = "_")
    
    
    # merge resistance & mutation data
    coding_df_res <- base::merge(
      x = coding_df,
      y = resistance,
      by = "change",
      all.x = T
    )
    
    if (all_muts == F) {
      coding_df_res2  = coding_df_res[grepl("[0-9]",coding_df_res$mutation_id),]
      coding_df_res2 = rbind(coding_df_res2,
                            coding_df_res[grepl("frameshift",coding_df_res$change),])
      coding_df_res = coding_df_res2
    }
    
    
    # annotate any frameshift mutations
    if(length(grep("frameshift", coding_df_res$change)) > 0){
      which.fs = grep("frameshift", coding_df_res$change)
      coding_df_res[which.fs,c("Aciclovir", "Ganciclovir", "Cidofovir", "Brincidofovir", "Pencyclovir")] = "Resistant"
      coding_df_res[which.fs,"note"] = "Frameshifts often result in premature stop codons, which heavily reduce antiviral efficacy"
      #coding_df_res <- cbind(resistance_site,coding_df)
    }
    return(coding_df_res)
    
  }
