
# # #debug
# inFile = list()
# inFile$datapath =  system.file("testdata", "A10.vcf", package = "cmvdrg")
# f.infile = as.character(inFile$datapath)
# session = list()
# session$token = "test"
# # end debug


shinyServer(function(input, output, session) {
  global = list()
  global$res_table = system.file("db", "hsvdrg-db.csv", package = "hsvdrg")
  global$date <- format(Sys.time(), "%Y-%m-%d")
  global$dir = ""
  global$genome = genome="NC_001806.2"
  global$path_gff3_file=system.file("ref", "NC_001806.2.gff3", package = "hsvdrg")
  global$path_fasta_file=system.file("ref", "NC_001806.2.fasta", package = "hsvdrg")
  global$path_txdb=system.file("ref", "NC_001806.2.sqlite", package = "hsvdrg")
  

  
  
  
  ###************  processes  ****************###  
  
  
  vcf.d.all <- reactive({
    # Returns annotated variants dataframe
    # in  - any of the accepted file formats from file load
    # out - all mutants table
    
    inFile <- input$vcf.file
    if (is.null(inFile)){
      return(NULL)}
    
    
    # handle input .tab .vcf .fasta
    # returns minimum intermediate table
    dat1 <- hsvdrg::read_input(inFile$datapath, global)
    
    ### annotate variants
    # used in optional processes, i.e. identifying syn / nonsynonymous mutations in resistance genes.
    dat2 <- hsvdrg::annotate_variants(f.dat = dat1, global)
    
    
    return(dat2)
  })
  
  vcf.d.res <- reactive({
    # This is the main table used in the shinyapp, dataframe of resistant mutants.e ssentially call_resistance(..., all_mutations = F)
    # out - only resistant data. used by most plots.

    inFile <- input$vcf.file
    if (is.null(inFile)){
      return(NULL)}
    
    dat2 = vcf.d.all()

    ### add res info
    dat3 <- hsvdrg::add_resistance_info(f.dat = dat2, resistance_table=global$res_table, all_muts = F, anecdotal = input$vcf.anecdotal)
    
    return(dat3)
  })
  

  
  
  
  
  ########################################### Outputs#######################################  
  output$vcf.table_clin <- DT::renderDataTable({
    if(is.null(vcf.d.res())){
      return(NULL)}
    out = hsvdrg::make_clin_table(vcf.d.res())  
    return(out)
  })
  
  
  
  output$vcf.table_res <- renderTable({
    # resistance mutation table
    if(is.null(vcf.d.res())){
      return(NULL)
    }else if (nrow(vcf.d.res()) == 0){
      dat <- data.frame(notes = "there was no resistance identified in your sample")
      return
    }else{
    dat <- vcf.d.res()
    #full data is saved and user can DL
    #we really only need to present some of these columns
    # dat = data.frame( gene = dat$gene, aa_change = dat$aa_change, freq = dat$freq, "ref_var_count" = paste(dat$RefCount,	dat$VarCount, sep = "_"),
    #                   ref_pos = dat$start, Aciclovir = dat$Aciclovir, Cidofovir = dat$Cidofovir, Foscarnet = dat$Foscarnet,
    #                   Brivudin = dat$Brivudin, Pencyclovir = dat$Pencyclovir,
    #                   reference = dat$ref_title, ref_link = dat$ref_link, test_method = dat$test_method, test_method_class = dat$tm_class,
    #                   co_gene = dat$co_gene, co_aa = dat$co_aa
    # )
    #dat = dat[,c(28,9,10,8,26,27,29,30,31,32:40,42,43,45,46,47,48)]
    #dat$reference = paste0("<a href='",  dat$ref_link, "' target='_blank'>",substr(dat$reference,1,20),"</a>")
    dat$ref_link= NULL
    #dat <- reshape2::dcast(dat, GENE + AA_CHANGE + freq ~ variable)
    }
    return(dat[,])
  },sanitize.text.function = function(x) x)
  

  output$res.plot.dbheatmap <- renderPlot({
    # plots as a heatmap, the amount of dat we have per gene, per drug.
    # //todo - currently records the number of data points, could record mean value etc
    resistance = utils::read.csv(global$res_table, header = TRUE,na.strings = c("", "NA"), stringsAsFactors = F)
    if (input$vcf.anecdotal == F) {
      resistance = resistance[resistance$tm_class == 'in vitro',]
    }
    resistance = reshape2::melt(resistance, measure.vars = colnames(resistance[,6:10]))
    resistance = resistance[resistance$value != "",]
    resistance = resistance[!is.na(resistance$value),]
    resistance$value = stringr::str_replace(resistance$value, ">", "")
    resistance$value = stringr::str_replace(resistance$value, ",.{1,5}", "")
    resistance = resistance[resistance$value > 1,]
    
    resistance = reshape2::dcast(resistance, gene ~ variable,fun.aggregate = length)
    resistance = reshape2::melt(resistance, id.vars = "gene")
    resistance$Number_of_entries= resistance$value
    resistance$Drug = resistance$variable
    g = ggplot2::ggplot(data = resistance) +
      ggplot2::geom_tile(ggplot2::aes(x = .data$Drug, y = .data$gene, fill = .data$Number_of_entries)) +
      ggplot2::theme_classic() +
      ggplot2::scale_fill_gradient(low="white", high="red") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,vjust = 0.5)) +
      ggplot2::xlab("Drug") +
      ggplot2::ylab("Gene") +
      ggplot2::guides(fill = ggplot2::guide_legend(title="Number of entries"))
    return(g)
  })
  

  # output$vcf.plot.res <- renderPlot({
  #   validate(
  #    need(input$vcf.file != "", "Please load a file")
  #   )
  #   #if (is.null(input$vcf.file)){
  #   #  return(NULL)}
  #   if (nrow(vcf.d.res()) == 0){ # if no resistance data was identified
  #     return(NULL)
  #   }
  #   # plot resistance muts
  #   coding_df_res_resistance <- vcf.d.res()
  #   
  #   if(is.null(coding_df_res_resistance$time_point)){
  #     g <- ggplot( coding_df_res_resistance ,aes(x = "mutations", y=freq,fill=factor(change)))
  #   }else{
  #     g <- ggplot( coding_df_res_resistance ,aes(x=time_point,y=freq,fill=factor(change)))
  #   }
  #   g <- g +
  #     geom_bar(position="dodge",stat="identity") +
  #     scale_fill_discrete(drop=FALSE) +
  #     scale_x_discrete(drop=FALSE) +
  #     theme_bw() +
  #     labs(fill = "Resistance")
  #   # ggsave(plot = g, filename = paste("session",global$date, session$token, "res_table.png", sep = "/"), device = "png")
  #   return(g)
  # })
  
  
  ### lollipops
  
  # ## ul30 ##
  output$vcf.title.UL30 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(base::grep(unique(mut$Gene), pattern = "U30", session)) > 0){
      paste("Resistance mutations in UL30 Gene")
    }
  })
    output$vcf.plot.lollipop.UL30 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- hsvdrg::plot_lollipop(mut, f.gene = "UL30",global = global)
    return(g)
  })
  ## ul23 ##
    output$vcf.title.UL23 <- renderText({
      mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
      if(length(base::grep(unique(mut$Gene), pattern = "UL23")) > 0){
        paste("Resistance mutations in UL23 Gene")
      }
    })
  output$vcf.plot.lollipop.UL23 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- hsvdrg::plot_lollipop(mut, f.gene = "UL23",global = global)
    return(g)
  })

  
  
  
  
  ###************  Outputs  ****************###  
  output$vcf.ip <- renderText({
    print(session$request$REMOTE_ADDR)
  })
  
  # downlaod buttons
  output$vcf.o.res1 <- renderUI({
    req(vcf.d.res())
    downloadButton("vcf.o.res", "download resistance data")
  })
  
  output$vcf.o.all1 <- renderUI({
    req(vcf.d.res())
    downloadButton("vcf.o.all", "download all variants data")
  })
  
  
  
  # debug in brownser app mode
  output$vcf.o.res <- downloadHandler(
    filename = function(){paste(global$date, "_resmuts.csv", sep="")},
    content = function(filename){
      dat <- vcf.d.res()
      utils::write.csv(x = dat, file = filename, row.names = F)
    }
  )
  
  # all mutations, synonymous and non synonymous
  output$vcf.o.all <- downloadHandler(
    filename = function(){paste(global$date, "_allmuts.csv", sep="")},
    content = function(filename){
      dat <- vcf.d.all()
      dat = hsvdrg::add_resistance_info(f.dat = dat, resistance_table=global$res_table, all_muts = T, anecdotal = input$vcf.anecdotal)
      utils::write.csv(x = dat, file = filename, row.names = F)
    }
  )

})
