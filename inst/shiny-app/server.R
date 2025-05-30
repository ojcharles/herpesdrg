
# # #debug
# inFile = list()
# inFile$datapath =  system.file("testdata", "A10.vcf", package = "cmvdrg")
# input = list()
# input$virus = "HCMV"
# f.infile = as.character(inFile$datapath)
# session = list()
# session$token = "test"
# # end debug


shinyServer(function(input, output, session) {
  
  global <- reactive({
    

    virus = input$virus
    
    global = list()
    global$res_table = system.file("herpesdrg-db", "herpesdrg-db.tsv", package = "herpesdrg")
    global$date <- format(Sys.time(), "%Y-%m-%d")
    global$dir = tempdir()
    global$virus_genome = utils::read.csv(system.file("", "virus-genome.csv", package = "herpesdrg"),stringsAsFactors = F)
    global$genome = global$virus_genome[global$virus_genome$virus == virus,2]
    global$path_gff3_file=system.file("ref", paste0(global$genome,".gff3"), package = "herpesdrg")
    global$path_fasta_file=system.file("ref", paste0(global$genome,".fasta"), package = "herpesdrg")
    global$path_txdb=system.file("ref", paste0(global$genome,".sqlite"), package = "herpesdrg")
  
    return(global)
  })
  
  
  const_restable = reactive({ 
    d = utils::read.delim(global()$res_table, header = TRUE,na.strings = c("", "NA"), stringsAsFactors = F, sep = "\t") 
    d = d[d$status == "A",]
    })
  
  
  
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
    dat1 <- herpesdrg::read_input(inFile$datapath, global())
    
    ### annotate variants
    # used in optional processes, i.e. identifying syn / nonsynonymous mutations in resistance genes.
    dat2 <- herpesdrg::annotate_variants(dat1, global())
    
    
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
    dat3 <- herpesdrg::add_resistance_info(f.dat = dat2, resistance_table=global()$res_table, all_muts = F, virus = input$virus)
    
    return(dat3)
  })
  

  
  
  
  
  ########################################### Outputs#######################################  
  output$vcf.table_clin <- DT::renderDataTable({
    if(is.null(vcf.d.res())){
      return( NULL )}
    if(nrow(vcf.d.res()) == 0){
      return( data.frame(msg = "no resistance data observed") )}
    # does data contain mutations with > 10% freq?
    freqs = stringr::str_extract(vcf.d.res()$freq, "[0-9]{1,3}.[0-9]{0,3}")
    freqs = as.numeric(gsub("%", "",freqs))
    if ( length( freqs[freqs > 10.00]) == 0 ){
      return( data.frame(msg = "no resistance data above 10% frequency observed") )
    }
    out = herpesdrg::make_clin_table(vcf.d.res())
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
    #                   Brivudin = dat$Brivudin, Penciclovir = dat$Penciclovir,
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
    resistance = const_restable()
    resistance = reshape2::melt(resistance, measure.vars = colnames(resistance[,6:17]))
    resistance = resistance[resistance$value != "",]
    resistance = resistance[!is.na(resistance$value),]
    resistance$value = stringr::str_replace(resistance$value, ">", "")
    resistance$value = stringr::str_replace(resistance$value, ",.{1,5}", "")
    #resistance = resistance[resistance$value > 1,]
    
    #input$dbmetric_virus = c("HCMV", "HSV1")
    #input$dbmetric_drug = c("Aciclovir", "letermovir")
    resistance_filt_dbmetric = resistance[resistance$virus %in% input$dbmetric_virus,]
    
    # collapse by gene, removes concept of virus.
    resistance_filt_dbmetric = reshape2::dcast(resistance_filt_dbmetric, gene ~ variable,fun.aggregate = length)
    resistance_filt_dbmetric = reshape2::melt(resistance_filt_dbmetric, id.vars = "gene")
    resistance_filt_dbmetric$Number_of_entries= resistance_filt_dbmetric$value
    resistance_filt_dbmetric$Drug = resistance_filt_dbmetric$variable
    
    resistance_filt_dbmetric = resistance_filt_dbmetric[resistance_filt_dbmetric$Drug %in% input$dbmetric_drug,]
    resistance_filt_dbmetric = resistance_filt_dbmetric[resistance_filt_dbmetric$Number_of_entries != 0,]
    
    {
     g = ggplot2::ggplot(data = resistance_filt_dbmetric) +
      ggplot2::geom_tile(ggplot2::aes(x = .data$Drug, y = .data$gene, fill = .data$Number_of_entries)) +
      ggplot2::theme_classic() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,vjust = 0.5)) +
      ggplot2::xlab("Drug") +
      ggplot2::ylab("Gene") +
      ggplot2::geom_label(ggplot2::aes(x = .data$Drug, y = .data$gene, label = .data$Number_of_entries))
      #ggplot2::guides(fill = ggplot2::guide_legend(title="Number of entries"))
    }
    return(g)
  })
  
  ### lollipops
  
  ## ul23 ##
  output$vcf.plot.lollipop.UL23 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- herpesdrg::plot_lollipop(mut, f.gene = "UL23")
    return(g)
  })
  # ## ul27 ##
output$vcf.plot.lollipop.UL27 <- renderPlot({
  mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
  g <- herpesdrg::plot_lollipop(mut, f.gene = "UL27")
  return(g)
})
  # ## ul30 ##
    output$vcf.plot.lollipop.UL30 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- herpesdrg::plot_lollipop(mut, f.gene = "UL30")
    return(g)
  })
    # ## ul51 ##
    output$vcf.plot.lollipop.UL51 <- renderPlot({
      mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
      g <- herpesdrg::plot_lollipop(mut, f.gene = "UL51")
      return(g)
    })
  # ## ul54 ##
  output$vcf.plot.lollipop.UL54 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- herpesdrg::plot_lollipop(mut, f.gene = "UL54")
    return(g)
  })
  # ## ul56 ##
  output$vcf.plot.lollipop.UL56 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- herpesdrg::plot_lollipop(mut, f.gene = "UL56")
    return(g)
  })
  # ## ul89 ##
  output$vcf.title.UL89 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(base::grep(unique(mut$Gene), pattern = "UL89", session)) > 0){
      paste("Resistance mutations in UL89 Gene")
    }
  })
  output$vcf.plot.lollipop.UL89 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- herpesdrg::plot_lollipop(mut, f.gene = "UL89")
    return(g)
  })
  # ## ul97 ##
  output$vcf.title.UL97 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(base::grep(unique(mut$Gene), pattern = "UL97", session)) > 0){
      paste("Resistance mutations in UL97 Gene")
    }
  })
  output$vcf.plot.lollipop.UL97 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- herpesdrg::plot_lollipop(mut, f.gene = "UL97")
    return(g)
  })
  # ## ORF28 ##
  output$vcf.title.ORF28 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(base::grep(unique(mut$Gene), pattern = "ORF28", session)) > 0){
      paste("Resistance mutations in ORF28 Gene")
    }
  })
  output$vcf.plot.lollipop.ORF28 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- herpesdrg::plot_lollipop(mut, f.gene = "ORF28")
    return(g)
  })
  # ## ORF36 ##
  output$vcf.title.ORF36 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(base::grep(unique(mut$Gene), pattern = "ORF36", session)) > 0){
      paste("Resistance mutations in ORF36 Gene")
    }
  })
  output$vcf.plot.lollipop.ORF36 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- herpesdrg::plot_lollipop(mut, f.gene = "ORF36")
    return(g)
  })
  # ## U38 ##
  output$vcf.title.U38 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(base::grep(unique(mut$Gene), pattern = "U38", session)) > 0){
      paste("Resistance mutations in U38 Gene")
    }
  })
  output$vcf.plot.lollipop.U38 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- herpesdrg::plot_lollipop(mut, f.gene = "U38")
    return(g)
  })
  # ## U69 ##
  output$vcf.title.U69 <- renderText({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    if(length(base::grep(unique(mut$Gene), pattern = "U69", session)) > 0){
      paste("Resistance mutations in U69 Gene")
    }
  })
  output$vcf.plot.lollipop.U69 <- renderPlot({
    mut <- data.frame(vcf.d.res(),stringsAsFactors = F)
    g <- herpesdrg::plot_lollipop(mut, f.gene = "U69")
    return(g)
  })
  



  
  
  
  
  ##### misc small outputs & downloads 
  output$browse_dbtable = DT::renderDT({
    d = const_restable()
    d$ref_title = paste0( substr(d$ref_titl, 1, 20), "..." )
    d$note = paste0( substr(d$note, 1, 30), "..." )
    d = d[d$status == "A",]
    #d$ref_title = paste0("<a href='",d$ref_title,"'>",d$ref_link,"</a>")
    #d$ref_title = sprintf(paste0('<a href="', URLdecode(d$ref_link),'" target="_blank">', d$ref_title ,'</a>'))
    d$ref_title = NULL
    d$ref_doi = NULL
    d$status = NULL
    colnames(d)[2:5] = c("virus_short_name", "genotype_name", "gene_name", "amino_acid_change")
    d[,-1]
    },
    , options = list(
      autoWidth = TRUE,
      #scrollX = TRUE,
      #scrollY = TRUE,
      dom = 'Bfrtip',
      buttons = c('csv', 'excel'),
      server = FALSE
      ),
    extensions = 'Buttons',
    filter = c('top')
  )
  
  output$vcf.ip <- renderText({
    print(session$request$REMOTE_ADDR)
  })
  
  # about db metric
  output$about_db_total = renderText({
    paste0("Total number of entries is ", as.character(nrow(const_restable())))
  })
  
  output$about_db_unique = renderText({
    uniques = unique( paste0(const_restable()$virus, const_restable()$gene, const_restable()$aa_change) )
    paste0("Number of unique mutations is ", as.character(length(uniques)))
  })
  
  output$about_db_most_recent = renderText({
    dates = as.Date(const_restable()$created_date, format = "%d/%m/%Y")
    most_recent_entry = max(dates)
    paste0("Most recent entry was created ", as.character(most_recent_entry))
  })
  
  # download buttons
  output$vcf.o.res1 <- renderUI({
    req(vcf.d.res())
    downloadButton("vcf.o.res", "download resistance data")
  })
  
  output$vcf.o.all1 <- renderUI({
    req(vcf.d.res())
    downloadButton("vcf.o.all", "download all variants data")
  })
  
  
  
  # debug in browser app mode
  output$vcf.o.res <- downloadHandler(
    filename = function(){paste(global()$date, "_resmuts.csv", sep="")},
    content = function(filename){
      dat <- vcf.d.res()
      utils::write.csv(x = dat, file = filename, row.names = F)
    }
  )
  
  # all mutations, synonymous and non synonymous
  output$vcf.o.all <- downloadHandler(
    filename = function(){paste(global()$date, "_allmuts.csv", sep="")},
    content = function(filename){
      dat <- vcf.d.all()
      dat = herpesdrg::add_resistance_info(f.dat = dat, resistance_table=global()$res_table, all_muts = T, virus = input$virus)
      utils::write.csv(x = dat, file = filename, row.names = F)
    }
  )

})
