# UI for CMV resistance database
# Author - Charles, OJ
# Date - April 2020
################
# design idea https://community.rstudio.com/t/background-images-in-shiny/12261/3
library(shiny)


# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    # Application title
    titlePanel("HerpesDRG - Herpesvirus Drug Resistance Genotyping"),
    h4("HSV1, HSV2, HCMV, VZV & HHV6b"),
    navlistPanel(widths = c(2,6),
       "Resistance Genotyping",
       tabPanel("Homepage",
                div("The below tool enables researchers to detect herpesvirus resistance mutations in common sequencing file formats. Check the help page for further details."),
                br(),
                h4("Upload a file to begin processing"),
                selectInput("virus", label = "Virus", choices = c("HCMV", "HSV1", "HSV2", "VZV", "HHV6b"),multiple = F),
                fileInput("vcf.file", label = "Accepts VCF, Varscan-tab, FASTA formats",multiple = F, width = 400),
                br(),
                div(style="display:inline-block",uiOutput("vcf.o.res1")),
                div(style="display:inline-block",uiOutput("vcf.o.all1")),
                #div(style="display:inline-block",downloadButton("vcf.o.res", "download resistance data")),
                #div(style="display:inline-block",downloadButton("vcf.o.all", "download all-mutants data")),
                br(),br(),h4("Clinical overview - mutations > 10% frequency"),br(),
                DT::dataTableOutput("vcf.table_clin"),
                br(),br(),h4("Identified resistant mutants"),br(),
                tableOutput("vcf.table_res")
                #plotOutput("vcf.plot.res")

       ), #tabpanel
       

       tabPanel("Genetics",
                p(""),
                br(),
                strong("Graphics showing location of identified resistance mutations, in context to all database mutations"),
                plotOutput("vcf.plot.lollipop.UL23"),
                plotOutput("vcf.plot.lollipop.UL27"),
                plotOutput("vcf.plot.lollipop.UL30"),
                plotOutput("vcf.plot.lollipop.UL51"),
                plotOutput("vcf.plot.lollipop.UL54"),
                plotOutput("vcf.plot.lollipop.UL56"),
                plotOutput("vcf.plot.lollipop.UL89"),
                plotOutput("vcf.plot.lollipop.UL97"),
                plotOutput("vcf.plot.lollipop.ORF28"),
                plotOutput("vcf.plot.lollipop.ORF36"),
                plotOutput("vcf.plot.lollipop.U38"),
                plotOutput("vcf.plot.lollipop.U69")
       ), #tabpanel
       tabPanel("Browse Mutations",
                strong("Browse database entries using the table below:"),
                DT::DTOutput("browse_dbtable")
       ),
       
       
       "Information",
       # tabPanel("Database-entries",
       #  # todo: make some fields mandatory? using shinyjs https://deanattali.com/2015/06/14/mimicking-google-form-shiny/
       #  # todo: handle submission, requires a data arhchitecture first.
       #   sidebarLayout(
       #     sidebarPanel(
       #       h3("Download resDB"),
       #       selectInput("db.db-version", "select db version",
       #                   c("v1.0", "v1.1-beta")),
       #       downloadButton("db.download", "Download Resistance CSV"),
       #       "",
       #       h3("Issues?"),
       #       textInput("db.gene", "Gene", placeholder = "e.g. UL54"),
       #       textInput("db.aamut", "Amino Acid Location & mutation",placeholder = "e.g. A571G"),
       #       selectInput("db.drugres", "Drug Resistance to assign / update",
       #                   c("", "GCV","CDV", "FOSarnate", "BCV", "LENtermovir", "TOMoglovir", "MARibavir", "other-in note"),multiple = T),
       #       textAreaInput("db.note", label = "Include in this note, what you would like to see updated",
       #                     placeholder = "e.g. UL54 Mutation A123S is assigned as resistant to GCV, it is missed as conferring resistance to FOS"),
       #       textInput("db.ref", "Online reference ",placeholder = "e.g. https://www.biorxiv.org/content/link/to/our/paper"),
       #       p(), #separete intput from output
       #       checkboxInput("db.contact", "I wish to be contacted abut this change", FALSE),
       #       textInput("db.contact_email", "email address for contact"),
       #       actionButton("db.submit", "Submit", class = "btn-primary")
       #     ),
       #     mainPanel(
       #       h3("Submit any noted issues with the CMV resistance database i.e. any missed mutations in literature,
       #          additional drugs a gene mutation confers resitance to"),
       #       h4("here we will show any current issues, noted conflicts in the literature for an AAmut")
       #     ) #main
       #   ) #sbl
       # ), #tabpanel
       
       
       tabPanel("Database-metrics",
                strong("Gene - drug resistance heatmap, showing number of entries for each combination"),
                p(""),
                br(),
                selectInput("dbmetric_virus", label = "Virus choice", choices = c("HCMV", "HSV1", "HSV2", "VZV", "HHV6b"),multiple = T, selected = c("HCMV", "HSV1")),
                selectInput("dbmetric_drug", label = "Antiviral choice", choices = c("Ganciclovir", "Aciclovir", "Cidofovir", "Foscarnet", "Brincidofovir", "Letermovir", "Brivudine", "Penciclovir", "Tomeglovir", "Maribavir", "Amenamevir"), selected = c("Aciclovir", "Letermovir"),multiple = T),
                plotOutput("res.plot.dbheatmap")
                
       ), #tabpanel
       
       tabPanel("Help",
                p(""),
                br(),
                strong("This page desciribes usage and results from HerpesDRG"),
                br(),
                p(""),
                strong("Inputs"),
                br(),
                p("First select the virus references you want to use in processing from the drop-down list."),
                div("Then use the file upload box to upload a single sequence / variant file in any of the accepted formats. Currently:"),
                div(img(src = "1.png"), style = "text-align: center;"),
                HTML("<ul><li>FASTA (Whole genome, or fragments i.e. UL54 & UL89</li><li>Varaint Call Format >= Ver4.0 </li><li>Varscan tabluated format</li></ul>"),
                br(),
                strong("Clinical Overview"),
                div("The ouputs include a concise summary of the expected resistance phenotype against all therepeutic options in the database. Cutoff values are in line with the recommendations from 'The Third International Consensus Guidelines on the Managemenet of Cytomegalovirus in Solid-Organ Transplantation' and are presented as follows:"),
                HTML("<ul><li>High level (any mutation confers a fold change above 15)</li><li>Moderate level (maximal fold change returned is between 5-15)</li><li>Low resistance (maximal fold change returned is between 2-5)</li><li>No Resistance (maximal fold change returned was less than 2, or only anecdotal 'Polymorphism' was returned)</li><li>Resistant, magnitude unknown ( only anecdotal resistant data was returned)</li><li>NA (no mutations returned)</li></ul>"),
                div(img(src = "3.png"), style = "text-align: center;"),
                br(),
                p(""),
                br(),
                strong("Identified resistant mutants"),
                div("Shows all genetic variants that ave mathcing records in the resistance database identified in the upload. An explanation is below:", tags$br(),
                "green: Genetic information. Gene is Amino-acid change, freq is percentage of variant_reads / (all reads), ref_var_count is absolute values, ref_pos is relative nucleotide position to virus ncbi reference.", tags$br(),
                "red: Antiviral sensitivity change. Each header is an antiviral drug recorded within the database, a single mutation can map to many entries & an entry can map to many drugs. Values are the EC50 fold change determined relative to a wild type virus. 'Resistant' or 'Polymorphism' in entries are typically associated with anecdotal or less defined entries & users should consult the references", tags$br(), tags$br(),
                "blue: Reference information. reference is a hyperlink to the original reference where available. test_method is a brief summary of the method used to derive the entry as descriped in the references article, test_method_class in even more brief.", tags$br(), tags$br(),
                "co_gene and co_aa: some database entries contain tested viruses with co-occuring marker transferred mutations. Where flagged these columns are filled."),
                div(img(src = "4.png"), style = "text-align: center;")
                
                
       ),
       
       tabPanel("Terms of Use",
                p(""),
                br(),
                strong("Disclaimer"),
                div("This database has been created and is maintained by the Breuer Lab at UCL as a benefit to the research and education community. This tool is provided on an 'as is, best endeavours' basis only, and without guarantee to its accuracy or reliability.", tags$br(),
                  "We clearly state here that the tool is not validated for use in a diagnostic or clinical care context.", tags$br(),
                  "We clearly state here that any file uploaded should be devoid of patient identifiable features, such as filename.", tags$br(),
                  "When you upload a file to use this web service, you are accepting these terms."),
                br(),
                
                strong("Privacy Policy"),
                p("Every time you access the hosted webserver at a UCL web address  we reserve the right to register the following data for statistical and troubleshooting purposes."),
                HTML("<ul><li>Date and time of access</li><li>IP adress.</li></ul>"),
                p("The contents of uploaded files are only temporarily stored, and not read manually."),
                br(),
                br(),
                br(),
                br(),
                p("last updated: Nov 2023")
                
                
                
                
       ), #tabpanel
       tabPanel("About",
                p(""),
                br(),
                p("
                The effective use of antivirals is crucial to the management of herpesvirus infections in immuno-compromised hosts.
                
                Genotypic resistance testing is a widely adopted part of this strategy as resistance has been identified for all licensed drugs.
                
                HerpesDRG is the first open knowledgebase linking all published herpes mutations to their in-vitro impact on antiviral sensitivity.
                  
                This webserver is a simple tool to enable resistance genotyping from common virus sequence data."),
                br(),
                
                strong("Database"),
                p(
                  "The tabular database has been manually curated by reading published literature and extracting key information.",
                  
                  "Each entry is a mutation containing the link between virus, gene, AAchange & EC50 fold change relative to a control.
                
                The data comes from both marker transfer in-vitro experiments and characterised clinical isolates.
                
                We update the database manually circa every quarter with new data."), 
                textOutput("about_db_total"),
                textOutput("about_db_unique"),
                textOutput("about_db_most_recent"),
                br(),
                
                strong("Method:"),
                div("Reports are presented relative to the viral reference NCBI strain.
                
                If providiing a VCF, ensure the data are mapped to the correct reference strain, these are also defined in ./inst/ref.  
                
                If providing a fasta file, 'MAFFT --add --keeplength' aligns and data to the reference, therefore whole genome or genetic fragments (i.e. polymerase only) are equally accepted."),
                br(),
                
                strong("Open Source"),
                p(
                  "Anyone is welcome to propose changes to the knowledgebase or R package. 
                
                i.e. If you spot mistakes or omissions in the  either contact oscar.charles.18@ucl.ac.uk or make a pull request."
                ),
                
                p(HTML(paste0("Link to the ",a("herpesdrg-db github reposisitory", href = "https://github.com/ojcharles/herpesdrg-db")))) ,
                
                p(HTML(paste0("Link to the ",a("herpesdrg R package github repository", href ="https://github.com/ucl-pathgenomics/herpesdrg")))),
                br(),br(),
                
                
                strong("Contact & Referencing"),
                p("The main developer is Oscar Charles who can be contacted at this", HTML(paste0(a("Email.  ", href = "mailto:ocar.charles.18@ucl.ac.uk"))),
                  "To reference HerpesDRG please cite ", a("this article.", href = "https://www.biorxiv.org/content/10.1101/2020.05.15.097907v2")), 
                
                strong("Ackowledgements"),
                p("We would like to thank the MRC and Wellcome Trust as sponsors of this work."),
                fluidRow( 
                  column(4, img(width = "90%", src = "ucl.png")),
                  column(4, img(width = "90%", src = "mrc.png")),
                  column(4, img(width = "90%", src = "wt.jpg")),
                )
                
       ) #tabpanel
    ) #navlistpanel
  ) #fluidpage
) #shinyui
