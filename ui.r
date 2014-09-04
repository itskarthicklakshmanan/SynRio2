library(shiny)
library(graph)
library(Rgraphviz)
library(shinyIncubator)
library(shinyBS)
library(circlize)
library(seqLogo)
shinyUI(navbarPage("SynRio",inverse = TRUE,
   tabPanel("Home",icon = icon("fa fa-home"),
    sidebarLayout(
      sidebarPanel( tags$style(type='text/css', ".span4 { max-width: 400px; }"),includeScript("www/binder.js"),includeScript("www/jquery-ui.min.js"),includeCSS("www/jquery-ui.css"),
		 wellPanel(),bsAlert(inputId = "alert_anchor")
	  		 
	  
	  ),
      mainPanel(progressInit()
	  
	  
	  ))
    
  ),
  navbarMenu("G-View",icon = icon("fa fa-dashboard"),
     tabPanel("Circos view",icon = icon("fa fa-bullseye"),
    sidebarLayout(
      sidebarPanel( tags$style(type='text/css', ".span4 { max-width: 400px; }"),includeScript("www/binder.js"),includeScript("www/jquery-ui.min.js"),includeCSS("www/jquery-ui.css"),
		 wellPanel(checkboxInput("track", "Gene Panel",TRUE),conditionalPanel(
    condition = "input.track",checkboxInput("cog", "COG Track",TRUE),checkboxInput("anno", "Pathway Track",TRUE))

        
      ),fluidRow(column(6,selectInput("circo_input", "",
                    c(paste("S",1:72,sep="")))),column(4,actionButton("circo1", "",icon = icon("fa fa-eye"))))
	  		 
	  
	  ),
      mainPanel(tags$head(tags$style(type='text/css', "#myImage { min-height:700px;}")),
      plotOutput("circo", width = "1200px", height="1200px")
	  
	  
	  
	  
	  ))
    
  ),
  
  
  
  tabPanel("Browser",icon = icon("fa fa-exchange"),
    sidebarLayout(
      sidebarPanel(tags$style(type='text/css', ".span4 { max-width: 400px; }"),
	  wellPanel(fluidRow(column(2,actionButton("searchhelp", "",icon = icon("fa fa-question-circle"))),column(6,textInput("gene_ser", "","sll1514")),column(4,actionButton("geneserch", "Search",icon = icon("fa fa-search"))))),
		 wellPanel(fluidRow(
		 radioButtons("regsel", "", choices=c("Range Select" = "oner","Genome Segments" = "twor")),conditionalPanel(
    condition = "input.regsel == 'oner'",
		 column(5,
      numericInput("instart", "Start:", min = -200, max = 3000000, value = 1000, step = 3000)),column(5,
	  numericInput("inend", "End:",  min = 4000, max = 3600000, value = 10000, step = 3000))),conditionalPanel(
    condition = "input.regsel == 'twor'", column(2,actionButton("pushorfl", "",icon = icon("fa fa-arrow-circle-left"))),column(8,uiOutput("slidersmain")),column(2,actionButton("pushorfr", "",icon = icon("fa fa-arrow-circle-right"))))
        
      )),wellPanel(radioButtons("inCheckboxGroup", "Mode",choices=
                           c("Default" = "option0","COG Color" = "option1",
                             "Kasuza" = "option2"),selected="option0")),tableOutput('track'),conditionalPanel(
    condition = "input.anno_drop",wellPanel(
							 radioButtons("tracklist", "Tracks", choices=c("Tata-box" = "tata","Promoter" = "tp", "Non-coding RNA"="ncrna", "Anti-sense RNA"="asrna","Intergenic region"="intern",
							 "UTR"="utr"), selected = NULL)))
                            
        
      
	  		 
	  
	  ),
      mainPanel(tags$head(tags$style(type='text/css', "#myImage { min-height:500px;background-color:'white';}"),tags$style(type='text/css', "#myImage2 { max-height:300px; top:100px;}"),tags$style(type='text/css', "#genebox {color: #999999; max-width: 650px; align:center; position:absolute; top:240px; width:600px;}")),uiOutput("message"), 
	   fluidRow(column(4,fluidRow(column(2,checkboxInput("gridr", " Grid",TRUE)),column(2,checkboxInput("finetune", " Zoom",TRUE)),column(2,checkboxInput("nav_drop", "NavI",TRUE)),column(2,checkboxInput("anno_drop", " Track",TRUE)),column(2,checkboxInput("user_track", " UTS",FALSE)))),
	   column(4,selectInput("gselect", "","",multiple = TRUE)),column(3,actionButton("add_basket", "+",icon = icon("fa fa-shopping-cart")))),
      bsAlert(inputId = "alert_gene"),plotOutput("myImage", clickId ="kar", hoverId="kar1", width = "1200px", height="100px"),conditionalPanel(
    condition = "input.anno_drop",
	  	bsCollapse(multiple = TRUE, id = "collapse1", 
        bsCollapsePanel("GC Content", fluidRow(column(1,selectInput("gcnum", "", choices=c(50,100,200,300,400,500,600,700,800,900,1000)))),plotOutput("myImage3", width = "1200px", height="300px"), id="col1", value="test1"), 
		bsCollapsePanel("Annotation tracks", plotOutput("myImage2", width = "1200px", height="200px"), id="col2", value="test2"))),conditionalPanel(
    condition = "input.user_track",
	    bsCollapse(multiple = TRUE, id = "collapse2", 
        bsCollapsePanel("User Track1",plotOutput("myImage4", width = "1200px", height="400px"),id="col21", value="test21"), 
		bsCollapsePanel("User Track2",plotOutput("myImage5", width = "1200px", height="200px"),id="col22", value="test22"))),
	  
	  tags$style(type='text/css', "#genebox1 {color: #999999; max-width: 950px; align:center; position:absolute; top:160px; width:600px;}"),uiOutput("message1"),conditionalPanel(
	  condition = "input.finetune",
	  fluidRow(
wellPanel(
	 column(6,
       tags$div(title="Fine tune to view sequence",uiOutput("sliders"))),conditionalPanel(
	  condition = "input.nav_drop", absolutePanel(top= 180, right = 20, width = 210, height= 600,draggable = TRUE,
wellPanel(p("Navigator"),bsActionButton("twoleft", "<<"),bsActionButton("oneleft", "<"),bsActionButton("oneright", ">"),bsActionButton("tworight", ">>"),bsActionButton("pushleft", "<-"),bsActionButton("pushright", "->")))
        
          
      ),column(2,
        bsActionButton("moleft", "left"),bsActionButton("moright", "Right"),actionButton("locusserch", "",icon = icon("fa fa-refresh"))
		 ) 
      )  
	  )))
    )
  )),
      navbarMenu("Data",icon = icon("fa fa-th-list"),
	  
    tabPanel("Locus datatable",icon = icon("fa fa-table"),
	 sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Choose Genomic data:", 
                  choices = c("Gene", "Promoter", "TATA box","Antisense-RNA","Noncoding-RNA","5'UTR","HIPs","Intergenic sequence","Transcription unit")),
      radioButtons("filetype", "File type:",
                   choices = c("csv", "txt")),
      downloadButton('downloadData', 'Download')
    ),
    mainPanel(
      dataTableOutput('glist')
    )
  )
  ), tabPanel("Extract Sequence",icon = icon("fa fa-scissors"),
	 sidebarLayout(
    sidebarPanel(
      wellPanel(radioButtons("seqrange", "Sequence Type:",choices = c("Browser_Range", "Ranger"))),conditionalPanel(
    condition = "input.seqrange=='Browser_Range'",wellPanel(checkboxInput("gen_seq","Split By genes",FALSE),checkboxInput("gen_seqr","Reverse complement",FALSE))),
	conditionalPanel(
    condition = "input.seqrange=='Ranger'",wellPanel(fluidRow(column(4,numericInput("range_start","Start",1)),column(4,numericInput("range_end","End",20000)),
	column(3,actionButton("infuse1", "",icon = icon("fa fa-scissors"))))))
    ),
    mainPanel(tags$style(type='text/css', "#seq_bse { max-height:350px;min-width:1300px;overflow:auto;}"),tags$style(type='text/css', "#seq_bse1 { max-height:350px;min-width:1300px;overflow:auto;}"),conditionalPanel(
    condition = "input.seqrange=='Browser_Range'",
      verbatimTextOutput('seq_bse')),conditionalPanel(
    condition = "input.seqrange=='Ranger'",
      verbatimTextOutput('seq_bse1'))
    )
  )
  ),
    tabPanel("Upload dataset",icon = icon("fa fa-upload"),
       sidebarLayout(
      sidebarPanel(tags$style(type='text/css', ".span4 { max-width: 400px; }"),
		                    wellPanel(

	HTML("<p><span style=\"color:#336666;font-size:16px\">
Add annotation file</span></p>"),
	fileInput('file1', '',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),

    checkboxInput('header', 'Header', TRUE),
    radioButtons('sep', 'Separator',
                 c(Comma=',',
                   Semicolon=';',
                   Tab='\t',
				   Space=' '),
                 selected=",")),uiOutput("idfers"),  fluidRow(

    column(4,
      selectInput("input_type", "Select task",
        c( "Annotate","subset"
        
      )
    )),

    column(4, 
      # This outputs the dynamic UI component
      uiOutput("ui")
    )),

				 
				 
				 
				 
				 actionButton("merged", "Run"),hr(),conditionalPanel(
    condition = "input.merged",wellPanel(p("Select genes"),checkboxInput("gselectckb", "Select All",FALSE),fluidRow(column(8,selectInput("gselect1", "","",multiple = TRUE)),column(3,actionButton("add_basket1", "+",icon = icon("fa fa-shopping-cart"))))))),
        mainPanel(tags$style(type='text/css', "#contents { max-height:750px;overflow:auto;}"),  tabsetPanel(type = "tabs", 
        tabPanel("User data", tableOutput('contents')),
      

		
        tabPanel("Annotated data", dataTableOutput('contents2'))
      )
      )
  )
  )
  
  
  
  
  
  
  
  ),
    navbarMenu("Analysis",icon = icon("fa fa-cogs"),
	
	tabPanel("Cluster Gene Set",icon = icon("fa fa-sitemap"),
       sidebarLayout(
      sidebarPanel(tags$style(type='text/css', ".span4 { max-width: 400px; }"),
		 wellPanel(radioButtons("cluster", "Cluster Types",choices=
                           c("Gene Distance"= "maktr","Abiotic Stress Cluster" = "clust1","Go Cluster" = "clust2",
                             "KEGG Pathway Cluster" = "clust3"),selected="clust1"),hr(),fluidRow(column(10,checkboxInput("kgset","Select Genes from Gene Cart",FALSE)))
      )),
        mainPanel(conditionalPanel(
    condition = "input.cluster=='maktr'", 
      tabsetPanel(type = "tabs",
        tabPanel("Data Table",tags$style(type='text/css', "#disttab{ min-height: 600px; }"),
		wellPanel(fluidRow(column(3,selectInput("lend","Intergenic Distance",c(100,200,300,400,500,600,700,800,900))),column(2,selectInput("distt","Select Gene cluster","Select Row")),column(2,actionButton("basket1", "+",icon = icon("fa fa-shopping-cart"))),
		column(2,actionButton("ranger1", "+",icon = icon("fa fa-scissors"))),column(2,actionButton("rangerplot1", "+",icon = icon("fa fa-graph"))))),tableOutput("disttab")))
      ),conditionalPanel(
    condition = "input.cluster=='clust1'", 
      tabsetPanel(type = "tabs",
        tabPanel("Abiotic Stress profiling", fluidRow(column(7,plotOutput("stress_hmap",width = "800px", height="700px"))),absolutePanel(top= 150, right = 150, width = 400, height= 400,draggable = TRUE,plotOutput("stress_hmap2",width = "400px", height="500px"))), 
        tabPanel("Data Table", tableOutput("marrtab")))
      ),(conditionalPanel(
    condition = "input.cluster=='clust2'", 
      tabsetPanel(type = "tabs",
        tabPanel("Go cluster",column(7,plotOutput("go_pie",width = "1200px", height="600px"))),
		tabPanel("PCA",  column(7,plotOutput("go_map",width = "600px", height="600px")),wellPanel(verbatimTextOutput("go_sum")),wellPanel(verbatimTextOutput("go_load"))))
		)
      ), (conditionalPanel(
    condition = "input.cluster=='clust3'", 
      tabsetPanel(type = "tabs",
        tabPanel("Kegg cluster",""), 
        tabPanel("DATA", uiOutput("")))
      ))
  ))
  ),
  	tabPanel("Pattern Search",icon = icon("fa fa-align-right"),
       sidebarLayout(
      sidebarPanel(tags$style(type='text/css', ".span4 { max-width: 400px; min-height:600px;}"),tags$style(type='text/css', "#seqlogoplot { max-width: 600px; min-height:600px;}"),
	  radioButtons("patty", "", choices=c("User Pattern" = "patty_1","Pattern Dataset" = "patty_2")),hr(),conditionalPanel(
    condition = "input.patty=='patty_1'",
		wellPanel(textInput("pattern_ser1", "","GGCGATCGCC")
	
	)),conditionalPanel(
    condition = "input.patty=='patty_2'",
		wellPanel(fluidRow(column(12,selectInput("pattern_ser2", "",choices=c("TANNNT","GGCGATCGCC","TCATNNNNNNNTTTG","TGTGANNNNNNTCACA","RWWTRGGYNNYY","CAACNNNNNNGTTG","GATNNNNNNNNTAC",
		"CTANNNNNNNNNCTA","TTANNNNNNNTAA","GTTTGTTTT","TGACAGGCCGA"
)))))),fluidRow(column(3,checkboxInput("ambig","Exact",value = TRUE)),column(2,checkboxInput("comp_std","CPattern",value = TRUE))),
	
	fluidRow(column(5,uiOutput("sliders3"))),
		 hr(),radioButtons("patt", "", choices=c("Whole Genome" = "patt_1","Upstream Region" = "patt_2")),hr(),conditionalPanel(
    condition = "input.patt=='patt_1'",wellPanel(fluidRow(column(4,numericInput("patt_genlen_f","Start",1)),column(4,numericInput("patt_genlen_t","End",3600000)),
	column(3,actionButton("infuse", " ",icon = icon("fa fa-search")))),actionButton("pathook", "",icon = icon("fa fa-search"))))
	
	,
		 
		 conditionalPanel(
    condition = "input.patt=='patt_2'",wellPanel(fluidRow(column(5,selectInput("patt_promlen_t","Basepair length",c(50,100,200,300,400,500)))),checkboxInput("pgset","Use List from Gene Cart",FALSE))),
		 actionButton("patsearch", "Search Pattern",icon = icon("fa fa-search"))),
		 
        mainPanel(conditionalPanel(
    condition = "input.patt=='patt_1'",
      tabsetPanel(type = "tabs", 
        tabPanel("Method", plotOutput("patplot")), 
        tabPanel("Result", fluidRow(column(5,tableOutput("pattable")),column(6,plotOutput("seqlogoplot")),column(5,tableOutput("seqlogomat"))))
      )),conditionalPanel(
    condition = "input.patt=='patt_2'",
      tabsetPanel(type = "tabs", 
       tabPanel("Result", dataTableOutput("pattable1")),fluidRow(column(3,selectInput("filetype1", "",choices = c("csv", "txt"))),column(3,downloadButton('prompttData', 'Download'))),
	   tabPanel("Density Plot", plotOutput("patplotden"))
	   )
	   
      )
  ))
  ),
  	tabPanel("Network Views",icon = icon("fa fa-code-fork"),
       sidebarLayout(
      sidebarPanel(tags$style(type='text/css', ".span4 { max-width: 400px; }"), selectInput("inSelect3", "Dataset:",choices=c("genedata","userdata")),selectInput("inSelect2", "Select genes:",
                    multiple = TRUE,
                    ""),hr(),radioButtons("ppi", "", choices=c("Protein Network" = "ppi_1","KEGG Categories" = "kegg", "GO Categories" = "go")),conditionalPanel(
    condition = "input.ppi=='ppi_1'",
		 wellPanel(selectInput("int_source", "Interaction Source:", 
                  choices = c("string", "psimap", "ipfam","y2h")),
				  selectInput("int_type", "Network Type:", 
                  choices = c("neato", "circo", "twopi","dot","fdp")), actionButton("gobut", "layout network"))
      ),conditionalPanel(
    condition = "input.ppi=='kegg'",
		 wellPanel(selectInput("int_type1", "Network Type:", 
                  choices = c("neato", "circo", "twopi","dot","fdp")), actionButton("gobutk", "layout network"))
      ),conditionalPanel(
    condition = "input.ppi=='go'",
		 wellPanel(selectInput("int_type2", "Network Type:", 
                  choices = c("neato", "circo", "twopi","dot","fdp"),selected="fdp"),radioButtons("gotype", "", choices=c("Goid" = "goid","Go term" = "goterm")), actionButton("gobutc", "layout network"))
      )),
        mainPanel(tags$style(type='text/css', "#int {min-height:800px;max-width: 1200px; overflow:auto;}"),
		tags$style(type='text/css', "#intk {min-height:800px;max-width: 1200px; overflow:auto;}"),tags$style(type='text/css', "#intc {min-height:800px;max-width: 1200px; overflow:auto;}")
,conditionalPanel(
    condition = "input.ppi=='ppi_1'",
      tabsetPanel(type = "tabs",	  
        tabPanel("Interaction Network",bsAlert(inputId = "alert_int"),plotOutput("int")), 
        tabPanel("Result Table", dataTableOutput("inttable"))
      )),conditionalPanel(
    condition = "input.ppi=='kegg'",tabsetPanel(type = "tabs",	  
        tabPanel("Interaction Network",plotOutput("intk")) 
       
      )),conditionalPanel(
    condition = "input.ppi=='go'",tabsetPanel(type = "tabs",	  
        tabPanel("Gene ontology Network",plotOutput("intc")), 
        tabPanel("Result Table", dataTableOutput("ctable"))
      ))
  ))
  ),tags$div( id="dialog",tags$div( id="res"))
 
  
  
  ),tabPanel("Gene Cart",icon = icon("fa fa-shopping-cart"),
    sidebarLayout(
    sidebarPanel(actionButton("merge_glist", "Refresh gene set",icon = icon("fa fa-refresh")),hr(),helpText("subset data"),selectInput("inSelectg", "","",multiple = TRUE)

    ),
    mainPanel(tags$style(type='text/css', "#userdata { min-height:750px;overflow:auto;}"), tabsetPanel(type = "tabs",id ="inTabset",position ="right",
      tabPanel("Consolidated Gene cart", tableOutput('userdata3')),tabPanel("Gene set from Browser", tableOutput('userdata')),tabPanel("Gene set from Dist cluster",tableOutput('userdata1')),
	  tabPanel("Gene set from upload data",tableOutput('userdata4')))
    )
  )
  ),tabPanel("Help",icon = icon("fa fa-question-circle"),

  navlistPanel(
    "Select the topic",
    tabPanel("Introduction",
      h3("Introduction")
    ),
    tabPanel("G-View",
      h3("G-View")
    ),
	    tabPanel("Data",
      h3("Data")
    ),
	    tabPanel("Analysis",
      h3("Analysis")
    ),
    "-----",
    tabPanel("Color codes",
      h3("Color codes")
    ),
	"-----",
    tabPanel("Sample data",
      h3("Sample data")
    )
  )
  )
  
  
))

