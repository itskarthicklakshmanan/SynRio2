########################################################################
# SYNRIO2 : updated version
# Script written and maintained by:
# Karthick L
# PhD student
# Bioinformatics Division
# National Facility of Marine Cyanobacteria (NFMC)
# Bharathidasan University
# Trichy
# email: karthick.lakshman@gmail.com
########################################################################
library(shiny)
library(graph)
library(Rgraphviz)
library(shinyIncubator)
library(shinyBS)
library(circlize)
library(seqLogo)
shinyUI(navbarPage("SynRio",inverse = TRUE,fluid = TRUE,
   tabPanel("Home",icon = icon("fa fa-home"),
    sidebarLayout(
      sidebarPanel( tags$head(tags$style(style="text/css", "
             div.busy {
  position:absolute;
  top: 40%;
  left: 150%;
  margin-top: -100px;
  margin-left: -50px;
  display:none;
  background: rgba(230, 230, 230, .8);
  text-align: center;
  padding-top: 20px;
  padding-left: 30px;
  padding-bottom: 40px;
  padding-right: 30px;
  border-radius: 5px;
}

          ")),tags$style(style='text/css', ".span2{ max-width: 200px; }"),includeScript("www/binder.js"),includeScript("www/jquery-ui.min.js"),includeCSS("www/jquery-ui.css"),
		 wellPanel(h5("SynRio - View Synechocystis 6803 chromosome using R and Shiny"),hr(),p("Cretaed at:"),a(href="http://nfmc.res.in","National Facility for Marine Cyanobacteria",target="_blank"),hr(),
		 fluidRow(column(4,img(src="/img/bionfmc.png",height = 150,width=150)),column(4,img(src="/img/nfmc.png",height = 150,width=150)),column(4,img(src="/img/bdulogo.png",height = 150,width=150)))
		 )
	  		 
	  
	  ),
      mainPanel(h3("About SynRio"),hr(),p("Welcome to Syn-R-io, an interactive R application based on the shiny package for visual exploration of Synechocystis 6803 chromosome with simple data extraction options")
	  ,p("View Synechocystis genome as Circo or as a linear arrangement of genes. 
	  There are additional gene specific information retrival options and also some minimal data analysis features."),br(),p("As a quick start, please refer the 'Help' section or click on the G-View tab and select 'Browser'"),
	  fluidRow(column(8,wellPanel(h6("Quick start SynRio => (step 1) Click 'G-View' tab  |  (step 2) select 'Browser' option  |  (step 3) proceed to Data or Analysis tab"))))
	  
	  
	  ))
    
  ),
  navbarMenu("G-View",icon = icon("fa fa-dashboard"),
     tabPanel("Circos view",icon = icon("fa fa-bullseye"),
    sidebarLayout(
      sidebarPanel( tags$style(style='text/css', ".span4 { max-width: 400px; }"),includeScript("www/binder.js"),includeScript("www/jquery-ui.min.js"),includeCSS("www/jquery-ui.css"),
		 wellPanel(checkboxInput("track", "Gene Panel",TRUE),conditionalPanel(
    condition = "input.track",checkboxInput("cog", "COG Track",TRUE),checkboxInput("anno", "Pathway Track",TRUE))

        
      ),fluidRow(column(6,selectInput("circo_input", " Select Range",
                    c(paste("S",1:72,sep="")))),column(4,actionButton("circo1", "Go",icon = icon("fa fa-eye")))),
					textInput("segrang", "Selected genomic range in Base pairs:", ""),
					hr(),p("Legend"),img(src="/img/cir.png")
	  		 
	  
	  ),
      mainPanel(tags$head(tags$style(style='text/css', "#myImage { min-height:700px;}")),tags$div(class = "intro", p("Genome plot is loading..."), img(src="loader.gif")),
	  div(class = "busy", p("Calculation in progress.."), img(src="loader.gif")),
	  	
		

      plotOutput("circo", width = "1200px", height="1200px")
	  
	  
	  
	  
	  ))
    
  ),
  
  
  
  tabPanel("Browser",icon = icon("fa fa-exchange"),
    sidebarLayout(
      sidebarPanel(tags$style(style='text/css', ".span4 { min-width: 300px; }"),
	  wellPanel(fluidRow(column(2,actionButton("searchhelp", "",icon = icon("fa fa-question-circle"))),column(6,textInput("gene_ser", label = "")),column(4,actionButton("geneserch", "Search",icon = icon("fa fa-search"))))
	  ),
		 wellPanel(fluidRow(
		 radioButtons("regsel", "", choices=c("Range Select" = "oner","Genome Segments" = "twor")),conditionalPanel(
    condition = "input.regsel == 'oner'",
		 column(5,
      numericInput("instart", "Start:", min = -200, max = 3000000, value = 1000, step = 3000)),column(5,
	  numericInput("inend", "End:",  min = 4000, max = 3600000, value = 10000, step = 3000))),conditionalPanel(
    condition = "input.regsel == 'twor'", column(2,actionButton("pushorfl", "",icon = icon("fa fa-arrow-circle-left"))),column(8,uiOutput("slidersmain")),column(2,actionButton("pushorfr", "",icon = icon("fa fa-arrow-circle-right"))),
	textInput("segrang2", "View 50000 base pair range:", ""))
        
      )),wellPanel(radioButtons("inCheckboxGroup", "Mode",choices=
                           c("Default" = "option0","COG Color" = "option1",
                             "Kasuza" = "option2"),selected="option0")),tableOutput('track'),conditionalPanel(
    condition = "input.anno_drop",wellPanel(
							 radioButtons("tracklist", "Tracks", choices=c("Tata-box" = "tata","Promoter" = "tp", "Non-coding RNA"="ncrna", "Anti-sense RNA"="asrna","Intergenic region"="intern",
							 "UTR"="utr"), selected = NULL))),div(class = "busy", p("Calculation in progress.."), img(src="loader.gif"))
                            
        
      
	  		 
	  
	  ),
      mainPanel(tags$head(tags$style(style='text/css', "#myImage { min-height:500px;background-color:'white';}"),tags$style(style='text/css', "#myImage2 { max-height:300px; top:100px;}"),tags$style(style='text/css', "#genebox {color: #999999; max-width: 650px; align:center; position:absolute; top:240px; width:600px;}")),uiOutput("content"), 
	   fluidRow(column(4,fluidRow(column(2,checkboxInput("gridr", " Grid",TRUE)),column(2,checkboxInput("finetune", " Zoom",TRUE)),column(2,checkboxInput("nav_drop", "NavI",FALSE)),column(2,checkboxInput("anno_drop", " Track",TRUE)),column(2,checkboxInput("user_track", " UTS",FALSE)))),
	   column(5,absolutePanel(top= -4, right = 0, width = 500, height= 40,draggable = TRUE,selectInput("gselect", "","",multiple = TRUE))),column(3,actionButton("add_basket", "+",icon = icon("fa fa-shopping-cart")))),
      tags$div(class = "intro", p("Browser is loading..."), img(src="loader.gif")),plotOutput("myImage", click ="kar", hover="kar1", width = "1200px", height="100px"),conditionalPanel(
    condition = "input.anno_drop",
	  	bsCollapse(multiple = TRUE, id = "collapse1", 
        bsCollapsePanel("GC Content", fluidRow(column(1,selectInput("gcnum", "", choices=c(50,100,200,300,400,500,600,700,800,900,1000)))),plotOutput("myImage3", width = "1200px", height="300px"), id="col1", value="test1"), 
		bsCollapsePanel("Annotation tracks", plotOutput("myImage2", width = "1200px", height="200px"), id="col2", value="test2"))),conditionalPanel(
    condition = "input.user_track",
	    bsCollapse(multiple = TRUE, id = "collapse2", 
        bsCollapsePanel("User Track1",plotOutput("myImage4", width = "1200px", height="400px"),id="col21", value="test21"))),
	  
	  tags$style(style='text/css', "#genebox1 {color: #999999; max-width: 950px; align:center; position:absolute; top:160px; width:600px;}"),uiOutput("content1"),conditionalPanel(
	  condition = "input.finetune",
	  fluidRow(
wellPanel(
	 column(6,
       tags$div(title="Fine tune to view sequence",uiOutput("sliders"))),conditionalPanel(
	  condition = "input.nav_drop", absolutePanel(top= 180, right = 20, width = 210, height= 100,draggable = TRUE,
p("Navigator"),bsButton("twoleft", "<<"),bsButton("oneleft", "<"),bsButton("oneright", ">"),bsButton("tworight", ">>"),bsButton("pushleft", "<-"),bsButton("pushright", "->"))
        
          
      ),column(2,
        bsButton("moleft", "left"),bsButton("moright", "Right"),actionButton("locusserch", "",icon = icon("fa fa-refresh"))
		 ) 
      )  
	  )))
    )
  ),tabPanel("compare Plot",icon = icon("fa fa-sort-amount-asc"),
    sidebarLayout(
    sidebarPanel(selectInput("variable", "Select Cyano Genome",
                choices = c("Acaryochloris marina MBIC11017_chr","Acaryochloris marina MBIC11017_pREB1","Acaryochloris marina MBIC11017_pREB2","Acaryochloris marina MBIC11017_pREB3","Acaryochloris marina MBIC11017_pREB4","Acaryochloris marina MBIC11017_pREB5","Acaryochloris marina MBIC11017_pREB6","Acaryochloris marina MBIC11017_pREB7","Acaryochloris marina MBIC11017_pREB8","Acaryochloris marina MBIC11017_pREB9","Anabaena cylindrica PCC 7122_chr","Anabaena cylindrica PCC 7122_pANACY.01","Anabaena cylindrica PCC 7122_pANACY.03","Anabaena cylindrica PCC 7122_pANACY.04","Anabaena cylindrica PCC 7122_pANACY.06","Anabaena cylindrica PCC 7122_pANACY.05","Anabaena cylindrica PCC 7122_pANACY.02","Anabaena sp. 90_cir","Anabaena sp. 90_pANA01","Anabaena sp. 90_pANA02","Anabaena sp. 90_chr","Anabaena sp. 90_pANA03","Anabaena variabilis ATCC 29413_plasmid A","Anabaena variabilis ATCC 29413_plasmid B","Anabaena variabilis ATCC 29413_plasmid C","Anabaena variabilis ATCC 29413_chr","Anabaena variabilis ATCC 29413_ins element","Arthrospira platensis NIES-39_chr","Calothrix sp. PCC 6303_pCAL6303.01","Calothrix sp. PCC 6303_pCAL6303.02","Calothrix sp. PCC 6303_chr","Calothrix sp. PCC 6303_pCAL6303.03","Calothrix sp. PCC 7507_chr","Chamaesiphon minutus PCC 6605_chr","Chamaesiphon minutus PCC 6605_pCHA6605.02","Chamaesiphon minutus PCC 6605_pCHA6605.01","Chroococcidiopsis thermalis PCC 7203_chr","Chroococcidiopsis thermalis PCC 7203_pCHRO.02","Chroococcidiopsis thermalis PCC 7203_pCHRO.01","Crinalium epipsammum PCC 9333_pCRI9333.01","Crinalium epipsammum PCC 9333_pCRI9333.02","Crinalium epipsammum PCC 9333_pCRI9333.04","Crinalium epipsammum PCC 9333_pCRI9333.05","Crinalium epipsammum PCC 9333_pCRI9333.07","Crinalium epipsammum PCC 9333_chr","Crinalium epipsammum PCC 9333_pCRI9333.03","Crinalium epipsammum PCC 9333_pCRI9333.06","Crinalium epipsammum PCC 9333_pCRI9333.08","Cyanobacterium aponinum PCC 10605_chr","Cyanobacterium aponinum PCC 10605_pCYAN10605.01","Cyanobacterium stanieri PCC 7202_chr","Cyanobium gracile PCC 6307_chr","Cyanothece sp. ATCC 51142_plasmid A","Cyanothece sp. ATCC 51142_plasmid B","Cyanothece sp. ATCC 51142_plasmid C","Cyanothece sp. ATCC 51142_plasmid D","Cyanothece sp. ATCC 51142_cir","Cyanothece sp. ATCC 51142_linear","Cyanothece sp. PCC 7424_chr","Cyanothece sp. PCC 7424_pP742403","Cyanothece sp. PCC 7424_Pp742404","Cyanothece sp. PCC 7424_pP742405","Cyanothece sp. PCC 7424_pP742406","Cyanothece sp. PCC 7424_Pp742402","Cyanothece sp. PCC 7424_pP742401","Cyanothece sp. PCC 7425_pP742501","Cyanothece sp. PCC 7425_pP742503","Cyanothece sp. PCC 7425_chr","Cyanothece sp. PCC 7425_pP742502","Cyanothece sp. PCC 7822_chr","Cyanothece sp. PCC 7822_Cy782203","Cyanothece sp. PCC 7822_Cy782204","Cyanothece sp. PCC 7822_Cy782205","Cyanothece sp. PCC 7822_Cy782201","Cyanothece sp. PCC 7822_Cy782202","Cyanothece sp. PCC 7822_Cy782206","Cyanothece sp. PCC 8801_pP880101","Cyanothece sp. PCC 8801_pP880102","Cyanothece sp. PCC 8801_chr","Cyanothece sp. PCC 8801_pP880103","Cyanothece sp. PCC 8802_pP880201","Cyanothece sp. PCC 8802_chr","Cyanothece sp. PCC 8802_pP880202","Cyanothece sp. PCC 8802_pP880203","Cyanothece sp. PCC 8802_pP880204","Cylindrospermum stagnale PCC 7417_pCYLST.02","Cylindrospermum stagnale PCC 7417_chr","Cylindrospermum stagnale PCC 7417_pCYLST.03","Cylindrospermum stagnale PCC 7417_pCYLST.01","Dactylococcopsis salina PCC 8305_chr","Geitlerinema sp. PCC 7407_chr","Gloeobacter kilaueensis JS1_chr","Gloeobacter violaceus PCC 7421_chr","Gloeocapsa sp. PCC 7428_chr","Gloeocapsa sp. PCC 7428_pGLO7428.01","Gloeocapsa sp. PCC 7428_pGLO7428.04","Gloeocapsa sp. PCC 7428_pGLO7428.03","Gloeocapsa sp. PCC 7428_pGLO7428.02","Halothece sp. PCC 7418_chr","Leptolyngbya sp. PCC 7376_chr","Microcoleus sp. PCC 7113_chr","Microcoleus sp. PCC 7113_pMIC7113.01","Microcoleus sp. PCC 7113_pMIC7113.03","Microcoleus sp. PCC 7113_pMIC7113.05","Microcoleus sp. PCC 7113_pMIC7113.06","Microcoleus sp. PCC 7113_pMIC7113.08","Microcoleus sp. PCC 7113_pMIC7113.02","Microcoleus sp. PCC 7113_pMIC7113.04","Microcoleus sp. PCC 7113_pMIC7113.07","Microcystis aeruginosa NIES-843_chr","'Nostoc azollae' 0708_chr","'Nostoc azollae' 0708_pAzo01","'Nostoc azollae' 0708_pAzo02","Nostoc punctiforme PCC 73102_chr","Nostoc punctiforme PCC 73102_pNPUN05","Nostoc punctiforme PCC 73102_pNPUN03","Nostoc punctiforme PCC 73102_pNPUN01","Nostoc punctiforme PCC 73102_pNPUN02","Nostoc punctiforme PCC 73102_pNPUN04","Nostoc sp. PCC 7107_chr","Nostoc sp. PCC 7120_pCC7120beta","Nostoc sp. PCC 7120_pCC7120zeta","Nostoc sp. PCC 7120_pCC7120gamma","Nostoc sp. PCC 7120_pCC7120epsilon","Nostoc sp. PCC 7120_chr","Nostoc sp. PCC 7120_pCC7120delta","Nostoc sp. PCC 7120_pCC7120alpha","Nostoc sp. PCC 7524_pNOS7524.01","Nostoc sp. PCC 7524_chr","Nostoc sp. PCC 7524_pNOS7524.02","Oscillatoria acuminata PCC 6304_chr","Oscillatoria acuminata PCC 6304_pOSCIL6304.02","Oscillatoria acuminata PCC 6304_pOSCIL6304.01","Oscillatoria nigro-viridis PCC 7112_chr","Oscillatoria nigro-viridis PCC 7112_pOSC7112.02","Oscillatoria nigro-viridis PCC 7112_pOSC7112.03","Oscillatoria nigro-viridis PCC 7112_pOSC7112.05","Oscillatoria nigro-viridis PCC 7112_pOSC7112.01","Oscillatoria nigro-viridis PCC 7112_pOSC7112.04","Pleurocapsa sp. PCC 7327_chr","Prochlorococcus marinus str. AS9601_chr","Prochlorococcus marinus str. MIT 9211_chr","Prochlorococcus marinus str. MIT 9215_chr","Prochlorococcus marinus str. MIT 9301_chr","Prochlorococcus marinus str. MIT 9303_chr","Prochlorococcus marinus str. MIT 9312_chr","Prochlorococcus marinus str. MIT 9313_chr","Prochlorococcus marinus str. MIT 9515_chr","Prochlorococcus marinus str. NATL1A_chr","Prochlorococcus marinus subsp. marinus str. CCMP1375_chr","Prochlorococcus marinus subsp. pastoris str. CCMP1986_chr","Pseudanabaena sp. PCC 7367_pPSE7367.01","Pseudanabaena sp. PCC 7367_chr","Rivularia sp. PCC 7116_chr","Rivularia sp. PCC 7116_pRIV7116.01","Rivularia sp. PCC 7116_pRIV7116.02","Stanieria cyanosphaera PCC 7437_chr","Stanieria cyanosphaera PCC 7437_pSTA7437.02","Stanieria cyanosphaera PCC 7437_pSTA7437.03","Stanieria cyanosphaera PCC 7437_pSTA7437.01","Stanieria cyanosphaera PCC 7437_pSTA7437.05","Stanieria cyanosphaera PCC 7437_pSTA7437.04","Synechococcus elongatus PCC 6301_chr","Synechococcus elongatus PCC 7942_plasmid 1","Synechococcus elongatus PCC 7942_chr","Synechococcus sp. CC9311_chr","Synechococcus sp. CC9605_chr","Synechococcus sp. CC9902_chr","Synechococcus sp. JA-2-3B'a(2-13)_chr","Synechococcus sp. JA-3-3Ab_chr","Synechococcus sp. PCC 6312_chr","Synechococcus sp. PCC 6312_pSYN6312.01","Synechococcus sp. PCC 7002_pAQ7","Synechococcus sp. PCC 7002_chr","Synechococcus sp. PCC 7002_pAQ1","Synechococcus sp. PCC 7002_pAQ3","Synechococcus sp. PCC 7002_pAQ4","Synechococcus sp. PCC 7002_pAQ5","Synechococcus sp. PCC 7002_pAQ6","Synechococcus sp. PCC 7502_pSYN7502.01","Synechococcus sp. PCC 7502_pSYN7502.02","Synechococcus sp. PCC 7502_chr","Synechococcus sp. RCC307_chr","Synechococcus sp. WH 7803_chr","Synechococcus sp. WH 8109_chr","Synechocystis sp. PCC 6803_PsysM","Synechocystis sp. PCC 6803_pSYSA","Synechocystis sp. PCC 6803_PsysG","Synechocystis sp. PCC 6803_PsysX","Synechocystis sp. PCC 6803_chr","Synechocystis sp. PCC 6803_chr","Synechocystis sp. PCC 6803_pSYSA_M","Synechocystis sp. PCC 6803_pSYSG_M","Synechocystis sp. PCC 6803_pCA2.4_M","Synechocystis sp. PCC 6803_pCC5.2_M","Synechocystis sp. PCC 6803_pSYSM_M","Synechocystis sp. PCC 6803_pSYSX_M","Synechocystis sp. PCC 6803_pCB2.4_M","Synechocystis sp. PCC 6803 substr. GT-I_chr","Synechocystis sp. PCC 6803 substr. PCC-N_chr","Synechocystis sp. PCC 6803 substr. PCC-P_chr","Thermosynechococcus elongatus BP-1_chr","Thermosynechococcus sp. NK55_chr","Trichodesmium erythraeum IMS101_chr"), selected="30")            
				,numericInput("cone", "Genome viewer start:",100),numericInput("ctwo", "Genome viewer end:",30000),br(),actionButton("comp_but", "View Genomic region")
				

    ),
    mainPanel(tags$style(style='text/css', "#comp_plot { min-height:700px;overflow:auto;}"),tabsetPanel(style = "tabs", 
        tabPanel("Comparison Plot", div(class = "busy", p("Calculation in progress.."), img(src="loader.gif")),plotOutput("comp_plot")),
      

		
        tabPanel("Genome comparison data", dataTableOutput("comp_table"))
      )
    )
  )
  )),
      navbarMenu("Data",icon = icon("fa fa-th-list"),
	  
    tabPanel("Locus datatable",icon = icon("fa fa-table"),
	 sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Choose Genomic data:", 
                  choices = c("Gene", "Promoter", "TATA box","Antisense-RNA","Noncoding-RNA","5'UTR","HIPs","Intergenic sequence","Transcription unit")),
      radioButtons("filestyle", "File style:",
                   choices = c("csv", "txt")),
      downloadButton('downloadData', 'Download')
    ),
    mainPanel(tags$div(class = "intro", p("Locus Table is loading..."), img(src="loader.gif")),
      dataTableOutput('glist')
    )
  )
  ),tabPanel("Pathway data",icon = icon("fa fa-retweet"),
    sidebarLayout(
    sidebarPanel(selectInput("pathselect", "Select Pathway :",""),radioButtons("filestylep", "File style:",
                   choices = c("csv", "txt")),
      downloadButton('downloadDatap', 'Download')

    ),
    mainPanel(tags$style(style='text/css', "#userdata { min-height:750px;overflow:auto;}"),tags$div(class = "intro", p("Pathway table is loading..."), img(src="loader.gif")), dataTableOutput("path_table")
    )
  )
  ), tabPanel("Extract Sequence",icon = icon("fa fa-scissors"),
	 sidebarLayout(
    sidebarPanel(
      wellPanel(radioButtons("seqrange", "Sequence style:",choices = c("Browser_Range", "Sequence Range","Promoter Range"))),conditionalPanel(
    condition = "input.seqrange=='Browser_Range'",wellPanel(checkboxInput("gen_seq","Split By genes",FALSE),checkboxInput("gen_seqr","Reverse complement",FALSE))),
	conditionalPanel(
    condition = "input.seqrange=='Sequence Range'",wellPanel(fluidRow(column(4,numericInput("range_start","Start",1)),column(4,numericInput("range_end","End",20000)),
	column(3,actionButton("infuse1", "",icon = icon("fa fa-scissors")))))),
	conditionalPanel(
    condition = "input.seqrange=='Promoter Range'",wellPanel(p("Sequence Upstream in BP"),fluidRow(column(4,numericInput("promrange_start","From",1)),column(4,numericInput("promrange_end","To",200))),
	checkboxInput("ckbinprom","Genes from Cart",FALSE),actionButton("prominfuse1", "Extract",icon = icon("fa fa-scissors"))),conditionalPanel(
    condition = "input.prominfuse1",
	fluidRow(column(3,radioButtons("filestyle2", "",choices = c("csv", "txt"))),column(5,downloadButton('promextractdown', 'Download')))))
    ),
    mainPanel(tags$style(style='text/css', "#seq_bse { max-height:350px;min-width:1300px;overflow:auto;}"),tags$style(style='text/css', "#seq_bse1 { max-height:350px;min-width:1300px;overflow:auto;}"),conditionalPanel(
    condition = "input.seqrange=='Browser_Range'",
      verbatimTextOutput('seq_bse')),conditionalPanel(
    condition = "input.seqrange=='Sequence Range'",
      verbatimTextOutput('seq_bse1')),conditionalPanel(
    condition = "input.seqrange=='Promoter Range'",
      dataTableOutput('promtable1'))
    )
  )
  ),
    tabPanel("Upload dataset",icon = icon("fa fa-upload"),
       sidebarLayout(
      sidebarPanel(tags$style(style='text/css', ".span4 { max-width: 400px; }"),
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
      selectInput("input_style", "Select task",
        c( "Annotate","subset"
        
      )
    )),

    column(4, 
      # This outputs the dynamic UI component
      uiOutput("ui")
    )),

				 
				 
				 
				 
				 actionButton("merged", "Run"),hr(),conditionalPanel(
    condition = "input.merged",wellPanel(p("Select genes"),fluidRow(column(8,selectInput("gselect1", "","",multiple = TRUE)),column(3,actionButton("add_basket1", "+",icon = icon("fa fa-shopping-cart"))))))),
        mainPanel(tags$style(style='text/css', "#contents { max-height:750px;overflow:auto;}"),  tabsetPanel(style = "tabs", 
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
      sidebarPanel(tags$style(style='text/css', ".span4 { max-width: 400px; }"),tags$style(style='text/css', "#go_maptab { max-height: 200px; max-width:1200px;overflow:auto;}"),tags$style(style='text/css', "#kegg_maptab { max-height: 200px; max-width:1200px;overflow:auto;}"),
		 wellPanel(radioButtons("cluster", "Cluster styles",choices=
                           c("Gene Distance"= "maktr","Abiotic Stress Cluster" = "clust1","Go Cluster" = "clust2",
                             "KEGG Pathway Cluster" = "clust3"),selected="clust1"),hr(),fluidRow(column(10,checkboxInput("kgset","Select Genes from Gene Cart",FALSE)))
      )),
        mainPanel(conditionalPanel(
    condition = "input.cluster=='maktr'", 
      tabsetPanel(style = "tabs",
        tabPanel("Data Table",tags$style(style='text/css', "#disttab{ min-height: 600px; }"), p("Gene distance analysis is functional if genomic range is > 10000 bp"),
		wellPanel(fluidRow(column(3,selectInput("lend","Intergenic Distance",c(100,200,300,400,500,600,700,800,900))),column(2,selectInput("distt","Select Gene cluster","Select Row")),column(2,actionButton("basket1", "+  Add to gene cart",icon = icon("fa fa-shopping-cart"))),
		column(4,actionButton("ranger1", "+  Add to seq extracter",icon = icon("fa fa-scissors"))))),tableOutput("disttab")))
      ),conditionalPanel(
    condition = "input.cluster=='clust1'", 
      tabsetPanel(style = "tabs",
        tabPanel("Abiotic Stress profiling", fluidRow(column(7,tags$div(class = "intro1", p("Please select Browser tab First...")),plotOutput("stress_hmap",width = "800px", height="700px"))),absolutePanel(top= 150, right = 150, width = 400, height= 400,draggable = TRUE,plotOutput("stress_hmap2",width = "400px", height="500px"))), 
        tabPanel("Data Table", tableOutput("marrtab")))
      ),(conditionalPanel(
    condition = "input.cluster=='clust2'", 
      tabsetPanel(style = "tabs",
	    tabPanel("Go cluster",  fluidRow(column(8,plotOutput("go_mapfirst",width = "600px", height="600px")),column(2,checkboxInput("fitchk","Group Cluster",FALSE)),conditionalPanel(
    condition = "input.fitchk", column(4,numericInput("fitnumb","",3)))),
		br(),wellPanel(tableOutput("go_maptab"))),
        tabPanel("Break cluster",radioButtons("gostylep", "Go style",choices=c("Biological Process"= "B","Molecular Function"= "M","Cellular Componenet"= "C"),selected="B"),column(7,plotOutput("go_pie",width = "1200px", height="600px"))),
		tabPanel("PCA Plots",  column(7,plotOutput("go_map1",width = "600px", height="600px")), column(7,plotOutput("go_map",width = "600px", height="600px"))),
		tabPanel("PCA Table",  wellPanel(verbatimTextOutput("go_sum")),wellPanel(verbatimTextOutput("go_fit")),wellPanel(verbatimTextOutput("go_load"))))
		)
      ), (conditionalPanel(
    condition = "input.cluster=='clust3'", 
      tabsetPanel(style = "tabs",
        tabPanel("Kegg cluster",fluidRow(column(6,plotOutput("kegg_map",width = "600px", height="600px")),column(5,plotOutput("kegg_mappie",width = "600px", height="600px"))),br(),wellPanel(tableOutput("kegg_maptab"))))
      ))
  ))
  ),
  	tabPanel("Pattern Search",icon = icon("fa fa-align-right"),
       sidebarLayout(
      sidebarPanel(tags$style(style='text/css', ".span4 { max-width: 400px; min-height:600px;}"),tags$style(style='text/css', "#seqlogoplot { max-width: 600px; min-height:600px;}"),
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
	column(3,actionButton("infuse", " ",icon = icon("fa fa-search")))),actionButton("pathook", "Add to User Track",icon = icon("fa fa-exchange")))),
		 
		 conditionalPanel(
    condition = "input.patt=='patt_2'",wellPanel(fluidRow(column(5,selectInput("patt_promlen_t","Basepair length",c(50,100,200,300,400,500)))),checkboxInput("pgset","Use List from Gene Cart",FALSE))),
		 actionButton("patsearch", "Search Pattern",icon = icon("fa fa-search"))),
		 
        mainPanel(conditionalPanel(
    condition = "input.patt=='patt_1'",
      tabsetPanel(style = "tabs", 
        tabPanel("Method", plotOutput("patplot")), 
        tabPanel("Result", fluidRow(column(5,tableOutput("pattable")),column(6,plotOutput("seqlogoplot")),column(5,tableOutput("seqlogomat"))))
      )),conditionalPanel(
    condition = "input.patt=='patt_2'",
      tabsetPanel(style = "tabs", 
       tabPanel("Result", dataTableOutput("pattable1")),conditionalPanel(
    condition = "input.patsearch",fluidRow(column(3,selectInput("filestyle1", "",choices = c("csv", "txt"))),column(3,downloadButton('prompttData', 'Download')))),
	   tabPanel("Density Plot", plotOutput("patplotden"))
	   )
	   
      )
  ))
  ),
  	tabPanel("Network Views",icon = icon("fa fa-code-fork"),
       sidebarLayout(
      sidebarPanel(tags$style(style='text/css', ".span4 { max-width: 400px; }"), selectInput("inSelect3", "Dataset:",choices=c("genedata","userdata")),selectInput("inSelect2", "Select genes:",
                    multiple = TRUE,
                    ""),hr(),radioButtons("ppi", "", choices=c("Protein Network" = "ppi_1","KEGG Categories" = "kegg", "GO Categories" = "go")),conditionalPanel(
    condition = "input.ppi=='ppi_1'",
		 wellPanel(selectInput("int_source", "Interaction Source:", 
                  choices = c("string", "psimap", "ipfam","y2h")),
				  selectInput("int_type", "Network style:", 
                  choices = c("neato", "circo", "twopi","dot","fdp")),radioButtons("ppi_col", "", choices=c("COG color" = "ppic_1","KEGG Color" = "ppic_2")), actionButton("gobut", "layout network"))
      ),conditionalPanel(
    condition = "input.ppi=='kegg'",
		 wellPanel(selectInput("int_type1", "Network style:", 
                  choices = c("neato", "circo", "twopi","dot","fdp")), actionButton("gobutk", "layout network"))
      ),conditionalPanel(
    condition = "input.ppi=='go'",
		 wellPanel(selectInput("int_type2", "Network style:", 
                  choices = c("neato", "circo", "twopi","dot","fdp"),selected="fdp"),radioButtons("gostyle", "", choices=c("Goid" = "goid","Go term" = "goterm")),
  actionButton("gobutc", "layout network"))
      )),
        mainPanel(tags$style(style='text/css', "#int {min-height:800px;max-width: 1200px; overflow:auto;}"),
		tags$style(style='text/css', "#intk {min-height:800px;max-width: 1200px; overflow:auto;}"),tags$style(style='text/css', "#intc {min-height:800px;max-width: 1200px; overflow:auto;}")
,conditionalPanel(
    condition = "input.ppi=='ppi_1'",
      tabsetPanel(style = "tabs",	  
        tabPanel("Interaction Network",tags$div(class = "intro1", p("Please select Browser tab First...")), plotOutput("int")), 
        tabPanel("Result Table", dataTableOutput("inttable"))
      )),conditionalPanel(
    condition = "input.ppi=='kegg'",tabsetPanel(style = "tabs",	  
        tabPanel("Interaction Network",plotOutput("intk")) 
       
      )),conditionalPanel(
    condition = "input.ppi=='go'",tabsetPanel(style = "tabs",	  
        tabPanel("Gene ontology Network",plotOutput("intc")), 
        tabPanel("Result Table", dataTableOutput("ctable"))
      ))
  ))
  ),tags$div( id="dialog",tags$div( id="res"))
 
  
  
  ),tabPanel("Gene Cart",icon = icon("fa fa-shopping-cart"),
    sidebarLayout(
    sidebarPanel(actionButton("merge_glist", "Refresh gene set",icon = icon("fa fa-refresh")),hr(),helpText("subset data"),selectInput("inSelectg", "","",multiple = TRUE)

    ),
    mainPanel(tags$style(style='text/css', "#userdata { min-height:750px;overflow:auto;}"), tabsetPanel(style = "tabs",id ="inTabset",position ="right",
      tabPanel("Consolidated Gene cart", tableOutput('userdata3')),tabPanel("Gene set from Browser", tableOutput('userdata')),tabPanel("Gene set from Dist cluster",tableOutput('userdata4')),
	  tabPanel("Gene set from upload data",tableOutput('userdata1')))
    )
  )
  ),tabPanel("Help",icon = icon("fa fa-question-circle"),

  navlistPanel(
    "",
    tabPanel("Introduction",
      h3("Introduction"),hr()
	  ,p("Syn-R-io could be used for visualization of the Synechocystis genome as Circo or as linear gene viewer. 
	  There are additional gene specific information retrival options and also some minimal data analysis features.
	  
	  
"),br(),img(src="/img/wf.png")
    ),
    tabPanel("G-View",
      h3("G-View"),hr(),h4("Circos View"),p("Synechocystis Genomic range with 50,000 BP size is viewable as Circos plot. Use the Dropdown menu to choose the range from S1-S72 and click the 'View' button"),br()
	  ,
	  img(src="/img/cir_1.png"),br(),br(),p("The circos plot gives a overview of genomic region highlighting plus and minus strands"),img(src="/img/cir.png"),img(src="/img/circo.png"),br(),br(),
	  p("The inner concentric circles represent COG and KEGG annotations with color codes as specified in the 'color codes' section(left side bar)"),hr(),h4("Browser View"),p(""),img(src="/img/cir_4.png"),
	  p("1) View Gene locus using Synechocystis gene id"),
	  p(" 2) Choose 'Range Select' option to specify Genome region(From- To) or choose 'Genome segments' option (s1-s72) to view a 50,000bp range"),
	  p("3) Choose gene color based on COG or Kasuza annotation"),
	  p("4) Select Annotation tract and view the tracks in 9"),
	  p("5) Select check box to view zoom panel(10), navigation panel(11) toggle annotation track(8,9) and view user track"),
	  p("6) select genes from the dropdown, based on the gene browser view click the button to add the selected genes to gene cart"),
	  p("7) Genome plot"),
	  p("8) plot GC content"),
	  p("9) Annotation tracks, changes as per the selection in 4"),
	  p("10) Range selecter to zoom in to genome plot panel"),
	  p("11) Navigate to left and right of the genome locus")
	  
    ),
	    tabPanel("Data",
      h3("Data"),hr(),h4("Locus Datatable"),img(src="/img/data.png"),p("1) Choose data style to download the options include, Gene, Promoter, Tatabox, Intergenomic region, Transcription unit, UTR, etc"),p("2) Download data as .csv ot .txt"),p("3) Data table"),hr(),h4("Extract Sequence"),img(src="/img/extr.png"),
	  p("1) Select genome range to extract sequence"),p("2) Specify genome range or upstream region in bp"),p("3) Extracted sequence"),hr(),h4("Upload Geneset"),img(src="/img/up.png"),
	  p("1) Browse data file containing a list of gene ids"),p("2) Specify data file structure"),p("3) Choose header column from the uploaded file"),p("4) Run to pick genes and add to gene cart")
    ),
	    tabPanel("Analysis",
      h3("Analysis"),hr(),h4("Cluster Geneset"),h5("Important: Please select genome browser range before doing the cluster analysis!"),img(src="/img/dist.png"),p("1) Select Gene Distance option to pick gene set based on intergenic distance between genes"),p("2) Specify Gene distance in BP"),p("3) Data output based on gene distance number from 2"),p("4) Choose a gene group"),p("5) Add to gene cart"),img(src="/img/array.png"),p("1) Choose Abiotic stress cluster to cluster gene set from genome browser range"),p("2) Heatmap plot based on stress profile"),img(src="/img/goc.png"),p("1) Choose Go cluster to group genes based on Go terms"),p("2) Distance plot output based on Go"),p("3) Break as number of groups"),img(src="/img/goc1.png"),p("Go groups based on Biological process"),img(src="/img/keggc.png"),p("1) Choose to group genes based on KEGG pathway"),p("2) Hirarchical plot based on KEGG pathways"),p("3) Pathway group"),hr(),h4("Pattern Search"),img(src="/img/patt.png"),
	  p("1) Choose between 'user pattern'  or 'Pattern dataset' to select a genome pattern"),
	  p("2) Confirm or modify pattern"),
	  p("3) Select Exact if pattern is exact or deselect if the pattern has ambiguous letters. Check CPattern to include complement strand in the search"),
	  p("4) Search patteren in whole genome or only in upstream region"),
	  p("5) If upstream region is choosen pease select the number of base pairs"),
	  p("7) Density plot of the pattern"),
	  p("8) Download data"),hr(),h4("Network Views"),img(src="/img/net.png"),
	  p("1) Choose between 'Genesets'  or 'user sets' to select a data source"),
	  p("2) Pick specific gene list"),
	  p("3) Choose protein network option to group genes based on PPI or KEGG categories to view KEGG based network or else select Go based associations"),
	  p("5) Choose network style"), p("6) Color genes based onCOG or KEGG pathways"),
	  p("7) Gene networs plot")
    ),
    "-----",
    tabPanel("Color codes",
      h3("Color codes"), hr(),h5("KEGG Color"),img(src="/img/keggclr.png"),br(),
	  hr(),h5("COG Color"),img(src="/img/cogclr.png")
    ),
	"-----",
    tabPanel("Sample data",
      h3("Sample data"), hr(),h4("Highly expressed abiotic stress genes in Synechosystis genome"),br(),
	  wellPanel(fluidRow(column(4, radioButtons("stressdataset", "Choose a stress style:", choices = c("COLD", "HEAT", "LIGHT", "HYPEROSMOTIC","NACL","H202","UV")),actionButton("stressbut", "Select list")),
	  column(5,radioButtons("filestyle3", "File style:",
                   choices = c("csv", "txt"))),
      column(5,downloadButton('downData', 'Download')))
    ),
	  dataTableOutput("stresstab")
    )
  )
  )
  
  
))

