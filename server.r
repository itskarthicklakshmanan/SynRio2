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
options(shiny.maxRequestSize=10*1024^2)
library(graph)
library(Rgraphviz)
library(shinyBS)
library(circlize)
require(Biostrings, quietly = TRUE)
library(shinyIncubator)
library(seqLogo)
source("helper.R")
source('all_sessions.R', local=TRUE)

shinyServer(function(input, output, session) {


observe({ 
	if (input$geneserch==0){ return(NULL)}
	input$geneserch
	isolate({
	gone<-subset(kar,kar$id==input$gene_ser)
	updateNumericInput(session=session, "instart", value=as.numeric(gone$start)-4000)
    updateNumericInput(session=session, "inend", value=as.numeric(gone$end)+4000)
	}) 
  }) 
  


  observe({ 
	isolate({	
    updateSelectInput(session=session, "sel_one", choices=paste(syn[,2],syn[,3],syn[,1],sep=" - "))
	}) 
  }) 
  
  
  
  
output$sliders <- renderUI({

    r1<-as.integer(input$instart) 
	r2<-as.integer(input$inend)
	    slide<-sliderInput("inslider", "",
                 min =r1 , max = r2, value = c(r1, r2), step = 50)
  
  })
  
  
  
  output$slidersmain <- renderUI({
	    slide<-sliderInput("inslider2", "",
                 min =1 , max = 71, value = 2, step = 1,ticks = TRUE)
  
  })

  
  
    observe({ 
	if (input$pushorfl==0){ return(NULL)}
	input$pushorf1
	isolate({
	updateSliderInput(session=session, "inslider2", value=input$inslider2-1)

	}) 
  })
      observe({ 
	if (input$pushorfr==0){ return(NULL)}
	input$pushorfr
	isolate({
	updateSliderInput(session=session,"inslider2", value=input$inslider2+1)

	}) 
  })

  observe({ 
	if (input$regsel==0){ return(NULL)}
	if(input$regsel=="twor")
	{
	updateTextInput(session=session, "instart", value=seg[input$inslider2,2])
    updateTextInput(session=session, "inend", value=seg[input$inslider2,3])
	
	} else{
	updateTextInput(session=session,"instart", value=1000)
    updateTextInput(session=session, "inend", value=10000)
	}
  }) 
  
  observe({ 
	if (input$locusserch==0){ return(NULL)}
	input$locusserch
	isolate({
	updateNumericInput(session=session, "instart", value=1000)
    updateNumericInput(session=session, "inend", value=10000)
	}) 
  }) 
#-----------------------------------------------------navigator--------------------------------------------------------------   
  
observe({ 
	if (input$twoleft==0){ return(NULL)}
	input$twoleft
	isolate({
	updateNumericInput(session=session, "instart", value=as.numeric(input$instart)-8000)
    updateNumericInput(session=session, "inend", value=as.numeric(input$inend)-8000)
	}) 
  }) 
 observe({ 
	if (input$oneleft==0){ return(NULL)}
	input$oneleft
	isolate({
	updateNumericInput(session=session,"instart", value=as.numeric(input$instart)-4000)
    updateNumericInput(session=session, "inend", value=as.numeric(input$inend)-4000)
	}) 
  }) 
observe({ 
	if (input$tworight==0){ return(NULL)}
	input$tworight
	isolate({
	updateNumericInput(session=session, "instart", value=as.numeric(input$instart)+8000)
    updateNumericInput(session=session, "inend", value=as.numeric(input$inend)+8000)
	}) 
  })  
 observe({ 
	if (input$oneright==0){ return(NULL)}
	input$oneright
	isolate({
	updateNumericInput(session=session, "instart", value=as.numeric(input$instart)+4000)
    updateNumericInput(session=session, "inend", value=as.numeric(input$inend)+4000)
	}) 
  }) 
  observe({ 
	if (input$pushright==0){ return(NULL)}
	input$pushright
	isolate({
    updateNumericInput(session=session, "inend", value=as.numeric(input$inend)+4000)
	}) 
  }) 
    observe({ 
	if (input$pushleft==0){ return(NULL)}
	input$pushleft
	isolate({
    updateNumericInput(session=session, "inend", value=as.numeric(input$inend)-4000)
	}) 
  }) 
  
  observe({ 
	if (input$instart<9000 | input$inend<16000){
    updateButton(session, "twoleft",  disabled = TRUE)}else{updateButton(session, "twoleft", disabled = FALSE)}
  }) 

    observe({ 
	if (input$instart<5000 |input$inend<6000){
    updateButton(session, "oneleft",  disabled = TRUE)}else{updateButton(session, "oneleft", disabled = FALSE)}

  }) 

    observe({ 
	if (input$instart<5000){
    updateButton(session, "pushleft",  disabled = TRUE)}else{updateButton(session, "pushleft", disabled = FALSE)}
  }) 

  
  observe({ 
    if (input$moleft==0){ return(NULL)}
	input$moleft
		isolate({
	 updateSliderInput(session=session, "inslider", value=as.numeric(input$inslider)-50)

  })  })
  
    observe({ 
    if (input$moright==0){ return(NULL)}
	input$moright
	isolate({

	 updateSliderInput(session=session, "inslider", value=as.numeric(input$inslider)+50)

  })  })
  

  
#-----------------------------------------------------navigator--------------------------------------------------------------  
#-----------------------------------------------------plot start-------------------------------------------------------------- 

plotdata<-reactive({

if(input$inCheckboxGroup=="option1"){
outputclr<-cogc}
else if(input$inCheckboxGroup=="option2"){
outputclr<-kasc} else{
outputclr<-cogc}
})

observe({

segd<-subset(seg,seg$V1==input$circo_input)
res<-paste(segd[,2],"-", segd[,3])
updateTextInput(session=session, "segrang", value=res)

})

observe({

segd<-subset(seg,seg$V1==paste("S",input$inslider2,sep=""))
res<-paste(segd[,2],"-", segd[,3])
updateTextInput(session=session, "segrang2", value=res)

})
	
output$circo <- renderPlot({
plotd<-plotdata()

	if (as.character(input$track)=="TRUE"){data<-get(paste("plotd"))}else{data<-NULL}
	if (as.character(input$cog)=="TRUE"){cog<-get(paste("cogc"))}else{cog<-NULL}
	if (as.character(input$anno)=="TRUE"){kasc<-get(paste("kasc"))}else{kasc<-NULL}
if (is.null(input$circo1)){ return(NULL)}
input$circo1
isolate({
	if (is.null(input$circo_input)){cin<-NULL}else{cin<-paste(input$circo_input)}

mplot<-circo_plot(data,seg,cin,kasc,cog,trunit)

return(mplot)

})

})





sel_gene<-reactive({
gres<-subset(kar,kar$id==input$gene_ser)
return(as.character(gres$name))
})






marr<-reactive({
	if(input$kgset==TRUE){
		kar<-plotdata()
		dfg<-data.frame(id=input$inSelectg)
		dfg<-data.frame(id=dfg[,1])
		kar1<-merge(kar,dfg,by="id")
	}else{
		kar<-plotdata()
		kar1<-subset(kar,kar$start>input$inslider[1] & kar$end < input$inslider[2])}
		mdat<-merge(mdata,kar1,by="id",sort=F)
		return(mdat)

})


distabb<-reactive({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(input$lend)){ return(NULL)}
		kar<-plotdata()
		kar1<-subset(kar,kar$start>input$instart & kar$end < input$inend)
		kar1<-kar1[order(kar1[,3]), ]
		if(input$inend - input$instart>10000){
		genecluster(kar1,input$lend)}
	})

output$disttab<-renderTable({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(input$lend)){ return(NULL)}
		distabb()

	})

observe({
	updateSelectInput(session=session, "distt", choices=rownames(distabb()))

	})



output$stress_hmap<-renderPlot({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		mdat<-marr()
		mdat1<- mdat[,2:9]
		mdat1<-as.matrix(mdat1,stringsAsFactors = TRUE)
		dimnames(mdat1)<-list(mdat[,1],colnames(mdat1[,1:8]))
		rc <- rainbow(nrow(mdat1), start = 0, end = .3)
		cc <- rainbow(ncol(mdat1), start = 0, end = .3)
		heatmap(mdat1, col = c("red","yellow","green"), scale = "column",
		RowSideColors = mdat[,16],
		xlab ="", ylab =  "",
		main = "")
})

output$stress_hmap2<-renderPlot({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
	mdat<-marr()
	mdat1<- mdat[,c(15,16)]
	mdat2<- table(mdat1[,1])
	df<-data.frame(cat=names(mdat2))
	colnames(df)<-"cat"
	col<-merge(df,mdat1,by="cat")
	pie(mdat2,labels=names(mdat2),col=unique(col[,2]), radius = 0.5,cex = 0.8)
})

output$marrtab1<- renderTable({
    mdat<-marr()
	mdat1<- unique(mdat[,c(15,16)])
	#colnames(mdat1)<-c("Genes","High Light","BvsG Light","Osmotic","UV","Temperature","H202","Cold","Salt","strand","Gene Description")
  })


output$marrtab<- renderTable({
    mdat<-marr()
	mdat1<- mdat[,c(1:9,13,14)]
	#colnames(mdat1)<-c("Genes","High Light","BvsG Light","Osmotic","UV","Temperature","H202","Cold","Salt","strand","Gene Description")
  })

 
keggs<-reactive({

	if(input$kgset==TRUE){
		kar<-plotdata()
		dfg<-data.frame(id=input$inSelectg)
		dfg<-data.frame(id=dfg[,1])
		kar1<-merge(kar,dfg,by="id")
	}else{
		kar<-plotdata()
		kar1<-subset(kar,kar$start>input$inslider[1] & kar$end < input$inslider[2])}
		mdat<-merge(path,kar1,by="id",sort=F)
		return(mdat)
	})




output$kegg_map<-renderPlot({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		gdat<-keggs()
		gdat1<- gdat[,c(1,2)]
		d<-dist(table(gdat1)) 
		fit <- hclust(d) 
		plot(fit, main="Go Cluster", xlab=NULL)
	})

output$kegg_maptab<-renderTable({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		gdat<-keggs()
		gdat1<- gdat[,c(1,3)]
		t(table(gdat1))
	})

output$kegg_mappie<-renderPlot({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		mdat<-keggs()
		mdat1<- mdat[,c(1,4)]
		mdat2<- table(mdat1[,2])
		pie(mdat2,labels=names(mdat2), radius = 0.5,cex = 0.8)
	})




 
gos<-reactive({

	if(input$kgset==TRUE){
		kar<-plotdata()
		dfg<-data.frame(id=input$inSelectg)
		dfg<-data.frame(id=dfg[,1])
		kar1<-merge(kar,dfg,by="id")
	}else{
		kar<-plotdata()
		kar1<-subset(kar,kar$start>input$inslider[1] & kar$end < input$inslider[2])}
		mdat<-merge(sgo,kar1,by="id",sort=F)
		return(mdat)
	})

output$go_mapfirst<-renderPlot({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		gdat<-gos()
		gdat1<- gdat[,c(1,4)]
		d<-dist(table(gdat1)) 
		fit <- hclust(d) 
		plot(fit, main="Go Cluster", xlab=NULL)
	if(input$fitchk==TRUE){
		rect.hclust(fit, k=input$fitnumb, border="red")}else{return()}
	})

output$go_maptab<-renderTable({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		gdat<-gos()
		gdat1<- gdat[,c(1,4)]
		t(table(gdat1))
	})



output$go_map<-renderPlot({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		gdat<-gos()
		gdat1<- gdat[,c(1,3)]
		gdat2<- princomp(t(table(gdat1)),cor=T)
		plot(gdat2,type="lines") 
	})

output$go_map1<-renderPlot({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		mdat<-gos()
		gdat<-subset(mdat,mdat[,2]==input$gotypep)
		gdat1<- gdat[,c(1,3)]
		gdat2<- princomp(t(table(gdat1)),cor=T)
		biplot(gdat2,cex =0.8,expand=1)
	})
	
output$go_sum<-renderPrint({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		mdat<-gos()
		gdat<-subset(mdat,mdat[,2]==input$gotypep)
		gdat1<- gdat[,c(1,3)]
		summary(princomp(t(table(gdat1)), cor = TRUE,scores = TRUE))
	})
	
output$go_fit<-renderPrint({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		mdat<-gos()
		gdat<-subset(mdat,mdat[,2]==input$gotypep)
		gdat1<- gdat[,c(1,3)]
		fit<-(princomp(t(table(gdat1)), cor = TRUE,scores = TRUE))
		fit$scores
	})
	
output$go_load<-renderPrint({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		mdat<-gos()
		gdat<-subset(mdat,mdat[,2]==input$gotypep)
		gdat1<- gdat[,c(1,3)]
		loadings(princomp(t(table(gdat1)), cor = TRUE,scores = TRUE))
	})

output$go_pie<-renderPlot({
	if (is.null(input$inslider)){ return(NULL)}
	if (is.null(marr())){ return(NULL)}
		mdat<-gos()
		mdat<-subset(mdat,mdat[,2]==input$gotypep)
		mdat1<- mdat[,c(1,3)]
		mdat2<- table(mdat1[,2])
		pie(mdat2,labels=names(mdat2), radius = 0.5,cex = 0.8)
	})


updateButton(session, "btngrp1", style = "info", size = "small")


output$myImage <- renderPlot({
data1<-plotdata()
if (is.null(input$inslider)){ return(NULL)}
else if(input$inslider[2]-input$inslider[1]>80000){
createAlert(session, anchorID = "alert_gene", content = "Try to select region with in 80000 range!",append = FALSE )}else{
input1 <-input$inslider[1]
input2 <-input$inslider[2]
	  
	if (is.null(input$inCheckboxGroup)){cogclr<-NULL}else{cogclr<-input$inCheckboxGroup}
	if (as.character(input$gridr)=="TRUE"){grir<-"TRUE"}else{grir<-"FALSE"}
mplot<-plotter(data1,trunit,input1,input2,0.2,1.9,"grey","grey",5,cogclr,grir)
observe({
     data1<-subset(data1,data1$start>input1 & data1$start<input2)
  	updateSelectInput(session=session, "gselect", choices=paste(data1$id,'----',data1$name, sep=""),selected=data1[1,])})
	  
	  
	      return(mplot)  
	  
	  
}})

output$tabval<- renderTable({


})




output$myImage2 <- renderPlot({
data1<-plotdata()
if (is.null(input$inslider)){ return(NULL)}
input1 <-input$inslider[1]
input2 <-input$inslider[2]
	if (is.null(input$tracklist)){trk<-NULL}else{trk<-get(paste(input$tracklist))}
    if (as.character(input$gridr)=="TRUE"){grir<-"TRUE"}else{grir<-"FALSE"}
mplot<-track(trk,input1,input2,5,grir)
return(mplot)

})

addPopover(session, id = "searchhelp", title = "Gene Search Box", content = "Use Synechocystic gene symbol. example: sll0528", placement = "right", trigger = "hover")



  gene1<-reactive({
    data1<-plotdata()
	if(is.null(input$kar)){return()}
  rownames(data1)<-1:nrow(data1)
if(input$kar[2]<1.1 &input$kar[2]>-1.1){
  gene1<-subset(data1,data1$start<input$kar[1] & data1$end >input$kar[1])
}else{return()}
return(gene1)
  })
  
 output$tabval<-renderText({
if(is.null(input$instart)){return()}
kar1<-subset(kar,kar$start>input$inslider[1] & kar$end < input$inslider[2])
unique(kar1[,1]) })
  
  
  
 output$content <- renderText({
gene1<-gene1()
  ou1<-paste("<div id='genebox' style='left:550px; top:180px; border:#c1c1c1; padding:10px; border:1px;'><b>",gene1[,c(1)],"</b>:",gene1[,c(6)],"</div>")
  updateTextInput(session=session, "geneid", value=gene1[,c(1)])
  return(ou1)
  })
  
  output$myImage3 <- renderPlot({
data1<-plotdata()
if (is.null(input$inslider)){ return(NULL)}
input1 <-input$inslider[1]
input2 <-input$inslider[2]
    if (as.character(input$gridr)=="TRUE"){grir<-"TRUE"}else{grir<-"FALSE"}
		
mplot<-gc_content(data1,input1,input2,grir,as.numeric(input$gcnum))
return(mplot)

})
  
  
  datac <- reactive({

    dist <- switch(input$variable,
"Acaryochloris marina MBIC11017_chr" =158333233,
"Acaryochloris marina MBIC11017_pREB1" =158339488,
"Acaryochloris marina MBIC11017_pREB2" =158339871,
"Acaryochloris marina MBIC11017_pREB3" =158340280,
"Acaryochloris marina MBIC11017_pREB4" =158340643,
"Acaryochloris marina MBIC11017_pREB5" =158340917,
"Acaryochloris marina MBIC11017_pREB6" =158341140,
"Acaryochloris marina MBIC11017_pREB7" =158341329,
"Acaryochloris marina MBIC11017_pREB8" =158341503,
"Acaryochloris marina MBIC11017_pREB9" =158341621,
"Anabaena cylindrica PCC 7122_chr" =440679730,
"Anabaena cylindrica PCC 7122_pANACY.01" =440685051,
"Anabaena cylindrica PCC 7122_pANACY.03" =440685202,
"Anabaena cylindrica PCC 7122_pANACY.04" =440685330,
"Anabaena cylindrica PCC 7122_pANACY.06" =440685390,
"Anabaena cylindrica PCC 7122_pANACY.05" =436843263,
"Anabaena cylindrica PCC 7122_pANACY.02" =443329695,
"Anabaena sp. 90_cir" =414075311,
"Anabaena sp. 90_pANA01" =414078960,
"Anabaena sp. 90_pANA02" =414079049,
"Anabaena sp. 90_chr" =414079097,
"Anabaena sp. 90_pANA03" =414079804,
"Anabaena variabilis ATCC 29413_plasmid A" =75812284,
"Anabaena variabilis ATCC 29413_plasmid B" =75812629,
"Anabaena variabilis ATCC 29413_plasmid C" =75812661,
"Anabaena variabilis ATCC 29413_chr" =75906225,
"Anabaena variabilis ATCC 29413_ins element" =292905252,
"Arthrospira platensis NIES-39_chr" =479126419,
"Calothrix sp. PCC 6303_pCAL6303.01" =428303471,
"Calothrix sp. PCC 6303_pCAL6303.02" =428303567,
"Calothrix sp. PCC 6303_chr" =428296779,
"Calothrix sp. PCC 6303_pCAL6303.03" =428303526,
"Calothrix sp. PCC 7507_chr" =427715354,
"Chamaesiphon minutus PCC 6605_chr" =434384225,
"Chamaesiphon minutus PCC 6605_pCHA6605.02" =434389900,
"Chamaesiphon minutus PCC 6605_pCHA6605.01" =436736818,
"Chroococcidiopsis thermalis PCC 7203_chr" =428205073,
"Chroococcidiopsis thermalis PCC 7203_pCHRO.02" =428210482,
"Chroococcidiopsis thermalis PCC 7203_pCHRO.01" =428204730,
"Crinalium epipsammum PCC 9333_pCRI9333.01" =428308248,
"Crinalium epipsammum PCC 9333_pCRI9333.02" =428308292,
"Crinalium epipsammum PCC 9333_pCRI9333.04" =428303616,
"Crinalium epipsammum PCC 9333_pCRI9333.05" =428308343,
"Crinalium epipsammum PCC 9333_pCRI9333.07" =428303659,
"Crinalium epipsammum PCC 9333_chr" =428303693,
"Crinalium epipsammum PCC 9333_pCRI9333.03" =428308168,
"Crinalium epipsammum PCC 9333_pCRI9333.06" =428308218,
"Crinalium epipsammum PCC 9333_pCRI9333.08" =428308241,
"Cyanobacterium aponinum PCC 10605_chr" =428768415,
"Cyanobacterium aponinum PCC 10605_pCYAN10605.01" =428771774,
"Cyanobacterium stanieri PCC 7202_chr" =428771848,
"Cyanobium gracile PCC 6307_chr" =427701340,
"Cyanothece sp. ATCC 51142_plasmid A" =172034820,
"Cyanothece sp. ATCC 51142_plasmid B" =172034859,
"Cyanothece sp. ATCC 51142_plasmid C" =172034885,
"Cyanothece sp. ATCC 51142_plasmid D" =172034905,
"Cyanothece sp. ATCC 51142_cir" =172034917,
"Cyanothece sp. ATCC 51142_linear" =172054853,
"Cyanothece sp. PCC 7424_chr" =218437013,
"Cyanothece sp. PCC 7424_pP742403" =218442241,
"Cyanothece sp. PCC 7424_Pp742404" =218442305,
"Cyanothece sp. PCC 7424_pP742405" =218442325,
"Cyanothece sp. PCC 7424_pP742406" =218442342,
"Cyanothece sp. PCC 7424_Pp742402" =218442435,
"Cyanothece sp. PCC 7424_pP742401" =218442611,
"Cyanothece sp. PCC 7425_pP742501" =219883096,
"Cyanothece sp. PCC 7425_pP742503" =219883409,
"Cyanothece sp. PCC 7425_chr" =220905643,
"Cyanothece sp. PCC 7425_pP742502" =220910610,
"Cyanothece sp. PCC 7822_chr" =307149945,
"Cyanothece sp. PCC 7822_Cy782203" =307149708,
"Cyanothece sp. PCC 7822_Cy782204" =307149679,
"Cyanothece sp. PCC 7822_Cy782205" =307149645,
"Cyanothece sp. PCC 7822_Cy782201" =307591951,
"Cyanothece sp. PCC 7822_Cy782202" =307591289,
"Cyanothece sp. PCC 7822_Cy782206" =307591243,
"Cyanothece sp. PCC 8801_pP880101" =218203844,
"Cyanothece sp. PCC 8801_pP880102" =218203933,
"Cyanothece sp. PCC 8801_chr" =218244892,
"Cyanothece sp. PCC 8801_pP880103" =218249153,
"Cyanothece sp. PCC 8802_pP880201" =256818682,
"Cyanothece sp. PCC 8802_chr" =257057919,
"Cyanothece sp. PCC 8802_pP880202" =257062240,
"Cyanothece sp. PCC 8802_pP880203" =256823861,
"Cyanothece sp. PCC 8802_pP880204" =256823891,
"Cylindrospermum stagnale PCC 7417_pCYLST.02" =434389956,
"Cylindrospermum stagnale PCC 7417_chr" =434402184,
"Cylindrospermum stagnale PCC 7417_pCYLST.03" =434408147,
"Cylindrospermum stagnale PCC 7417_pCYLST.01" =436669979,
"Dactylococcopsis salina PCC 8305_chr" =428778395,
"Geitlerinema sp. PCC 7407_chr" =428223461,
"Gloeobacter kilaueensis JS1_chr" =554634310,
"Gloeobacter violaceus PCC 7421_chr" =37519569,
"Gloeocapsa sp. PCC 7428_chr" =434390831,
"Gloeocapsa sp. PCC 7428_pGLO7428.01" =434395721,
"Gloeocapsa sp. PCC 7428_pGLO7428.04" =434395856,
"Gloeocapsa sp. PCC 7428_pGLO7428.03" =434408186,
"Gloeocapsa sp. PCC 7428_pGLO7428.02" =436735872,
"Halothece sp. PCC 7418_chr" =428774686,
"Leptolyngbya sp. PCC 7376_chr" =427722021,
"Microcoleus sp. PCC 7113_chr" =428308377,
"Microcoleus sp. PCC 7113_pMIC7113.01" =428314405,
"Microcoleus sp. PCC 7113_pMIC7113.03" =428314593,
"Microcoleus sp. PCC 7113_pMIC7113.05" =428314670,
"Microcoleus sp. PCC 7113_pMIC7113.06" =428314511,
"Microcoleus sp. PCC 7113_pMIC7113.08" =428314710,
"Microcoleus sp. PCC 7113_pMIC7113.02" =428314740,
"Microcoleus sp. PCC 7113_pMIC7113.04" =428314544,
"Microcoleus sp. PCC 7113_pMIC7113.07" =428314808,
"Microcystis aeruginosa NIES-843_chr" =166362741,
"'Nostoc azollae' 0708_chr" =298489614,
"'Nostoc azollae' 0708_pAzo01" =298501383,
"'Nostoc azollae' 0708_pAzo02" =298501435,
"Nostoc punctiforme PCC 73102_chr" =186680550,
"Nostoc punctiforme PCC 73102_pNPUN05" =186686638,
"Nostoc punctiforme PCC 73102_pNPUN03" =186686664,
"Nostoc punctiforme PCC 73102_pNPUN01" =186686738,
"Nostoc punctiforme PCC 73102_pNPUN02" =186687048,
"Nostoc punctiforme PCC 73102_pNPUN04" =186686999,
"Nostoc sp. PCC 7107_chr" =427705465,
"Nostoc sp. PCC 7120_pCC7120beta" =17158637,
"Nostoc sp. PCC 7120_pCC7120zeta" =17158061,
"Nostoc sp. PCC 7120_pCC7120gamma" =17227374,
"Nostoc sp. PCC 7120_pCC7120epsilon" =17227465,
"Nostoc sp. PCC 7120_chr" =17227497,
"Nostoc sp. PCC 7120_pCC7120delta" =17232874,
"Nostoc sp. PCC 7120_pCC7120alpha" =17233017,
"Nostoc sp. PCC 7524_pNOS7524.01" =427726395,
"Nostoc sp. PCC 7524_chr" =427727289,
"Nostoc sp. PCC 7524_pNOS7524.02" =427726360,
"Oscillatoria acuminata PCC 6304_chr" =428210537,
"Oscillatoria acuminata PCC 6304_pOSCIL6304.02" =428216242,
"Oscillatoria acuminata PCC 6304_pOSCIL6304.01" =428210485,
"Oscillatoria nigro-viridis PCC 7112_chr" =428315182,
"Oscillatoria nigro-viridis PCC 7112_pOSC7112.02" =428314827,
"Oscillatoria nigro-viridis PCC 7112_pOSC7112.03" =428315019,
"Oscillatoria nigro-viridis PCC 7112_pOSC7112.05" =428320964,
"Oscillatoria nigro-viridis PCC 7112_pOSC7112.01" =428320975,
"Oscillatoria nigro-viridis PCC 7112_pOSC7112.04" =428315127,
"Pleurocapsa sp. PCC 7327_chr" =428200461,
"Prochlorococcus marinus str. AS9601_chr" =123967536,
"Prochlorococcus marinus str. MIT 9211_chr" =159902540,
"Prochlorococcus marinus str. MIT 9215_chr" =157412338,
"Prochlorococcus marinus str. MIT 9301_chr" =126695337,
"Prochlorococcus marinus str. MIT 9303_chr" =124021714,
"Prochlorococcus marinus str. MIT 9312_chr" =78778385,
"Prochlorococcus marinus str. MIT 9313_chr" =33862273,
"Prochlorococcus marinus str. MIT 9515_chr" =123965234,
"Prochlorococcus marinus str. NATL1A_chr" =124024712,
"Prochlorococcus marinus subsp. marinus str. CCMP1375_chr" =33239452,
"Prochlorococcus marinus subsp. pastoris str. CCMP1986_chr" =33860560,
"Pseudanabaena sp. PCC 7367_pPSE7367.01" =428219846,
"Pseudanabaena sp. PCC 7367_chr" =428216284,
"Rivularia sp. PCC 7116_chr" =427733619,
"Rivularia sp. PCC 7116_pRIV7116.01" =427740400,
"Rivularia sp. PCC 7116_pRIV7116.02" =427740454,
"Stanieria cyanosphaera PCC 7437_chr" =434396591,
"Stanieria cyanosphaera PCC 7437_pSTA7437.02" =434401063,
"Stanieria cyanosphaera PCC 7437_pSTA7437.03" =434401190,
"Stanieria cyanosphaera PCC 7437_pSTA7437.01" =434408264,
"Stanieria cyanosphaera PCC 7437_pSTA7437.05" =434408556,
"Stanieria cyanosphaera PCC 7437_pSTA7437.04" =436736068,
"Synechococcus elongatus PCC 6301_chr" =56750010,
"Synechococcus elongatus PCC 7942_plasmid 1" =81230333,
"Synechococcus elongatus PCC 7942_chr" =81298811,
"Synechococcus sp. CC9311_chr" =113952711,
"Synechococcus sp. CC9605_chr" =78211558,
"Synechococcus sp. CC9902_chr" =78183584,
"Synechococcus sp. JA-2-3B'a(2-13)_chr" =86607503,
"Synechococcus sp. JA-3-3Ab_chr" =86604733,
"Synechococcus sp. PCC 6312_chr" =427711179,
"Synechococcus sp. PCC 6312_pSYN6312.01" =427714745,
"Synechococcus sp. PCC 7002_pAQ7" =170076470,
"Synechococcus sp. PCC 7002_chr" =170076636,
"Synechococcus sp. PCC 7002_pAQ1" =170079460,
"Synechococcus sp. PCC 7002_pAQ3" =170079464,
"Synechococcus sp. PCC 7002_pAQ4" =170079482,
"Synechococcus sp. PCC 7002_pAQ5" =170079513,
"Synechococcus sp. PCC 7002_pAQ6" =170079553,
"Synechococcus sp. PCC 7502_pSYN7502.01" =428223389,
"Synechococcus sp. PCC 7502_pSYN7502.02" =428223448,
"Synechococcus sp. PCC 7502_chr" =428220140,
"Synechococcus sp. RCC307_chr" =148241099,
"Synechococcus sp. WH 7803_chr" =148238336,
"Synechococcus sp. WH 8109_chr" =33864539,
"Synechocystis sp. PCC 6803_PsysM" =38505535,
"Synechocystis sp. PCC 6803_pSYSA" =38505668,
"Synechocystis sp. PCC 6803_PsysG" =38505775,
"Synechocystis sp. PCC 6803_PsysX" =38505825,
"Synechocystis sp. PCC 6803_chr" =384435229,
"Synechocystis sp. PCC 6803_chr" =451813329,
"Synechocystis sp. PCC 6803_pSYSA_M" =451816673,
"Synechocystis sp. PCC 6803_pSYSG_M" =451816489,
"Synechocystis sp. PCC 6803_pCA2.4_M" =451816539,
"Synechocystis sp. PCC 6803_pCC5.2_M" =451816781,
"Synechocystis sp. PCC 6803_pSYSM_M" =451816542,
"Synechocystis sp. PCC 6803_pSYSX_M" =451816787,
"Synechocystis sp. PCC 6803_pCB2.4_M" =451816777,
"Synechocystis sp. PCC 6803 substr. GT-I_chr" =383489963,
"Synechocystis sp. PCC 6803 substr. PCC-N_chr" =383324079,
"Synechocystis sp. PCC 6803 substr. PCC-P_chr" =383320909,
"Thermosynechococcus elongatus BP-1_chr" =22297544,
"Thermosynechococcus sp. NK55_chr" =571026788,
"Trichodesmium erythraeum IMS101_chr" =113473942)
  dat1<-subset(sid,sid$gi==dist)
  dat2<- dat1[order(dat1$start2),]
  return(dat2)
  })
  
output$comp_plot<-renderPlot({
if (input$comp_but==0){ return(NULL)}
	input$comp_but
	isolate({
plot(x=1, y=1, type="h", col="white", xlim=c(as.numeric(input$cone),as.numeric(input$ctwo)), ylim=c(-7,8), xlab="", ylab="", yaxt="n",xaxt="n",bty="n",main=paste("Genome level comparison between Synechocystis Genome VS", input$variable), sub=paste("Synechocystis Genome VS", input$variable))
segments(as.numeric(input$cone), 2, as.numeric(input$ctwo), 2, col="#D00000", lwd=1)
segments(as.numeric(input$cone), -2, as.numeric(input$ctwo), -2, col="#3333CC", lwd=1)
datac<-datac()
rect(datac$start2, -1.5, datac$end2, 1.5, xpd=FALSE , border="#ffffff", col="#cccccc20")

rect(datac$start2, 1.8, datac$end2, 2.5, xpd=FALSE , border="#ffffff", col="#D0000070")
rect(datac$start2, -2.5, datac$end2, -1.8, xpd=FALSE , border="#ffffff", col="#3333CC70")

if(as.numeric(input$ctwo)-as.numeric(input$cone)>8000){
text((datac$start2+datac$end2)/2,-5.2,labels=paste(datac$start1, "-", datac$end1,sep=" "),cex = 1, srt = 90,col = "#666666")
text((datac$start2+datac$end2)/2,5.2,labels=paste(datac$start2, "-", datac$end2,sep=" "),cex = 1, srt = 90,col = "#666666")}
else{
text((datac$start2+datac$end2)/2,-5.2,labels=paste(datac$start1, "-", datac$end1,sep=" "),cex = 1, srt = 0,col = "#666666")
text((datac$start2+datac$end2)/2,5.2,labels=paste(datac$start2, "-", datac$end2,sep=" "),cex = 1, srt = 0,col = "#666666")
}

legend("topright",legend=c("Synechocystis 6803_chr", input$variable),lty=1,lwd=2,col=c("#D00000","#3333CC"),ncol=2,bty="n",cex=1.2,text.col=c("black"),inset=0.01)
})
  })
    output$comp_table<- renderDataTable({ 
  datac<-datac()[,c(2,3,5,6,7,8,10)]
   
  
    },options = list(aLengthMenu = c(5, 10, 15,20,25,50), iDisplayLength = 10)) 
  
  
observe({
     pathdata<-unique(path$specificpath)
  	updateSelectInput(session=session, "pathselect", choices=pathdata)

})


output$path_table<- renderDataTable({ 
  pathdata1<-subset(path,path$specificpath==input$pathselect)
  pathm<-merge(metapath,pathdata1,by="id")
  pathm[,c(1,2,4,5)]
    
    },options = list(aLengthMenu = c(5, 10, 15,20,25,50), iDisplayLength = 10))  
	
	
  

 blastdata<-reactive({
 	if (input$blast_data==0){ return(NULL)}
	input$blast_data
	
	isolate({
	bdat<-read.table("bdata.txt", as.is=T, header=T)
  
 })
   
 }) 
  
  
output$blastplot<-renderPlot({
if (input$blast_data==0){ return(NULL)}
bdata<- blastdata()
plot(x=1, y=1, type="h", col="white", xlim=c(as.numeric(bdata$s_start[1]),as.numeric(bdata$s_end[1])), ylim=c(nrow(bdata),0), xlab="", ylab="", yaxt="n",xaxt="n",bty="n")

for (i in 1:nrow(bdata)){
segments(as.numeric(bdata$q_start[1]), i+1, as.numeric(bdata$q_end[1]), i+1, col="#cccccc80", lwd=1)
}

colr<-list()
for (i in 1:nrow(bdata)){
if(bdata$bit_score[i] <40){colr[i]<-"#000000"}
else if(bdata$bit_score[i] >40 & bdata$bit_score[i] <50){colr[i]<-"#3366FF"}
else if(bdata$bit_score[i] > 50 & bdata$bit_score[i] < 80){colr[i]<-"#336633"}
else if(bdata$bit_score[i] > 80 & bdata$bit_score[i] < 200){colr[i]<-"#CC3399"}
else if(bdata$bit_score[i] >200){colr[i]<-"#CC0033"}
}
colr<-do.call(rbind,colr)

for (i in 1:nrow(bdata)){
rect(bdata$q_start[i], i-0.3+1, bdata$q_end[i], i+0.3+1, xpd=FALSE , border="#ffffff", col=colr[i,])
}
  })

  
#-----------------------------------------------------plot end--------------------------------------------------------------   


  output$blasttable<- renderDataTable({ 
  blastdata()

  
    },options = list(aLengthMenu = c(5, 10, 15,20,25,50), iDisplayLength = 10)) 




output$seq_bse<- renderText({
        kar<-plotdata()
 		frag <- Views(fs[[1]], as.numeric(input$inslider[1]) ,as.numeric(input$inslider[2]))
		if(input$gen_seq==FALSE & input$gen_seqr==FALSE){
		seq_data1<-as.character(frag)}
		else if(input$gen_seq==TRUE & input$gen_seqr==FALSE){
		seq_data1<-as.character(frag)
		kar1<-subset(kar,kar$start>input$inslider[1] & kar$end<input$inslider[2])
		kar2<-data.frame(id=kar1$id,start=kar1$start-input$inslider[1],end=kar1$end-input$inslider[1],strand=kar1$strand,name=kar1$name)
		paste(">",kar2$id,":",kar2$strand,":",kar2$name,"\n",substring(seq_data1, kar2$start, kar2$end),"\n")}
	else if(input$gen_seq==TRUE & input$gen_seqr==TRUE){
		seq_data1<-as.character(frag)
		kar1<-subset(kar,kar$start>input$inslider[1] & kar$end<input$inslider[2])
		kar2<-data.frame(id=kar1$id,start=kar1$start-input$inslider[1],end=kar1$end-input$inslider[1],strand=kar1$strand,name=kar1$name)
		kar3<-list()
		kar4<-list()
		for(i in 1:nrow(kar2)){if(kar2[i,4]=="+"){
		kar3[i]<-paste(">",kar2[i,1],":",kar2[i,4],":",kar2[i,5],"\n",substring(seq_data1, kar2[i,2], kar2[i,3]),"\n")}
		else{kar4[i]<-paste(">",kar2[i,1],":",kar2[i,4],":",kar2[i,5],"\n",Biostrings::reverseComplement(DNAString(substring(seq_data1, kar2[i,2], kar2[i,3]))),"\n")}}
		kar3<-do.call(rbind,kar3)
		kar4<-do.call(rbind,kar4)
		rbind(kar3,kar4)}

 })
 
 output$seq_bse1<- renderText({
        kar<-plotdata()
		if(input$infuse1==0){return()}
		input$infuse1
		isolate({
 		frag <- Views(fs[[1]], as.numeric(input$range_start) ,as.numeric(input$range_end))
		seq_data1<-as.character(frag)
	 })
  })
 

 promextractdata<-reactive({
if (input$prominfuse1==0){ return(NULL)}
	input$prominfuse1
if(input$ckbinprom==FALSE){
	data<-kar}else{
	gdb<-data.frame(id=input$inSelectg)
	data=merge(kar,gdb,by="id")}
	isolate({
plus<-subset(data,data$strand=="+")
plus<-subset(plus,plus$id!="slr0611")
minus<-subset(data,data$strand=="-")
pp<-data.frame(id=plus$id,strand=plus$strand,start=plus$start-as.numeric(input$promrange_end),end=plus$start-as.numeric(input$promrange_start))
mp<-data.frame(id=minus$id,strand=minus$strand,start=minus$start+as.numeric(input$promrange_start),end=minus$start+as.numeric(input$promrange_end))
promp<-Views(fs[[1]],pp$start,pp$end)
promn<-Views(fs[[1]],mp$start,mp$end)

pp1<-data.frame(pp,seq=as.character(promp))
mp1<-data.frame(mp,seq=as.character(promn))

prom<-rbind(pp1,mp1)
op<-prom[order(prom[,3]), ]


return(op)
   })
  
 }) 

  output$promtable1<- renderDataTable({ 
  promextractdata()

  
    },options = list(aLengthMenu = c(5, 10, 15,20,25,50), iDisplayLength = 10)) 
 
   output$promextractdown<- downloadHandler(

    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
		  paste(input$dataset, input$filetype2, sep = ".")
	  },

    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      sep <- switch(input$filetype2, "csv" = ",", "txt" = "\t")

      # Write to a file specified by the 'file' argument
      write.table(promextractdata(), file, sep = sep,
        row.names = FALSE)
    }
  )
#--------------------------------------------------------------- sequence extract end------------------------------
#--------------------------------------------------------------- Pattern search start------------------------------
 
 
 
 
 
 
 

plotCov <- function(cov, start, end, mycol="red") {
               plot.new()
               plot.window(c(start, end), c(0, 3))
               axis(1)
               lines(start:end, cov[start:end], type="l", col=mycol)} 
 
 gpromploter<-reactive({
 	if (input$patsearch==0){ return(NULL)}
	input$patsearch
	
	isolate({
	pat<-as.character(input$pattern_ser1)
mysearch <-matchPattern(pat, fs[[1]] , fixed=input$ambig,max.mismatch=as.numeric(input$mm1), min.mismatch=as.numeric(input$mm1))

mycov <- as.integer(coverage(mysearch, 1, length(fs[[1]])))
  
 })
   
 })
  gpromploter1<-reactive({
 	if (input$patsearch==0){ return(NULL)}
	input$patsearch
	
	isolate({
	pat<-Biostrings::complement(DNAString(as.character(input$pattern_ser1)))
mysearch <-matchPattern(pat, fs[[1]] , fixed=input$ambig,max.mismatch=as.numeric(input$mm1), min.mismatch=as.numeric(input$mm1))

mycov <- as.integer(coverage(mysearch, 1, length(fs[[1]])))
  
 })
   
 })
 output$patplot<- renderPlot({
 	if (input$patsearch==0){ return(NULL)}
	input$patsearch
	
	isolate({
	mycov<- gpromploter()
	mycov1<- gpromploter1()
if(input$comp_std==FALSE){
plotCov(mycov, start=input$patt_genlen_f, end=input$patt_genlen_t, mycol="blue")}
else{
op <- par(mfrow=c(2,1))
plotCov(mycov, start=input$patt_genlen_f, end=input$patt_genlen_t, mycol="blue")
plotCov(mycov1, start=input$patt_genlen_f, end=input$patt_genlen_t, mycol="red")
par(op)
}

  })
  })
  
  output$myImage4<- renderPlot({
 	if (input$pathook==0){ return(NULL)}
	input$pathook
	
	isolate({
	mycov<- gpromploter()
	mycov1<- gpromploter1()
if(input$comp_std==FALSE){
plotCov(mycov, start=input$inslider[1], end=input$inslider[2], mycol="blue")}
else{
par(mfrow=c(2,1))
plotCov(mycov, start=input$patt_genlen_f, end=input$patt_genlen_t, mycol="blue")
plotCov(mycov1, start=input$patt_genlen_f, end=input$patt_genlen_t, mycol="red")

}


  })
  }) 
  
 observe({
if(input$patty=="patty_2")
updateTextInput(session=session,"pattern_ser1", value=input$pattern_ser2)
})  
  
 
patttab1<-reactive({
if (input$patsearch==0){ return(NULL)}
	input$patsearch
	
	isolate({
	if(input$comp_std==FALSE){
	pat<-as.character(input$pattern_ser1)
mysearch <-matchPattern(pat, fs[[1]] , fixed=input$ambig,max.mismatch=as.numeric(input$mm1), min.mismatch=as.numeric(input$mm1))
  rang<-data.frame(start=unlist(start(mysearch)),end=unlist(end(mysearch)), pattern=as.character(mysearch))
  subset(rang,rang$start>input$patt_genlen_f & rang$start<input$patt_genlen_t)}
  else{
  pat<-as.character(input$pattern_ser1)
  pat1<-Biostrings::complement(DNAString(as.character(input$pattern_ser1)))
  mysearch <-matchPattern(pat, fs[[1]] , fixed=input$ambig,max.mismatch=as.numeric(input$mm1), min.mismatch=as.numeric(input$mm1))
  rang<-data.frame(start=unlist(start(mysearch)),end=unlist(end(mysearch)), pattern=as.character(mysearch))
  po<-subset(rang,rang$start>input$patt_genlen_f & rang$start<input$patt_genlen_t)
  mysearch1 <-matchPattern(pat1, fs[[1]] , fixed=input$ambig,max.mismatch=as.numeric(input$mm1), min.mismatch=as.numeric(input$mm1))
  rang1<-data.frame(start=unlist(start(mysearch1)),end=unlist(end(mysearch1)), pattern=as.character(mysearch1))
  pt<-subset(rang1,rang1$start>input$patt_genlen_f & rang1$start<input$patt_genlen_t)
  rbind(po,pt)
  }
  
  })


}) 
  
     output$pattable<- renderTable({ 

     pattab<-patttab1()

    })
	
	output$seqlogoplot<- renderPlot({ 
     
     pattab<-patttab1()
	 
     patterns<-DNAStringSet(pattab$pattern)
	 if (nrow(pattab)==0){ return(NULL)}
	 pwm <- PWM(patterns)
	 
	 seqLogo(t(t(pwm) * 1/colSums(pwm)), ic.scale=FALSE)
    })
	
	output$seqlogomat<- renderTable({ 
     pattab<-patttab1()
	 
     patterns<-DNAStringSet(pattab$pattern)
	 if (nrow(pattab)==0){ return(NULL)}
	 pwm <- PWM(patterns)
	 data.frame(t(abs(pwm)))
	 
    })
	
	
  
 promdata<-reactive({
if (input$patsearch==0){ return(NULL)}
	input$patsearch
if(input$pgset=="FALSE"){
	data<-kar}else{
	gdb<-data.frame(id=input$inSelectg)
	data=merge(kar,gdb,by="id")}
	isolate({
pat<-as.character(input$pattern_ser1)
pat1<-Biostrings::complement(DNAString(as.character(input$pattern_ser1)))
plus<-subset(data,data$strand=="+")
minus<-subset(data,data$strand=="-")
pp<-data.frame(id=plus$id,start=plus$start-as.numeric(input$patt_promlen_t),end=plus$start)
mp<-data.frame(id=minus$id,start=minus$start,end=minus$start+as.numeric(input$patt_promlen_t))
prom<-rbind(pp,mp)
op<-prom[order(prom[,2]), ]
proms<-Views(fs[[1]],op$start,op$end)
seq<-data.frame(seq=as.character(proms[2:length(proms),]))
subject<-DNAStringSet(seq$seq)
names(subject)<-data[2:nrow(data),1]
if(input$comp_std==FALSE){

mysearch <-vmatchPattern(pat, subject , fixed=input$ambig)
count_index <- countIndex(mysearch)
with_match <- which(count_index >= 1)

res1<-data.frame(unlist(mysearch))
res2<-data.frame(names=names(subject[with_match]),seq=as.character(subject[with_match]),pat=rep(as.character(pat),length(with_match)))
res<-merge(res1,res2,by="names")
}
else{
mysearch <-vmatchPattern(pat, subject , fixed=input$ambig)
count_index <- countIndex(mysearch)
with_match <- which(count_index >= 1)
res1<-data.frame(unlist(mysearch))
res2<-data.frame(names=names(subject[with_match]),seq=as.character(subject[with_match]),pat=rep(as.character(pat),length(with_match)))
pats1<-merge(res1,res2,by="names")
pat1<-Biostrings::complement(DNAString(as.character(input$pattern_ser1)))
mysearch1 <-vmatchPattern(pat1, subject , fixed=input$ambig)
count_index <- countIndex(mysearch1)
with_match <- which(count_index >= 1)
res4<-data.frame(unlist(mysearch))
res5<-data.frame(names=names(subject[with_match]),seq=as.character(subject[with_match]),pat=rep(as.character(pat1),length(with_match)))
pats2<-merge(res4,res5,by="names")

res<-rbind(pats1,pats2)

}
return(res)
   })
  
 }) 

  output$pattable1<- renderDataTable({ 
 promdata()

  
    },options = list(aLengthMenu = c(5, 10, 15,20,25,50), iDisplayLength = 10))
	
  output$patplotden<- renderPlot({ 
 denp<-promdata()
d<-density(-(denp$start+denp$end)/2)
  plot(d,main="Density of Pattern")
  polygon(d, col="grey")
  rug(-(denp$start+denp$end)/2, col="red")
    })	
	

	
  output$prompttData <- downloadHandler(

    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
		  paste(input$dataset, input$filetype1, sep = ".")
	  },

    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      sep <- switch(input$filetype1, "csv" = ",", "txt" = "\t")

      # Write to a file specified by the 'file' argument
      write.table(promdata(), file, sep = sep,
        row.names = FALSE)
    }
  )
	
    observe({ 
	if (input$infuse==0){ return(NULL)}
	input$infuse
	isolate({
	updateNumericInput(session=session, "patt_genlen_f", value=input$inslider[1])
	updateNumericInput(session=session, "patt_genlen_t", value=input$inslider[2])
	}) 
  })	
  
  output$sliders3 <- renderUI({
  sliderInput("mm1",label="",min=0,max=5,value=0,step=1)
  
  })
 
 

#-----------------------------------------------------pattern serach end-------------------------------------------------------------- 

#-----------------------------------------------------dataset start-------------------------------------------------------------- 

  
    uploaddata<-reactive({
	inFile <- input$file1

    if (is.null(inFile))
      return(NULL)
    
   dat1<-read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
   as.data.frame(sapply(dat1, function(x) gsub("\"", "", x)))
	
	})
  
  
  
  
  output$contents <- renderTable({
    
uploaddata()
  })




datasetInput <- reactive({

    karr<-switch(input$dataset,
           "Gene" = kar,
           "Promoter" = tpr,
           "TATA box" = tata,
		   "Antisense-RNA" = asrna,
		   "Noncoding-RNA"=ncrna,
		   "HIPs"=p1,
		   "5'UTR"=utr,"Intergenic sequence"=intern,"Transcription unit"=trunit)
res<-subset(karr, karr[,3]>input$instart & karr[,4]<input$inend )
return(res)
  })
  
  output$glist <- renderDataTable({
    datasetInput()
	
  },options = list(aLengthMenu = c(5, 10, 15,20,25,50), iDisplayLength = 10))
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  #   it should write out data to that filename.
  output$downloadData <- downloadHandler(

    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
		  paste(input$dataset, input$filetype, sep = ".")
	  },

    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")

      # Write to a file specified by the 'file' argument
      write.table(datasetInput(), file, sep = sep,
        row.names = FALSE)
    }
  )
  
u_data<- reactive({ 
if (input$add_basket==0){ return(NULL)}
if (is.null(input$gselect)){ return(NULL)}
input$add_basket
isolate({
udat<-data.frame(id=input$gselect)
gser<-list()
for(j in 1:nrow(udat)){
gser[j]<-strsplit(as.character(udat[j,1]),"----")}
gser1<-do.call(rbind,gser)
colnames(gser1)<-c("id","name")
udata<-merge(kar, gser1, by="id")
udata$id
})
})

u_data1<- reactive({ 
if (input$add_basket1==0){ return(NULL)}
if (is.null(input$gselect1)){ return(NULL)}
input$add_basket1
isolate({
udat<-data.frame(id=input$gselect1)
as.character(udat$id)
})
})

u_data2<- reactive({ 
distab<-distabb()
numb<-which(rownames(distab)==input$distt)

udat<-data.frame(id=distab[numb,])
udat$id

})

output$userdata<- renderTable({ 
data.frame(id=u_data())

})  
  
  
output$userdata1<- renderTable({ 
data.frame(id=u_data1())

})    
  
output$userdata4<- renderTable({ 
data.frame(id=u_data2())

})  
  
mglist<-reactive({
if (input$merge_glist==0){ return(NULL)}
input$merge_glist
isolate({
done<-data.frame(id=u_data())
dtwo<-data.frame(id=u_data1())
dthr<-data.frame(id=u_data2())
rid<-rbind(done,dtwo,dthr)
resu<-merge(rid,kar,by="id")
unique(resu)
})
})
  
  output$userdata2<- renderTable({ 
data.frame(id=input$inSelectg)

})  

  output$userdata3<- renderTable({ 
if (input$merge_glist==0){ return(NULL)}
lit1<-mglist()[,c(1,3,4,5,6)]
lit2<-data.frame(id=input$inSelectg)
merge(lit1,lit2,by="id")

})
 observe({
if (input$merge_glist==0){ return(NULL)}
input$merge_glist
isolate({ 
 updateSelectInput(session, "inSelectg", choices = as.character(mglist()[,1]),selected=as.character(mglist()[,1])) 
 })
 })
  
 observe({
if (input$basket1==0){ return(NULL)}
input$basket1
isolate({   
updateTabsetPanel(session, "inTabset", selected = "User Analysis Basket")
  })})
  
  observe({
if (input$ranger1==0){ return(NULL)}
input$ranger1

isolate({   

dfr<-data.frame(id=as.character(u_data2()))
u_data2<-merge(dfr,kar,by="id")
updateNumericInput(session, "range_start", value = as.numeric(min(u_data2$start)))
updateNumericInput(session, "range_end", value = as.numeric(max(u_data2$end)))
  })}) 

  #-----------------------------------------------------dataset end-------------------------------------------------------------- 
observe({
if(input$inSelect3=="userdata"){
if(is.null(mglist())){return(NULL)}else{
updateSelectInput(session, "inSelect2", choices =input$inSelectg)
}}else if(input$inSelect3=="genedata"){
intdata<-subset(kar, kar[,3]>input$instart & kar[,4]<input$inend )
idata<-data.frame(bait=intdata$id)
updateSelectInput(session, "inSelect2", choices = as.character(idata[,1]), selected=as.character(idata[1,]))}
})

karp<-reactive({

idata<-data.frame(bait=input$inSelect2)
interact<-merge(interact,idata,by="bait")
interact<-subset(interact,interact$source==input$int_source)
bind <- cbind(interact$bait, interact$prey)
karp <- bind[, c(1,2)]
return(karp)

})

karpc<-reactive({
karp<-karp()
if(input$ppi_col=="ppic_1"){colr<-cogc}else{colr<-kasc}

colnames(karp)<-c("id","prey")
one<-merge(karp,colr,by="id")
kk<-data.frame(id=one$id,clr=one$clr)
colnames(karp)<-c("bait","id")
two<-merge(karp,colr,by="id")
kk1<-data.frame(id=two$id,clr=two$clr)
kk2<-rbind(kk,kk1)
#res<-cbind(one,two)
return(kk2)
})




 output$int <- renderPlot({
if (input$gobut==0){ return(NULL)}
input$gobut
isolate({
karp<-karp()
nAttrs<-list()
gf<-list(as.character(karpc()$clr))
names(gf[[1]])<-karpc()$id
idat<-data.frame(bait=input$inSelect2)
idata<-list(as.character(rep("yellow",nrow(idat))))
names(idata[[1]])<-idat$bait
nAttrs$fillcolor <- idata[[1]]
nAttrs$fontcolor <-gf[[1]]
nAttrs$color <- gf[[1]]
g <- ftM2graphNEL(karp)
g1 <- layoutGraph(g)
nodes <- buildNodeList(g1)
nodelen<-length(nodes)
if(nodelen[1]>1){
renderGraph(g1)

plot(g1, input$int_type, nodeAttrs = nAttrs, attrs=list(node= list(lwd=5),
           edge=list(color="#cccccc"),
           graph=list(rankdir="LR")))}else{stop(paste("Please add more nodes or try other interaction resource"))}
})
})


  output$inttable <- renderDataTable({
  if (input$gobut==0){ return(NULL)}
input$gobut
isolate({  
karp<-karp()
	colnames(karp)<-c("id","id2")
	res1<-merge(kar,karp,by="id",sort=F)
	res1<-data.frame(res1[,c("id","name")])
	colnames(karp)<-c("id1","id")
	res2<-merge(kar,karp,by="id",sort=F)
	res2<-data.frame(res2[,c("id","name")])
	res<-cbind(res1,res2)
	

	 })
  },options = list(aLengthMenu = c(5, 10, 15,20,25,50), iDisplayLength = 10))
  
  
  
keggp<-reactive({

idata<-data.frame(id=input$inSelect2)
interact<-merge(kasc,idata,by="id")
karp <- interact[, c(1,7,8)]
return(karp)

}) 
  
 output$intk <- renderPlot({ 
if (input$gobutk==0){ return(NULL)}
input$gobutk
isolate({
karp<-keggp()
g <- ftM2graphNEL(as.matrix(karp[,c(1,2)]),edgemode = "directed")
nAgo<-makeNodeAttrs(g)
ag.obj <-agopen(g, recipEdges="combined",nodeAttrs=nAgo, name="",layoutType=input$int_type1)
fill.colors<-karp[,3]
for(i in 1:nrow(karp))
{ag.obj@AgNode[[i]]@fillcolor<-fill.colors[i]
ag.obj@AgNode[[i]]@shape<-"ellipse"
ag.obj@AgNode[[i]]@style<-"dotted"
}
nodes <- buildNodeList(g)
nodelen<-length(nodes)
if(nodelen[1]>1){
plot(ag.obj)}else{stop(paste("Please add more nodes or try other interaction resource"))}
})
})  
  
cogp<-reactive({

idata<-data.frame(id=input$inSelect2)
interact<-merge(sgo,idata,by="id")
karp <- interact[, c(1,2,3,4)]
return(karp)
})  
  
output$intc <- renderPlot({ 
if (input$gobutc==0){ return(NULL)}
input$gobutc
isolate({
karp<-cogp()
g <- ftM2graphNEL(as.matrix(karp[,c("id",input$gotype)]),edgemode="directed")

g1 <- layoutGraph(g)
nodes <- buildNodeList(g1)
nodelen<-length(nodes)

if(nodelen[1]>1){

renderGraph(g1)

plot(g1, input$int_type2, attrs=list(node=list(color=c("lightgreen"),fixedsize=FALSE,shape=c("box")),
           edge=list(color="#cccccc"),
           graph=list(rankdir="LR")))}else{stop(paste("Please add more nodes or try other interaction resource"))}
})
})   
  
  
  output$ctable <- renderDataTable({
  if (input$gobutc==0){ return(NULL)}
input$gobutc
isolate({  
karp<-cogp()
	 })
  },options = list(aLengthMenu = c(5, 10, 15,20,25,50), iDisplayLength = 10))  
  
  
  
  
keggpath<-reactive({

idata<-data.frame(id=input$inSelect2)
interact<-merge(sgo,idata,by="id")
karp <- interact[, c(1,2,3,4)]
return(karp)
})  
  
output$intc <- renderPlot({ 
if (input$gobutc==0){ return(NULL)}
input$gobutc
isolate({
karp<-cogp()
g <- ftM2graphNEL(as.matrix(karp[,c("id",input$gotype)]),edgemode="directed")

g1 <- layoutGraph(g)
nodes <- buildNodeList(g1)
nodelen<-length(nodes)

if(nodelen[1]>1){

renderGraph(g1)

plot(g1, input$int_type2, attrs=list(node=list(color=c("lightgreen"),fixedsize=FALSE,shape=c("box")),
           edge=list(color="#cccccc"),
           graph=list(rankdir="LR")))}else{stop(paste("Please add more nodes or try other interaction resource"))}
})
}) 
  
#-----------------------------------------annotate-------------------------------------
     
	 adata<-reactive({
	 udata<-uploaddata()
     numb<-which(names(udata)==input$dynamic)
     newdf<-data.frame(id=udata[numb])})

	 annodata<-reactive({
	 udata<-uploaddata()
     numb<-which(names(udata)==input$dynamic)
     newdf<-data.frame(id=udata[numb])
	 colnames(newdf)[1] <- "id"
	  restad<-merge(newdf,kar, by='id')
	  return(restad)

	 })



 output$ui <- renderUI({
    if (is.null(input$input_type))
      return()
    udata<- uploaddata()
    # Depending on input$input_style, we'll generate a different
    # UI component and send it to the client.
    switch(input$input_style,
	     
      "Annotate" = selectInput("dynamic", "Choose header",
        choices = colnames(udata)
		),      "subset" = selectInput("dynamic", "Select Genes",
        choices = colnames(udata),
        selected = colnames(udata),
        multiple = TRUE
      )
    )

  })
  
 
    observe({
	  if (input$merged==0){ return(NULL)}
	  input$merged
	    isolate({

 	updateSelectInput(session=session, "gselect1", choices=annodata()[,1])})    })
  
   observe({
  if (input$merged==0){ return(NULL)}
  if(input$input_type=="Annotate"){
  
output$contents2 <- renderDataTable({

input$merged

    isolate({
annodata()

  })

  })}
    })

datasetInputs <- reactive({
if (input$stressbut==0){ return(NULL)}
input$stressbut
	    isolate({
sres<-subset(stressdata,stressdata$stress==input$stressdataset)
  }) })
  
  output$stresstab<- renderDataTable({ 
     datasetInputs()
  },options = list(aLengthMenu = c(5, 10, 15,20,25,50), iDisplayLength = 5)) 
  
  
  output$downData <- downloadHandler(
    filename = function() {
		  paste(input$stressdataset, input$filetype3, sep = ".")
	  },

    content = function(file) {
      sep <- switch(input$filetype3, "csv" = ",", "txt" = "\t")

      write.table(datasetInputs(), file, sep = sep,
        row.names = FALSE)
    }
  )

    output$downDatamerge <- downloadHandler(
    filename = function() {
		  paste("Userdata", input$filetypemgd, sep = ".")
	  },

    content = function(file) {
      sep <- switch(input$filetypemgd, "csv" = ",", "txt" = "\t")
      lit1<-mglist()[,c(1,3,4,5,6)]
lit2<-data.frame(id=input$inSelectg)
mres<-merge(lit1,lit2,by="id")
      write.table(mres, file, sep = sep,
        row.names = FALSE)
    }
  )
  
      output$downloadDatap <- downloadHandler(
    filename = function() {
		  paste("Pathwaydata", input$filetypep, sep = ".")
	  },

    content = function(file) {
      sep <- switch(input$filetypep, "csv" = ",", "txt" = "\t")
  pathdata1<-subset(path,path$specificpath==input$pathselect)
  pathm<-merge(metapath,pathdata1,by="id")
  pres<-pathm[,c(1,2,4,5)]
      write.table(pres, file, sep = sep,
        row.names = FALSE)
    }
  )
    
})
