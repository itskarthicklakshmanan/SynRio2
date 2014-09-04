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
mdata<-read.delim("farray.txt",as.is=T,header=T)
path<-read.delim("path.txt",as.is=T,header=T)
fs <- readDNAStringSet("NC_000911.fna", "fasta")
sgo<-read.delim("gosyn.txt",as.is=T,header=T)
seg<-read.delim("seg.txt",as.is=T,header=F)
dataall <- read.delim("dataalln.txt", as.is=T, header=T)
interact <- read.delim("interactome.txt", as.is=T, header=T)
cogc <- read.delim("syncogclr.txt", as.is=T, header=T)
kasc <- read.delim("kasuzaclr.txt", as.is=T, header=T)
kar <-subset(dataall,dataall[,2]=="gene")
tata <-subset(dataall,dataall[,2]=="TATA")
asrna <-subset(dataall,dataall[,2]=="asRNA")
tp <-subset(dataall,dataall[,2]=="TP")
tpr <-read.delim("promo.txt", as.is=T, header=T)
ncrna <-subset(dataall,dataall[,2]=="ncRNA")
utr <-subset(dataall,dataall[,2]=="5'UTR")
p1 <-subset(dataall,dataall[,2]=="hip1")
intern <-subset(dataall,dataall[,2]=="intern")
trunit <- read.delim("trunit.txt", as.is=T, header=T)

shinyServer(function(input, output, session) {


observe({ 
	if (input$geneserch==0){ return(NULL)}
	input$geneserch
	isolate({
	gone<-subset(kar,kar$id==input$gene_ser)
	updateNumericInput(session=session, inputId="instart", value=as.numeric(gone$start)-4000)
    updateNumericInput(session=session, inputId="inend", value=as.numeric(gone$end)+4000)
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
	updateSliderInput(session=session, inputId="inslider2", value=input$inslider2-1)

	}) 
  })
      observe({ 
	if (input$pushorfr==0){ return(NULL)}
	input$pushorfr
	isolate({
	updateSliderInput(session=session, inputId="inslider2", value=input$inslider2+1)

	}) 
  })

  observe({ 
	if (input$regsel==0){ return(NULL)}
	if(input$regsel=="twor")
	{
	updateTextInput(session=session, inputId="instart", value=seg[input$inslider2,2])
    updateTextInput(session=session, inputId="inend", value=seg[input$inslider2,3])
	
	} else{
	updateTextInput(session=session, inputId="instart", value=1000)
    updateTextInput(session=session, inputId="inend", value=10000)
	}
  }) 
  
  observe({ 
	if (input$locusserch==0){ return(NULL)}
	input$locusserch
	isolate({
	updateNumericInput(session=session, inputId="instart", value=1000)
    updateNumericInput(session=session, inputId="inend", value=10000)
	}) 
  }) 
#-----------------------------------------------------navigator--------------------------------------------------------------   
  
observe({ 
	if (input$twoleft==0){ return(NULL)}
	input$twoleft
	isolate({
	updateNumericInput(session=session, inputId="instart", value=as.numeric(input$instart)-8000)
    updateNumericInput(session=session, inputId="inend", value=as.numeric(input$inend)-8000)
	}) 
  }) 
 observe({ 
	if (input$oneleft==0){ return(NULL)}
	input$oneleft
	isolate({
	updateNumericInput(session=session, inputId="instart", value=as.numeric(input$instart)-4000)
    updateNumericInput(session=session, inputId="inend", value=as.numeric(input$inend)-4000)
	}) 
  }) 
observe({ 
	if (input$tworight==0){ return(NULL)}
	input$tworight
	isolate({
	updateNumericInput(session=session, inputId="instart", value=as.numeric(input$instart)+8000)
    updateNumericInput(session=session, inputId="inend", value=as.numeric(input$inend)+8000)
	}) 
  })  
 observe({ 
	if (input$oneright==0){ return(NULL)}
	input$oneright
	isolate({
	updateNumericInput(session=session, inputId="instart", value=as.numeric(input$instart)+4000)
    updateNumericInput(session=session, inputId="inend", value=as.numeric(input$inend)+4000)
	}) 
  }) 
  observe({ 
	if (input$pushright==0){ return(NULL)}
	input$pushright
	isolate({
    updateNumericInput(session=session, inputId="inend", value=as.numeric(input$inend)+4000)
	}) 
  }) 
    observe({ 
	if (input$pushleft==0){ return(NULL)}
	input$pushleft
	isolate({
    updateNumericInput(session=session, inputId="inend", value=as.numeric(input$inend)-4000)
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
	if (input$inend<5000){
    updateButton(session, "pushleft",  disabled = TRUE)}else{updateButton(session, "pushleft", disabled = FALSE)}
  }) 

  
  observe({ 
    if (input$moleft==0){ return(NULL)}
	input$moleft
		isolate({
	 updateSliderInput(session=session, inputId="inslider", value=as.numeric(input$inslider)-50)

  })  })
  
    observe({ 
    if (input$moright==0){ return(NULL)}
	input$moright
	isolate({

	 updateSliderInput(session=session, inputId="inslider", value=as.numeric(input$inslider)+50)

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
    withProgress(session, min=1, max=12, expr={
      for(i in 1:12) {
        setProgress(message = 'Circo layout in progress',
                    detail = 'Almost done...',
                    value=i)
        Sys.sleep(0.1)

      }})

return(mplot)

})

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


distab<-reactive({
if (is.null(input$inslider)){ return(NULL)}
if (is.null(input$lend)){ return(NULL)}
kar<-plotdata()
kar1<-subset(kar,kar$start>input$instart & kar$end < input$inend)
kar1<-kar1[order(kar1[,3]), ]
genecluster(kar1,input$lend)
})

output$disttab<-renderTable({
if (is.null(input$inslider)){ return(NULL)}
if (is.null(input$lend)){ return(NULL)}
distab()

})

observe({
updateSelectInput(session=session, inputId="distt", choices=rownames(distab()))

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

  
gos<-reactive({
kar<-plotdata()
kar1<-subset(kar,kar$start>input$inslider[1] & kar$end < input$inslider[2])
mdat<-merge(sgo,kar1,by="id",sort=F)

})
 
  
output$go_map<-renderPlot({
if (is.null(input$inslider)){ return(NULL)}
if (is.null(marr())){ return(NULL)}
gdat<-gos()
gdat1<- gdat[,c(1,3)]
gdat2<- princomp(t(table(gdat1)),cor=T)

biplot(gdat2,cex =0.8,expand=1)
})
output$go_sum<-renderPrint({
if (is.null(input$inslider)){ return(NULL)}
if (is.null(marr())){ return(NULL)}
gdat<-gos()
gdat1<- gdat[,c(1,3)]
summary(princomp(t(table(gdat1)), cor = TRUE,scores = TRUE))
})
output$go_load<-renderPrint({
if (is.null(input$inslider)){ return(NULL)}
if (is.null(marr())){ return(NULL)}
gdat<-gos()
gdat1<- gdat[,c(1,3)]
loadings(princomp(t(table(gdat1)), cor = TRUE,scores = TRUE))
})

output$go_pie<-renderPlot({
if (is.null(input$inslider)){ return(NULL)}
if (is.null(marr())){ return(NULL)}
mdat<-gos()
mdat<-subset(mdat,mdat[,2]=="M")
mdat1<- mdat[,c(1,3)]
mdat2<- table(mdat1[,2])

pie(mdat2,labels=names(mdat2), radius = 0.5,cex = 0.8)

})







updateButtonGroup(session, "btngrp1", toggle = "checkbox", style = "info", size = "small")
createAlert(session, inputId = "alert_anchor", message = "MESSAGE GOES HERE", type = "info", dismiss = TRUE, block = FALSE, append = FALSE )

output$myImage <- renderPlot({
data1<-plotdata()
if (is.null(input$inslider)){ return(NULL)}
else if(input$inslider[2]-input$inslider[1]>80000){
createAlert(session, inputId = "alert_gene", message = "Try to select region with in 80000 range!", type = "warning", dismiss = TRUE, block = FALSE, append = FALSE )}else{
input1 <-input$inslider[1]
input2 <-input$inslider[2]
	  
	if (is.null(input$inCheckboxGroup)){cogclr<-NULL}else{cogclr<-input$inCheckboxGroup}
	if (as.character(input$gridr)=="TRUE"){grir<-"TRUE"}else{grir<-"FALSE"}
mplot<-plotter(data1,trunit,input1,input2,0.2,1.9,"grey","grey",5,cogclr,grir)
observe({
     data1<-subset(data1,data1$start>input1 & data1$start<input2)
  	updateSelectInput(session=session, inputId="gselect", choices=paste(data1$id,'----',data1$name, sep=""),selected=data1[1,])})
	
	    withProgress(session, min=1, max=6, expr={
      for(i in 1:6) {
        setProgress(message = 'Synbrowser layout in progress',
                    detail = 'Almost done...',
                    value=i)
        Sys.sleep(0.1)
      }})  
	  
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
  
  
  
 output$message <- renderText({
gene1<-gene1()
  ou1<-paste("<div id='genebox' style='left:550px; top:180px; border:#c1c1c1; padding:10px; border:1px;'><b>",gene1[,c(1)],"</b>:",gene1[,c(6)],"</div>")
  updateTextInput(session=session, inputId="geneid", value=gene1[,c(1)])
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
  
  
  
  
  
  
  
  
  
  
#-----------------------------------------------------plot end--------------------------------------------------------------   
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
updateTextInput(session=session,inputId="pattern_ser1", value=input$pattern_ser2)
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
	 data.frame(t(pwm))
	 
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
withProgress(session, min=1, max=12, expr={
      for(i in 1:12) {
        setProgress(message = 'Pattern search initiated',
                    detail = 'Calculation started...',
                    value=i)
        Sys.sleep(0.1)

      }})
mysearch <-vmatchPattern(pat, subject , fixed=input$ambig)
count_index <- countIndex(mysearch)
with_match <- which(count_index >= 1)
 
withProgress(session, min=1, max=12, expr={
      for(i in 1:12) {
        setProgress(message = 'Pattern search progressed',
                    detail = 'Almost done...',
                    value=i)
        Sys.sleep(0.1)

      }})

res1<-data.frame(unlist(mysearch))
res2<-data.frame(names=names(subject[with_match]),seq=as.character(subject[with_match]),pat=rep(as.character(pat),length(with_match)))
res<-merge(res1,res2,by="names")
}
else{
withProgress(session, min=1, max=12, expr={
      for(i in 1:12) {
        setProgress(message = 'Pattern search initiated',
                    detail = 'Calculation started...',
                    value=i)
        Sys.sleep(0.1)

      }})
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
withProgress(session, min=1, max=12, expr={
      for(i in 1:12) {
        setProgress(message = 'Pattern search progressed',
                    detail = 'Almost done...',
                    value=i)
        Sys.sleep(0.1)

      }})
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
	updateNumericInput(session=session, inputId="patt_genlen_f", value=input$inslider[1])
	updateNumericInput(session=session, inputId="patt_genlen_t", value=input$inslider[2])
	}) 
  })	
  
  output$sliders3 <- renderUI({
  sliderInput(inputId="mm1",label="",min=0,max=5,value=0,step=1)
  
  })
 
 

#-----------------------------------------------------sequence end-------------------------------------------------------------- 

#-----------------------------------------------------dataset start-------------------------------------------------------------- 

  
    uploaddata<-reactive({
	inFile <- input$file1

    if (is.null(inFile))
      return(NULL)
    
   dat1<-read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
   dat1
	
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
udat$id
})
})

u_data2<- reactive({ 
distab<-distab()
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
mglist()[,c(1,3,4,5,6)]

})
 observe({
if (input$merge_glist==0){ return(NULL)}
input$merge_glist
isolate({ 
 updateSelectInput(session, "inSelectg", choices = mglist()[,1],selected=mglist()[,1]) 
 })
 })
  
 observe({
if (input$basket1==0){ return(NULL)}
input$basket1
isolate({   
updateTabsetPanel(session, inputId="inTabset", selected = "User Analysis Basket")
  })})
  
  observe({
if (input$ranger1==0){ return(NULL)}
input$ranger1

isolate({   

dfr<-data.frame(id=as.character(u_data2()))
u_data2<-merge(dfr,kar,by="id")
updateNumericInput(session, inputId="range_start", value = as.numeric(min(u_data2$start)))
updateNumericInput(session, inputId="range_end", value = as.numeric(max(u_data2$end)))
withProgress(session, min=1, max=12, expr={
      for(i in 1:12) {
        setProgress(message = 'Gene sequence range added',
                    detail = 'Open Sequence Extract tab...',
                    value=i)
        Sys.sleep(0.1)

      }})
  })}) 

  #-----------------------------------------------------dataset end-------------------------------------------------------------- 
observe({
if(input$inSelect3=="userdata"){
if(is.null(mglist())){return(NULL)}else{
updateSelectInput(session, "inSelect2", choices =input$inSelectg)
}}else if(input$inSelect3=="genedata"){
intdata<-subset(kar, kar[,3]>input$instart & kar[,4]<input$inend )
idata<-data.frame(bait=intdata$id)
updateSelectInput(session, "inSelect2", choices = idata[,1], selected=idata[1,])}
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
colnames(karp)<-c("id","prey")
one<-merge(karp,kasc,by="id")
kk<-data.frame(id=one$id,clr=one$clr)
colnames(karp)<-c("bait","id")
two<-merge(karp,kasc,by="id")
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
    # Depending on input$input_type, we'll generate a different
    # UI component and send it to the client.
    switch(input$input_type,
	     
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

 	updateSelectInput(session=session, inputId="gselect1", choices=annodata()[,1])})    })
  
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

})
