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


subs <- function(x, n) {
  sst <- strsplit(x, '')[[1]]
  m <- matrix('', nrow=n, ncol=(length(sst)+n-1)%/%n)
  m[seq_along(sst)] <- sst
  apply(m, 2, paste, collapse='')
}



plotter<- function(data=data,trunit=trunit,rstart=rstart,rend=rend, genepos=genepos, labelpos=labelpos,gcol=gcol,lcol=lcol, size=size,genecol,grir){
  
kar<-data
if (is.null(rstart) & is.null(rend)){ return(NULL)}
else{
plot(x=1, y=1, type="h", col="white", xlim=c(as.numeric(rstart),as.numeric(rend)), ylim=c(-2,3.5), xlab="", ylab="", yaxt="n",xaxt="n",bty="n")
box(lwd=1,col="#cccccc")
if(grir=="TRUE"){
grid(nx = 10, ny = NULL, col = "cornsilk2", lty = 6 ,lwd = 0.5)}

axis(3, at=seq(as.numeric(rstart),as.numeric(rend),by=100),col='#cccccc');
axis(1, at=seq(as.numeric(rstart),as.numeric(rend),by=100),col='#cccccc');

kar1 <- subset(kar,kar[,5]=="+")
kar2 <- subset(kar,kar[,5]=="-")
clr <- c(rgb(1, 0, 0, 0.1), rgb(0, 0, 1, 0.1),rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5))


if(!is.null(trunit)){
tu1 <- subset(trunit,trunit[,5]=="+")
tu2 <- subset(trunit,trunit[,5]=="-")

if(as.numeric(rend)-as.numeric(rstart)<20000){
rect(tu1$start, 0, tu1$end, genepos+0.001, xpd=FALSE , border="#ffffff", col="#00FF0020")
rect(tu2$start, -genepos-0.001, tu2$end, 0, xpd=FALSE , border="#ffffff", col="#FF000020")
}}

if(genecol=="option1"){
for(i in 1:nrow(data)){
if(data[i,5]=="+"){
polygon(c(data[i,4]-100,data[i,3],data[i,3],data[i,4]-100,data[i,4]),c(0.1,0.1,0.9,0.9,0.5),col =data[i,8],border = "#cccccc")
if(data[i,4]-data[i,3]>=2){
text((data[i,4]+data[i,3])/2,1.5,labels=data[i,1],cex = 1, srt = 90,col = "#666666")
}else{NULL}}
else{
polygon(c(data[i,3],data[i,3]+100,data[i,4],data[i,4],data[i,3]+100),c(-0.5,-0.1,-0.1,-0.9,-0.9),col =data[i,8],border = "#cccccc")
if(data[i,4]-data[i,3]>=500){
text((data[i,4]+data[i,3])/2,-1.5,labels=data[i,1],cex = 1, srt = 90,col = "#666666")
}else{text((data[i,4]+data[i,3])/2,-1.5,labels="",cex = 1, srt = 90,col = "#666666")}}

}

}else if(genecol=="option2"){
for(i in 1:nrow(data)){
if(data[i,5]=="+"){
polygon(c(data[i,4]-100,data[i,3],data[i,3],data[i,4]-100,data[i,4]),c(0.1,0.1,0.9,0.9,0.5),col =data[i,8],border = "#cccccc")
if(data[i,4]-data[i,3]>=2){
text((data[i,4]+data[i,3])/2,1.5,labels=data[i,1],cex = 1, srt = 90,col = "#666666")
}else{NULL}}
else{
polygon(c(data[i,3],data[i,3]+100,data[i,4],data[i,4],data[i,3]+100),c(-0.5,-0.1,-0.1,-0.9,-0.9),col =data[i,8],border = "#cccccc")
if(data[i,4]-data[i,3]>=500){
text((data[i,4]+data[i,3])/2,-1.5,labels=data[i,1],cex = 1, srt = 90,col = "#666666")
}else{text((data[i,4]+data[i,3])/2,-1.5,labels="",cex = 1, srt = 90,col = "#666666")}}

}
}

else{
for(i in 1:nrow(data)){
if(data[i,5]=="+"){
polygon(c(data[i,4]-100,data[i,3],data[i,3],data[i,4]-100,data[i,4]),c(0.1,0.1,0.9,0.9,0.5),col ="#f1f1f1",border = "#cccccc")
if(data[i,4]-data[i,3]>=2){
text((data[i,4]+data[i,3])/2,1.5,labels=data[i,1],cex = 1, srt = 90,col = "#666666")
}else{NULL}}
else{
polygon(c(data[i,3],data[i,3]+100,data[i,4],data[i,4],data[i,3]+100),c(-0.5,-0.1,-0.1,-0.9,-0.9),col ="#f1f1f1",border = "#cccccc")
if(data[i,4]-data[i,3]>=500){
text((data[i,4]+data[i,3])/2,-1.5,labels=data[i,1],cex = 1, srt = 90,col = "#666666")
}else{text((data[i,4]+data[i,3])/2,-1.5,labels="",cex = 1, srt = 90,col = "#666666")}}

}
}



segments(rstart, 0.0, rend, 0.0, col="#333333", lwd=1)

require(Biostrings, quietly = TRUE)
fs <- readDNAStringSet("NC_000911.fna", "fasta")

		glist<-subset(kar,kar$start>as.numeric(rstart) & kar$end <as.numeric(rend))
		frag <- Views(fs[[1]], as.numeric(rstart) ,as.numeric(rend))
		seq_data1<-as.character(frag)
if(as.numeric(rend)-as.numeric(rstart)==200){
text(as.numeric(rstart):as.numeric(rend),2, labels=unlist(strsplit(seq_data1,'')),cex = 0.8)
segments(kar1$start,2,kar1$end,2, col = clr[1], lwd = 20)
segments(kar2$start,2,kar2$end,2, col = clr[2], lwd = 20)
}		else if(as.numeric(rend)-as.numeric(rstart)==150){
text(as.numeric(rstart):as.numeric(rend),2, labels=unlist(strsplit(seq_data1,'')),cex = 1)
segments(kar1$start,2,kar1$end,2, col = clr[1], lwd = 22)
segments(kar2$start,2,kar2$end,2, col = clr[2], lwd = 22)
}
		else if(as.numeric(rend)-as.numeric(rstart)==100){
text(as.numeric(rstart):as.numeric(rend),2, labels=unlist(strsplit(seq_data1,'')),cex = 1.4)
segments(kar1$start,2,kar1$end,2, col = clr[1], lwd = 24)
segments(kar2$start,2,kar2$end,2, col = clr[2], lwd = 24)}


mtext("\n\n\n\n\n+\n\n\n\n-", side=2, line=1, cex.lab=1.2,las=1, col="black")
mtext("\n\n\n\n\n+\n\n\n\n-", side=4, line=1, cex.lab=1.2,las=1, col="black")
#title(paste("Synechocystis Genome view"))
#legend("topright",legend=c("Gene","TATA box","asRNA","ncRNA","Promoter","5'UTR"),lty=1,lwd=2,col=c("Grey","red","blue","black","green","#FF6633"),ncol=2,bty="n",cex=1.2,text.col=c("black"),inset=0.01)


}
}
  
track<- function(trk=trk,rstart=rstart,rend=rend,size=size,grir){
  
if (is.null(rstart) & is.null(rend)){ return(NULL)}

plot(x=1, y=1, type="h", col="white", xlim=c(as.numeric(rstart),as.numeric(rend)), ylim=c(-1,1), xlab="", ylab="", yaxt="n",xaxt="n",bty="n")
if(grir=="TRUE"){
grid(nx = 10, ny = NA, col = "cornsilk2", lty = 6 ,lwd = 0.5)}
axis(1, at=seq(as.numeric(rstart),as.numeric(rend),by=100),col='#C1C1C1');

if(!is.null(trk)){
trkd<-trk
trk1 <- subset(trkd,trkd[,5]=="+")
trk2 <- subset(trkd,trkd[,5]=="-")

cdf<-data.frame(type=c("TATA","ncRNA","asRNA","intern","5'UTR","TP"),colr=c("red","blue","black","green","#FF6633","pink"))
cd1<-subset(cdf,cdf[,1]==trk1[1,1])
cd1<-as.character(cd1$colr)

rect(rstart-1000, -0.8, rend+1000, 0.8, col="#ededed",xpd=FALSE , border="#cccccc")
segments(rstart, 0, rend, 0, col="grey", lwd=1)
segments(trk1$start, 0.5, trk1$end, 0.5, col=as.character(cd1),lwd=4,lend=1)
segments(trk2$start, -0.5, trk2$end, -0.5, col=as.character(cd1),lwd=4,lend=1)
if(as.numeric(rend)-as.numeric(rstart)<5000){
text((trk1$end)+100,0.5,labels=trk1$name,cex = 1, srt =0,col = "#666666")
text((trk2$end)+100,-0.5,labels=trk1$name,cex = 1, srt =0,col = "#666666")}

mtext("Track", side=2, line=1, cex.lab=1,las=1, col="black")

}


} 
  
  

gc_content<-function(data=data,rstart=rstart,rend=rend,grir,wpane){
kar<-data
require(Biostrings, quietly = TRUE)
fs <- readDNAStringSet("NC_000911.fna", "fasta")

wpane <- wpane
		glist<-subset(kar,kar$start>as.numeric(rstart) & kar$end <as.numeric(rend))
		frag <- Views(fs[[1]], as.numeric(rstart) ,as.numeric(rend))
		seq_data1<-as.character(frag)

gc <- rowSums(letterFrequencyInSlidingView(DNAString(seq_data1), wpane, c("G","C")))/wpane

plot(gc,typ='l', ann=F, bty="n",col="#009999")
if(grir=="TRUE"){
grid(nx = 10, ny = NA, col = "cornsilk2", lty = 6 ,lwd = 0.5)}
mtext("GC content", side=3, adj=0, line=1.2, cex=1, font=2)
box(lwd=1,col="#f1f1f1")
lines(lowess(x=gc), col=2)

} 
   
  
  
circo_plot<-function(data=data,seg=seg,cin=cin,kasc,cogc,trunit){
if (is.null(seg) & is.null(cin)){ return(NULL)}
da1<-data[,c(1,3,4,5)]
sect_sel<-subset(seg,seg$V1==cin)
da1<-subset(da1,da1$start>sect_sel[1,2] & da1$end<sect_sel[1,3])
da2<-subset(kasc,kasc$start>sect_sel[1,2] & kasc$end<sect_sel[1,3])
da3<-subset(cogc,cogc$start>sect_sel[1,2] & cogc$end<sect_sel[1,3])
da4<-subset(trunit,trunit$start>sect_sel[1,2] & trunit$end<sect_sel[1,3])

sect<-subset(seg,seg[,2]==sect_sel[1,2] & seg[,3]==sect_sel[1,3])
circos.clear()
par(mar = c(0.5, 0.5, 0.5, 0.5))
circos.par(points.overflow.warning = FALSE)
circos.par("canvas.xlim" = c(-1,1), "canvas.ylim" = c(-1, 1),gap.degree =1,start.degree = 0)
factors = seg[,1]
factors<-ordered(factors, levels = c(paste("S",1:72,sep="")))
circos.initialize(factors = factors, xlim = c(0,1))
circos.trackPlotRegion(factors = factors, ylim = c(0,1), bg.border = "GREY",bg.col = "#f1f1f1",track.height = 0.1,panel.fun = function(x, y,...) {
   #select details of current sector
   name = get.cell.meta.data("sector.index")
   i = get.cell.meta.data("sector.numeric.index")
   xlim = get.cell.meta.data("xlim")
   #plot labels
   circos.text(x=mean(xlim), y=0.5, labels=name, facing = "outside", cex=1.5,col="#cccccc")
})
circos.updatePlotRegion(sector.index = sect[1,1], bg.col = "RED",bg.border = "GREY")
circos.text(0.5, 0.5, sect[1,1],cex = 1.5,facing = "outside", col="white")
#circos.link("S1",1,"S12",1,col="#f1f1f1", border = "white")
#circos.link("S1",1,"S42",1,col="#f1f1f1", border = "white")

circos.clear()
if(!is.null(data)){
par(new = TRUE)
circos.par("canvas.xlim" = c(-1.4,1.4), "canvas.ylim" = c(-1.4, 1.4))
circos.genomicInitialize(da1,sector.names = NULL,plotType = c("label"))

ggcol<-list()
for(i in 1:nrow(da1)){
if(da1[i,4]=="+"){ggcol[i]<-"lightblue"}else{ggcol[i]<-"#F5DEB3"}
}
gcol<-do.call(rbind,ggcol)

circos.genomicTrackPlotRegion(da1, ylim = c(0,1),bg.border = "grey",bg.col = gcol, track.height = 0.1,panel.fun = function(x, y,...) {
   #select details of current sector
   name = get.cell.meta.data("sector.index")
   i = get.cell.meta.data("sector.numeric.index")
   xlim = get.cell.meta.data("xlim")
   #plot labels
   circos.text(x=mean(xlim), y=0, labels=name, facing = "outside", cex=1.5)
circos.axis(sector.index = name,labels.facing = "clockwise",labels.cex = 1.2, major.tick.percentage = 0.2, labels.away.percentage = 0.3, minor.ticks = 0)
})

text(0, 0, paste("Coordinate\n",sect[,2],"-",sect[,3],"bp" ), cex = 1.5)
col = c("#FF000020", "#00FF0020", "#0000FF20")

#circos.genomicLink(da1[c(1,6,8),],da1[c(11,12,20),],col="#0000FF20", border = "white")
}
circos.clear()



if(!is.null(cogc)){
par(new = TRUE)
circos.par("canvas.xlim" = c(-1.7,1.7), "canvas.ylim" = c(-1.7, 1.7))
circos.genomicInitialize(da1[,c(1,2,3)],sector.names = NULL,plotType = c("label"))
circos.genomicTrackPlotRegion(da1[,c(1,2,3)], ylim = c(0,0.2),bg.border = "grey",bg.col = paste(data[,8],"80",sep=""), track.height = 0.05,panel.fun = function(x, y,...) {
   #select details of current sector
   name = get.cell.meta.data("sector.index")
   i = get.cell.meta.data("sector.numeric.index")
   xlim = get.cell.meta.data("xlim")
})
}

if(!is.null(kasc)){
circos.clear()
par(new = TRUE)
circos.par("canvas.xlim" = c(-1.9,1.9), "canvas.ylim" = c(-1.9, 1.9))
circos.genomicInitialize(da1[,c(1,2,3)],sector.names = NULL,plotType = c("label"))
circos.genomicTrackPlotRegion(da1[,c(1,2,3)], ylim = c(0,0.2),bg.border = "grey",bg.col = paste(da2[,8],"80",sep=""), track.height = 0.05,panel.fun = function(x, y,...) {
   #select details of current sector
   name = get.cell.meta.data("sector.index")
   i = get.cell.meta.data("sector.numeric.index")
   xlim = get.cell.meta.data("xlim")
   #circos.text(x=mean(xlim), y=-1, labels=da3[i,7], facing = "clockwise", cex=1.2)
})

}


circos.clear()


}  
  
  
  
  
  
  
  
  
genecluster<-function(data=data, len=len){ 
kar<-data

karpp<-subset(kar,kar[,5]=="+")
karmm <-subset(kar,kar[,5]=="-")

plusg <- data.frame(karpp[1:nrow(karpp)-1,1],karpp[2:nrow(karpp),1],karpp[2:nrow(karpp),3]-karpp[1:nrow(karpp)-1,4])
minusg <- data.frame(karmm[1:nrow(karmm)-1,1],karmm[2:nrow(karmm),1],karmm[2:nrow(karmm),3]-karmm[1:nrow(karmm)-1,4])

colnames(plusg)<- c("gene1","gene2","dis")
colnames(minusg)<- c("gene1","gene2","dis")

objectsp <- c(as.character(plusg[1, "gene1"]), as.character(plusg[, "gene2"]))
scoresp <- c(0, plusg[, "dis"]) 
plus1 <- split(objectsp, cumsum(scoresp > len))

objectsm <- c(as.character(minusg[1, "gene1"]), as.character(minusg[, "gene2"]))
scoresm <- c(0, minusg[, "dis"]) 
minus1 <- split(objectsm, cumsum(scoresm > len))
names(plus1) <- paste("Pset", 1:length(	plus1), sep='')
names(minus1) <- paste("Mset", 1:length(minus1), sep='')

com <- c(plus1,minus1)
 
n <- max(sapply(com, length))

res<-do.call(rbind, lapply(com, `[`, seq_len(n)))

unique(res)

}
