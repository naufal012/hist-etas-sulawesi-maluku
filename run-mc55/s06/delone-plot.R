#   plots the Delaunay triangles
 source("drawmap.R")
 source("f2.R")
#
pdf(file="delone-plot.pdf", height=6, width=6)
infilename <- "delone2.out"

dataNumbers <-scan(infilename, nlines=1)
np <- dataNumbers[2]
nt <- dataNumbers[3]
dataPoint <-matrix(scan(infilename, skip=1, nlines=np),ncol=8, byrow=T)
delTriangle <-scan(infilename, skip = np +1, nlines=nt, what=list(SeqNo=0,
   Id1=0, Id2=0, Id3=0, ''))
delTriangle <- cbind( delTriangle$SeqNo, delTriangle$Id1,
                      delTriangle$Id2, delTriangle$Id3) 
dataPoint[,2]=dataPoint[,2]+118.0    #128. is map origin of lonitude
dataPoint[,3]=dataPoint[,3]+-6.0     #30. is map origin of latitude

#   remove rubbish axis labels
plot(dataPoint[,2:3],pch='.', cex=(dataPoint[,4]-2)/3, yaxs="i", xaxs="i",
     xlab="", ylab="", axes=FALSE)
                    drawmap.japan(lwd=2.5)
# axis(1, at=seq(0, 18, 2), lab=seq(164, 182, 2))
# axis(2, at=seq(0, 12, 2), lab=seq(-48, -36, 2))
axis(1)
axis(2)
box()

par(lwd=0.5)
for(i in 1:nrow(delTriangle))
    polygon(dataPoint[delTriangle[i,2:4],2], dataPoint[delTriangle[i,2:4],3])

dev.off()

