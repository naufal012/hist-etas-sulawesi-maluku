#Preparation for environments
source("drawmap.R")
source("f2.R")
#
#num<-scan('workb.config',nline=1); eqnum <- num[1]; lonw<- num[4]; latw<-num[5]
 num<- scan('delone2.out',nline=1); eqnum <- num[1]; lonw<- num[4]; latw<-num[5]
#eqs<- matrix(scan('workb.config',skip=1,nlines=eqnum),ncol=8,byrow=T)
 eqs<- matrix(scan('delone2.out',skip=1,nlines=eqnum),ncol=8,byrow=T)
#eq.val <- matrix(scan('alljpM50nBvalues.out', skip=2), ncol=1)
eq.val <- matrix(scan('delo2d-bvalues.omap', skip=1), ncol=1)
#RefVal <- 1 #scan('delopara.test', skip=3, nlines=2)[c(1,2,4,5,7)]
#RefVal <- 1 #scan('delo2d-bvalues.omap', skip=3, nlines=2)[c(1,2,4,5,7)]
eq.val[,1]<-exp(eq.val[,1])
#eq.val[,3]<-exp(eq.val[,3])*.3

#postscript("colorpoint.ps",paper='special',horizontal=T, height=6.0, width=6.9)
pdf("delo2d-bvalues.pdf",paper='special',height=6.0, width=6.9)

for(i in 1:1){

temp <- ecdf(eq.val[,i])(eq.val[,i])
zr1 <- range(pretty(temp))

zlevels <- 40 

col <- rev(rainbow(50)[1:40])
#col <- rainbow(50)[1:40]

lon<-eqs[,2]; lat<-eqs[,3]; lon0<-118.0; lat0<- -6.0
#lon1=141.5; lon2=144.5; lat1=39; lat2=41
 lon1=lon0; lon2=lon0+14.0; lat1=lat0; lat2=lat0+12.0

#filled2.contour(seq(lon0,lon0+lonw, length=10), seq(lat0,lat0+latw, length=10),
filled2.contour(seq(lon1,lon2, length=10), seq(lat1,lat2, length=10),
                matrix(-1e-10,ncol=10, nrow=10),
      #This value is useless is make it blank
                col=c(rgb(1,1,1),col), 
                levels =c(-1e-9,seq(0,1, 1/zlevels)),
                other.comm={
                    drawmap.japan(F)
                    drawmap.japan(T,lwd=0.25,col='gray')
                     points(lon+lon0, lat+lat0, pch=16, 
                            col=col[(temp-zr1[1])/(zr1[2]-zr1[1])*zlevels+1],
                            cex=0.75)
#                           cex=1.0*(eqs[,4]-5.5))
                },

               key.axes={
                  axis(4, at= ecdf(eq.val[,i])(pretty(eq.val[,i], n=10)),
                       label=(pretty(eq.val[,i], n=10)))
                 },

             key.frame=F
 )
}

dev.off()
