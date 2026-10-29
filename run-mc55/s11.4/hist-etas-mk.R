#Preparation for environments
source("drawmap.R")
source("f2.R")

#read data
 eq.file <- "delone2.out"
#eqs <- read.table(eq.file, skip=1)
#
 num<- scan(eq.file,nline=1); eqnum <- num[1]; lonw<- num[4]; latw<-num[5]
eqs<- matrix(scan(eq.file,skip=1,nlines=eqnum),ncol=8,byrow=T)
#
eq.val <- matrix(scan('hist-etas-mk.omap', skip=3), ncol=2)
RefVal <- scan('hist-etas-mk.omap', skip=1, nlines=2)[c(1,2,4,5,7)]
 
pdf("hist-etas-mk.pdf",paper='special',height=6.5, width=8)

for(i in 1:2){

temp <- ecdf(eq.val[,i])(eq.val[,i])

zr1 <- (range(pretty(temp)))

zlevels <- 40 

col <- rev(rainbow(50)[1:40])

lon0<- 118.0; lat0<- -6.0
lon<-eqs[,2]+lon0
lat<-eqs[,3]+lat0
filled2.contour(seq(lon0,lon0+lonw, length=10),
                seq(lat0,lat0+latw, length=10),
                matrix(-1e-10,ncol=10, nrow=10),
                col=c(rgb(1,1,1),col), 
                levels =c(-1e-9,seq(0,1, 1/zlevels)),
                other.comm={
#                    drawmap.calcoast(lwd=2, col=gray(0.5))
                    #drawmap.world()
                     drawmap.japan()
#                    drawmap.calfault(col='black', lwd=0.5)
#points(eqs[,2], eqs[,3], pch='.', 
#                    points(lon, lat, pch='.', 
                     points(lon, lat, pch=16, 
                            col=col[(temp-zr1[1])/(zr1[2]-zr1[1])*zlevels+1],
                            cex=1.25)
                },
               key.axes={
                  # axis(4, at= ecdf(eq.val[,i])(pretty(eq.val[,i], n=10)),
                  #       label=(pretty(eq.val[,i], n=10)))

                 if(i>2) axis(4, at= ecdf(eq.val[,i])(log(pretty(exp(eq.val[,i])*RefVal[i], n=10)/RefVal[i])),
                 label=(pretty(exp(eq.val[,i])*RefVal[i], n=10)))

                if(i<=2)axis(4, at=ecdf(eq.val[,i])(seq(-10,10, 1)*log(10)-log(RefVal[i])), 
                 label=c(expression(10^{-10}), expression(10^{-9}), 
                         expression(10^{-8}), 
                         expression(10^{-7}), expression(10^{-6}), 
                         expression(10^{-5}), expression(10^{-4}), 
                         expression(10^{-3}), expression(10^{-2}),
                         0.1, 1, 10, expression(10^2), 
                         expression(10^3),  expression(10^4), 
                         expression(10^5),  expression(10^6),  
                         expression(10^7),  expression(10^8),  
                         expression(10^9),  expression(10^10))) 
                  },
             key.frame=F
 )
}

dev.off()
