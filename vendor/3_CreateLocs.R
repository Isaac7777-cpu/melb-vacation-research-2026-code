setwd("/Users/tingjinc/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/A1_Spatial/1-RGeo/locs")

library(spatstat)

rate  = 0.5
n = 800                              # sample size
lbase = 10
l=lbase*(n/50)^rate                     # spatial domain of interest [-l/2,l/2]^2
lsize = (lbase/10)*(n/50)^(rate-1/2)    # 0.2*lsize: the minimum distance of two locations 

set.seed(2019+1000*rate)
loc = rSSI(0.2*lsize, n, win=as.owin(c(-l/2,l/2,-l/2,l/2)))


nameRdata=strwrap(paste("loc", n, ".Rdata", sep=""))
save(loc, file=nameRdata)





save(loc, file="loc100rate5.Rdata")

  
rm(list = ls(all = TRUE))
setwd("~/Dropbox/A2-1Mixed/1-Code")
#setwd("C:/Users/tingjinc/Dropbox/A2-1Mixed/1-Code")

library(spatstat)
##   Ratio
cov.model = "Exp"
n = 800
qq = 0.25
lbase = 10
lrate = 7
if (cov.model == "Exp"){
  sigmasq=4/2; tausq=1/2; r=3
  theta0 = c(r, tausq, sigmasq)
  dq = log(1.25/qq)*(theta0[1])
  theta.ini = c(2.4,0.3)
} else if (cov.model=="Gau"){
  sigmasq=4/2; tausq=1/2; r=5
  theta0 = c(r, tausq, sigmasq)
  dq = sqrt(log(1.25/qq))*(theta0[1])
  theta.ini = c(4.2, 0.3)
} else if (cov.model=="Cau"){
  theta0 = c(4, 1/2, 2)
  pkap = 0.8
  dq = theta0[1]*sqrt((1.25/qq)^(5/4)-1)
  theta.ini = c(3.5,0.25)
} else if (cov.model=="PExp"){
  theta0 = c(3, 1/2, 2, 1.5)
  theta.ini = c(2.4, 0.3, 1.2)
}
pr = dq
l = pr*lrate 
lsize = min(dq*3/8, (lbase/10)*(n/100)^(0.4-1/2)*lrate*0.8)  

set.seed(2019+1000*lrate)
loc = rSSI(0.2*lsize, n, win=as.owin(c(-l/2,l/2,-l/2,l/2)))

nameRdata=strwrap(paste("rloc/loc", cov.model, lrate, ".Rdata", sep=""))
save(loc, file=nameRdata)


