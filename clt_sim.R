
source("co_laplace_program_clt.R")
############################
#CLT

setwd("/Users/Treemanbaby/Desktop/zihanpaper/co laplace transform of volatility/program")


res1 = repfun(N=3000, n=23400, u=1, v=1.5, alpha1=0.5, alpha2=0.5,  omega1=0.01, omega2=0.01, theta=1/3,  gam=0)
write.table(res1$estclt, file="example clt V est omega1=omega2=0.01 u1v15.txt")


res2 = repfun(N=3000, n=23400, u=1, v=1.5, alpha1=0.5, alpha2=0.5,  omega1=0.03, omega2=0.03, theta=1/3,  gam=0)
write.table(res2$estclt, file="example clt V est omega1=omega2=0.03 u1v15.txt")

res3 = repfun(N=3000, n=23400, u=1, v=1.5, alpha1=0.5, alpha2=0.5,  omega1=0.1, omega2=0.1, theta=1/3,  gam=0)
write.table(res3$estclt, file="example clt V est omega1=omega2=0.1 u1v15.txt")




res1 = read.table("example clt V est omega1=omega2=0.01 u1v15.txt", head=TRUE)
res2 = read.table("example clt V est omega1=omega2=0.03 u1v15.txt", head=TRUE)
res3 = read.table("example clt V est omega1=omega2=0.1 u1v15.txt", head=TRUE)


x11 = seq(min(res1[,1]), max(res1[,1]), length=500)
x12 = seq(min(res1[,1]), max(res1[,1]), length=50)
x21 = seq(min(res2[,1]), max(res2[,1]), length=500)
x22 = seq(min(res2[,1]), max(res2[,1]), length=50)
x31 = seq(min(res3[,1]), max(res3[,1]), length=500)
x32 = seq(min(res3[,1]), max(res3[,1]), length=50)


windows(heigh=8, width=6, points=12)
# png(file="qqdensity2.png", width=6, height=8, pointsize = 10, units = "in", res=300)  
# par(mfrow=c(4,2), mar=c(4,4,2,1), font=6) 
par(mfrow=c(3,2), mar=c(4,4,2,1), font=6) 

hist(res1[,1], freq=FALSE, col=rainbow(70), breaks=x12, border=FALSE, axes=FALSE, xlim=c(-4, 4), 
     ylim=c(0,0.75), main="", xlab="", ylab="Density of studentized Vn", font=6, cex.lab=1.2)
lines(x11, dnorm(x11, mean=0, sd=1), col=2, lwd=2)
axis(1, at=seq(-4,4,length=5), labels=seq(-4,4,length=5), font.main=6)
axis(2, at=seq(0, 0.75, length=4), labels=seq(0, 0.75, length=4))
# abline(v=0, lty=2, col=2)
box()
mtext(expression(paste(alpha, "= 0.01")), side = 3, line = 0.5, lwd=2,font=6, cex=1)

qqplot(qnorm(1:1000/1000, mean=0, sd=1), quantile(res1[,1], 1:1000/1000), pch="+", 
       xlab="Normal Quantile", ylab="Sample Quantile", col=4, xlim=c(-4, 4), ylim=c(-4, 4), cex=2 ,cex.lab=1.2)
abline(a=0, b=1, col=2)


hist(res2[,1], freq=FALSE, col=rainbow(70), breaks=x22, border=FALSE, axes=FALSE, xlim=c(-4, 4), 
     ylim=c(0,0.75), main="", xlab="", ylab="Density of studentized Vn", font=6, cex.lab=1.2)
lines(x21, dnorm(x21, mean=0, sd=1), col=2, lwd=2)
axis(1, at=seq(-4,4,length=5), labels=seq(-4,4,length=5), font.main=6)
axis(2, at=seq(0, 0.75, length=4), labels=seq(0, 0.75, length=4))
# abline(v=0, lty=2, col=2)
box()
mtext(expression(paste(omega==0.03)), side = 3, line = 0.5, lwd=2,font=6, cex=1)

qqplot(qnorm(1:1000/1000, mean=0, sd=1), quantile(res2[,1], 1:1000/1000), pch="+", 
       xlab="Normal Quantile", ylab="Sample Quantile", col=4, xlim=c(-4, 4), ylim=c(-4, 4), cex=2 , cex.lab=1.2)
abline(a=0, b=1, col=2)


hist(res3[,1], freq=FALSE, col=rainbow(70), breaks=x32, border=FALSE, axes=FALSE, xlim=c(-4, 4), 
     ylim=c(0,0.75), main="", xlab="", ylab="Density of studentized Vn", font=6, cex.lab=1.2)
lines(x31, dnorm(x31, mean=0, sd=1), col=2, lwd=2)
axis(1, at=seq(-4,4,length=5), labels=seq(-4,4,length=5), font.main=6)
axis(2, at=seq(0, 0.75, length=4), labels=seq(0, 0.75, length=4))
# abline(v=0, lty=2, col=2)
box()
mtext(expression(paste(omega==0.1)), side = 3, line = 0.5, lwd=2,font=6, cex=1)

qqplot(qnorm(1:1000/1000, mean=0, sd=1), quantile(res3[,1], 1:1000/1000), pch="+", 
       xlab="Normal Quantile", ylab="Sample Quantile", col=4, xlim=c(-4, 4), ylim=c(-4, 4), cex=2, cex.lab=1.2 )
abline(a=0, b=1, col=2)





#######
#CLT
res1 = repfun(N=3000, n=23400, u=0.5, v=0.5, alpha1=0.5, alpha2=0.5,  omega1=0.01, omega2=0.01, theta=1/3,  gam=0)
write.table(res1$estclt, file="example clt V est omega1=omega2=0.01 uv05.txt")


res2 = repfun(N=3000, n=23400, u=0.5, v=0.5, alpha1=0.5, alpha2=0.5,  omega1=0.03, omega2=0.03, theta=1/3,  gam=0)
write.table(res2$estclt, file="example clt V est omega1=omega2=0.03 uv05.txt")

res3 = repfun(N=3000, n=23400, u=0.5, v=0.5, alpha1=0.5, alpha2=0.5,  omega1=0.1, omega2=0.1, theta=1/3,  gam=0)
write.table(res3$estclt, file="example clt V est omega1=omega2=0.1 uv05.txt")


res1 = read.table("example clt V est omega1=omega2=0.01 uv05.txt", head=TRUE)
res2 = read.table("example clt V est omega1=omega2=0.03 uv05.txt", head=TRUE)
res3 = read.table("example clt V est omega1=omega2=0.1 uv05.txt", head=TRUE)


x11 = seq(min(res1[,1]), max(res1[,1]), length=500)
x12 = seq(min(res1[,1]), max(res1[,1]), length=50)
x21 = seq(min(res2[,1]), max(res2[,1]), length=500)
x22 = seq(min(res2[,1]), max(res2[,1]), length=50)
x31 = seq(min(res3[,1]), max(res3[,1]), length=500)
x32 = seq(min(res3[,1]), max(res3[,1]), length=50)


windows(heigh=8, width=6, points=12)
# png(file="qqdensity2.png", width=6, height=8, pointsize = 10, units = "in", res=300)  
# par(mfrow=c(4,2), mar=c(4,4,2,1), font=6) 
par(mfrow=c(3,2), mar=c(4,4,2,1), font=6) 

hist(res1[,1], freq=FALSE, col=rainbow(70), breaks=x12, border=FALSE, axes=FALSE, xlim=c(-4, 4), 
     ylim=c(0,0.75), main="", xlab="", ylab="Density of co-Laplace transform of volatility", font=6, cex.lab=1.2)
lines(x11, dnorm(x11, mean=0, sd=1), col=2, lwd=2)
axis(1, at=seq(-4,4,length=5), labels=seq(-4,4,length=5), font.main=6)
axis(2, at=seq(0, 0.75, length=4), labels=seq(0, 0.75, length=4))
# abline(v=0, lty=2, col=2)
box()
mtext(expression(paste(omega==0.01)), side = 3, line = 0.5, lwd=2,font=6, cex=1)

qqplot(qnorm(1:1000/1000, mean=0, sd=1), quantile(res1[,1], 1:1000/1000), pch="+", 
       xlab="Normal Quantile", ylab="Sample Quantile", col=4, xlim=c(-4, 4), ylim=c(-4, 4), cex=2 ,cex.lab=1.2)
abline(a=0, b=1, col=2)


hist(res2[,1], freq=FALSE, col=rainbow(70), breaks=x22, border=FALSE, axes=FALSE, xlim=c(-4, 4), 
     ylim=c(0,0.75), main="", xlab="", ylab="Density of PRLT", font=6, cex.lab=1.2)
lines(x21, dnorm(x21, mean=0, sd=1), col=2, lwd=2)
axis(1, at=seq(-4,4,length=5), labels=seq(-4,4,length=5), font.main=6)
axis(2, at=seq(0, 0.75, length=4), labels=seq(0, 0.75, length=4))
# abline(v=0, lty=2, col=2)
box()
mtext(expression(paste(omega==0.03)), side = 3, line = 0.5, lwd=2,font=6, cex=1)

qqplot(qnorm(1:1000/1000, mean=0, sd=1), quantile(res2[,1], 1:1000/1000), pch="+", 
       xlab="Normal Quantile", ylab="Sample Quantile", col=4, xlim=c(-4, 4), ylim=c(-4, 4), cex=2 , cex.lab=1.2)
abline(a=0, b=1, col=2)


hist(res3[,1], freq=FALSE, col=rainbow(70), breaks=x32, border=FALSE, axes=FALSE, xlim=c(-4, 4), 
     ylim=c(0,0.75), main="", xlab="", ylab="Density of PRLT", font=6, cex.lab=1.2)
lines(x31, dnorm(x31, mean=0, sd=1), col=2, lwd=2)
axis(1, at=seq(-4,4,length=5), labels=seq(-4,4,length=5), font.main=6)
axis(2, at=seq(0, 0.75, length=4), labels=seq(0, 0.75, length=4))
# abline(v=0, lty=2, col=2)
box()
mtext(expression(paste(omega==0.1)), side = 3, line = 0.5, lwd=2,font=6, cex=1)

qqplot(qnorm(1:1000/1000, mean=0, sd=1), quantile(res3[,1], 1:1000/1000), pch="+", 
       xlab="Normal Quantile", ylab="Sample Quantile", col=4, xlim=c(-4, 4), ylim=c(-4, 4), cex=2, cex.lab=1.2 )
abline(a=0, b=1, col=2)









############################
#consistency
res111 = repfun(N=100, n=23400, u=0.2, v=0.2, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,  rho=0.3, gam=0.3)

res121 = repfun(N=5000, n=23400, u=0.2,  alpha=0.5, omega=0.03, theta1=1/3, theta2=1/3, rho=0)
res131 = repfun(N=5000, n=23400, u=0.2,  alpha=0.5, omega=0.05, theta1=1/3, theta2=1/3, rho=0)

res.omegac1 = rbind(res111$output, res121$output, res131$output)
write.table(res.omegac1, "example results for omega theta c1.txt")


##test theta take value 1
res1 = repfun(N=1000, n=23400, u=0.1, v=0.1, alpha1=0.5, alpha2=0.5,  omega1=0.01, omega2=0.01, theta=1,  gam=0)
write.table(res1$estclt, file="example clt V est omega1=omega2=0.01 uv0.1 theta1.txt")


res2 = repfun(N=1000, n=23400, u=0.1, v=0.1, alpha1=0.5, alpha2=0.5,  omega1=0.03, omega2=0.03, theta=1,  gam=0)
write.table(res2$estclt, file="example clt V est omega1=omega2=0.03 uv0.1 theta1.txt")

res3 = repfun(N=1000, n=23400, u=0.1, v=0.1, alpha1=0.5, alpha2=0.5,  omega1=0.1, omega2=0.1, theta=1,  gam=0)
write.table(res3$estclt, file="example clt V est omega1=omega2=0.1 uv0.1 theta1.txt")

res1 = read.table("example clt V est omega1=omega2=0.01 uv0.1 theta1.txt", head=TRUE)



x11 = seq(min(res1[,1]), max(res1[,1]), length=500)
x12 = seq(min(res1[,1]), max(res1[,1]), length=50)


windows(heigh=8, width=6, points=12)
# png(file="qqdensity2.png", width=6, height=8, pointsize = 10, units = "in", res=300)  
# par(mfrow=c(4,2), mar=c(4,4,2,1), font=6) 
par(mfrow=c(3,2), mar=c(4,4,2,1), font=6) 

hist(res1[,1], freq=FALSE, col=rainbow(70), breaks=x12, border=FALSE, axes=FALSE, xlim=c(-4, 4), 
     ylim=c(0,0.75), main="", xlab="", ylab="Density of co-Laplace transform of volatility", font=6, cex.lab=1.2)
lines(x11, dnorm(x11, mean=0, sd=1), col=2, lwd=2)
axis(1, at=seq(-4,4,length=5), labels=seq(-4,4,length=5), font.main=6)
axis(2, at=seq(0, 0.75, length=4), labels=seq(0, 0.75, length=4))
# abline(v=0, lty=2, col=2)
box()
mtext(expression(omega==0.01), side = 3, line = 0.5, lwd=2,font=6, cex=1)

qqplot(qnorm(1:1000/1000, mean=0, sd=1), quantile(res1[,1], 1:1000/1000), pch="+", 
       xlab="Normal Quantile", ylab="Sample Quantile", col=4, xlim=c(-4, 4), ylim=c(-4, 4), cex=2 ,cex.lab=1.2)
abline(a=0, b=1, col=2)

