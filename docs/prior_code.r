layout(matrix(1:6,nrow=2))
par(mgp=c(1.5,.4,0),tcl=.2,mar=c(1.5,2.5,2,.1))

plot(c(0,20),c(0,.2),xlab="",ylab="Prob",main="Uniform: min, max",type="n")
lines(0:100/5,dunif(0:100/5,6,14),col="red",lwd=2)
text(10,.14,"{U:6,14}",col="red")
lines(0:100/5,dunif(0:100/5,2,18),col="blue",lwd=2)
text(18,.08,"{U:2,18}",col="blue")

plot(c(0,20),c(0,.2),xlab="",ylab="Prob",main="Normal: mean, sd",type="n")
lines(0:100/5,dnorm(0:100/5,5,2),col="red",lwd=2)
text(10,.14,"{N:5,2}",col="red")
lines(0:100/5,dnorm(0:100/5,12,7),col="blue",lwd=2)
text(18,.07,"{U:12,6}",col="blue")

plot(c(0,20),c(0,.2),xlab="",ylab="Prob",main=expression(paste("Gamma: k, ",theta)),type="n")
lines(0:100/5,dgamma(0:100/5,shape=1,scale=5),col="red",lwd=2)
text(5,.14,"{G:1,5}",col="red")
lines(0:100/5,dgamma(0:100/5,shape=4,scale=2),col="blue",lwd=2)
text(14,.08,"{G:4,2}",col="blue")

plot(c(0,20),c(0,.2),xlab="",ylab="Prob",main="Exponential: scale",type="n")
lines(0:100/5,dexp(0:100/5,1/.5),col="red",lwd=2)
text(5,.14,"{E:0.5}",col="red")
lines(0:100/5,dexp(0:100/5,1/10),col="blue",lwd=2)
text(15,.05,"{E:10}",col="blue")

plot(c(0,20),c(0,.2),xlab="",ylab="Prob",main="Geometric: p",type="n")
lines(0:100/5,(.5)^(floor(0:100/5))*.5,col="red",lwd=2)
text(7,.14,"{M:0.5}",col="red")
lines(0:100/5,(.92)^(floor(0:100/5))*.08,col="blue",lwd=2)
text(18,.06,"{M:0.08}",col="blue")

