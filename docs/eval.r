#load these functions into R by typing: source("c:/docs/ssc/eval.r")
plot.log<-function(logfile) {
	#displays traces of the populations and # of lineages from logfiles
	read.table(logfile,head=T)->lg
	c(0,which(diff(lg$Sim)==1),nrow(lg))+1->grp
	ncol(lg)/2-1->npop
	par(ask=T,mar=c(3,3,2,3.2),mgp=c(2,.9,0))
	for(i in 1:(length(grp)-1)) {
		lg[grp[i]:(grp[i+1]-1),]->cur
		max(cur[,1:npop*2+1])->maxy
		maxy/max(cur[,1:npop*2+2])->scl
		plot(0,0,xlim=range(cur$Time),ylim=c(0,maxy),main=paste("Simulation",i),type="n",xlab="time",ylab="Ne")
		axis(4,at=0:(maxy/scl)*scl,labels=0:(maxy/scl))
		mtext("Lineages",4,1.7)
		for(j in 1:npop) {
			lines(cur$Time,cur[,j*2+1],lwd=3,col=rainbow(npop)[j])
			lines(cur$Time,cur[,j*2+2]*scl,lty="dashed",col=rainbow(npop)[j])
		}
	}
	par(ask=F)
}

reject<-function(statfile,cols=NA,obs=NA,Beaumont=F,accept.perc=NA,lmrm=F) {
#fast rejection: works only with strictly numeric values, and only if the first 10K prior values are representative of the whole set
require(locfit)
post<-function(x,cdf,q) {
  res=0
  for(i in 1:length(q)) {
    ix=sum(cdf<q[i])+0:1
    res[i]=ifelse(ix[1]==0,NA,weighted.mean(x[ix],w=abs(cdf[rev(ix)]-q[i])))
  }
  names(res)<-paste(q*100,"%",sep="")
  res
}
edgeless.density<-function(x,refl=2) {
	#refl: # of bandwidths to extend data
	range(x)->z
	2*bw.nrd0(sf[,pr[i]])->b #should be (x), right?
	which(x<z[1]+b)->ix.l
	which(x>z[2]-b)->ix.h
	density(c(2*z[1]-x[ix.l],x,2*z[2]-x[ix.h]))->r
	which(r$x<z[1])->ix.l
	ix.l[2:length(ix.l)-1]->ix.l
	which(r$x>z[2])[-1]->ix.h
	r$x[-c(ix.l,ix.h)]->r$x
	r$y[-c(ix.l,ix.h)]->r$y
	r
}

if(lmrm) read.csv(statfile,as.is=T,nrows=10000)->sf
else read.csv(statfile,as.is=T)->sf
for(i in 1:ncol(sf)) sf[,i]=as.numeric(sf[,i])
print(colnames(sf))
if(is.na(cols[1])) {
	readline("Which column/s (eg 4,23,27)? ")->cols
	as.numeric(strsplit(cols,',')[[1]])->cols
}
if(is.na(obs[1])) {
	readline("Observed values: ")->obs
	as.numeric(strsplit(obs,',')[[1]])->obs
}
euc=0; sds=0
for(i in 1:length(cols))
	euc=euc+((sf[,cols[i]]-obs[i])/sd(sf[,cols[i]],na.rm=T))^2
sqrt(euc)->euc
which(is.na(euc))->ix; if(length(ix)) {euc[-ix]->euc; sf[-ix,]->sf }
hist(log10(euc),xlab="Log euclidean distance",main="Accuracy of Simulations (pre-rejection)")
print("To retain the following % of your simulations, set delta to the following values")
print(quantile(euc,c(.001,.01,.05,.1,.25),na.rm=T))
if(is.na(accept.perc)) { as.numeric(readline("Please set delta (rejection) value: "))->delta
} else quantile(euc,accept.perc,na.rm=T)->delta

pdf(paste(substr(statfile,1,nchar(statfile)-8),"post.pdf",sep=""),w=11,h=8.5)
hist(log10(euc),xlab="Log euclidean distance",main="Accuracy of Simulations (pre-rejection)")
points(log10(delta),0,pch="*",col="red",cex=3)
legend("topright",paste("delta = ",signif(delta,5)),col="red",pch="*")

which(colnames(sf)=="PRIORS")->pr
pr=c((pr+1):ncol(sf))
prd=list(0)
for(i in 1:length(pr)) prd[[i]]=edgeless.density(sf[,pr[i]],2) #prior pdfs

if(lmrm) {
	sc=paste(cols,collapse=","); so=paste(obs,collapse=",")
	scall=paste("c:/Docs/SSC/LMRM/LMRM.exe",statfile,sc,so,delta)
	system(scall)
	read.csv(paste(substr(statfile,1,nchar(statfile)-8),"acc.csv",sep=""))->sf
} else sf[euc<delta,]->sf
if(nrow(sf)==0) stop("No simulations accepted (delta value too low)")
which(na.per.row(as.matrix(sf[,cols]))>0)->dr
if(length(dr)>0) {
	print(paste("Dropped",length(dr),"sims due to missing values."))
	sf[-dr,]->sf
}
euc=0
for(i in 1:length(cols)) 
	if((sd(sf[,cols[i]],na.rm=T)>0)->x) euc=euc+((sf[,cols[i]]-obs[i])/x)^2
sqrt(euc)->euc
pdd=list(0)

for(i in 1:length(pr)) {
prwt=approx(prd[[i]]$x,prd[[i]]$y,sf[,pr[i]])$y
prwt[which(is.na(prwt) | prwt<mean(prwt,na.rm=T)/1000)]=mean(prwt,na.rm=T)/1000 #Consider alternatives
if(Beaumont) { w=(1-(euc/max(euc))^2)/prwt } else w=1/prwt
density(sf[,pr[i]],weights=w/sum(w),cut=0)->fit
dx=10^(log10(diff(range(fit$x)))-3)
fit$x=c(min(fit$x)-dx,fit$x,max(fit$x)+dx)
fit$y=c(0,fit$y,0)
plot(prd[[i]],type="n",main=paste("Posterior",i,colnames(sf)[pr[i]]),xlab="Parameter value",ylim=c(0,max(fit$y)))
polygon(fit,col="plum")
lines(prd[[i]]$x,prd[[i]]$y,lty="dashed")
fn=cumsum(fit$y/sum(fit$y))
x=seq(min(fit$x),max(fit$x),length.out=length(fit$x))
sval=c(MLE=fit$x[which.max(fit$y)],post(x,fn,c(.5,.05,.95))) #MLE,etc
legend(x=ifelse(sval[1]<max(x),"topright","topleft"),leg=paste(c("MLE","50%","5%","95%"),"=",signif(sval,3)))
pdd[[i]]=cbind(theta=x,cdf=fn)
}
cat("Posterior graphics output to ",paste(substr(statfile,1,nchar(statfile)-8),"post.pdf",sep=""),"\n\n")
while(names(dev.cur())!="pdf") dev.next() 
dev.off()
return(list(accept.sims=sf[,c(cols,pr)],pdd=pdd))
}


mle<-function(z) {
post<-function(x,cdf,q) {
  res=0
  for(i in 1:length(q)) {
    ix=sum(cdf<q[i])+0:1
    res[i]=ifelse(ix[1]==0,NA,weighted.mean(x[ix],w=abs(cdf[rev(ix)]-q[i])))
  }
  names(res)<-paste(q*100,"%",sep="")
  res
}
	nrow(z$accept)->n
	for(i in 1:length(z$pdd)) {
		z$pdd[[i]]->x
		y=post(x[,1],x[,2],c(.5,.025,.975))
		x[-1,"cdf"]=diff(x[,"cdf"])
		if(i==1) res=c(MLE=x[which.max(x[,"cdf"]),1],y,n=n)
		else res=rbind(res,c(MLE=x[which.max(x[,"cdf"]),1],y,n=n))
	}
	rownames(res)<-colnames(z$accept)[ncol(z$accept)-length(z$pdd):1+1]
	res
}

aic.ssc<-function(statfile,cols,obs,params) {
	read.csv(statfile)->sf
	prob=1
	unitnorm<-function(x) x=(x-min(x,na.rm=T))/diff(range(x,na.rm=T))
	floor(sqrt(length(cols)))->flc
	layout(matrix(1:(flc*ceiling(length(cols)/flc)),nc=flc))
	par(mar=c(.2,0,0,0))
	for(i in 1:length(cols)) {
		density(unitnorm(sf[,cols[i]]),na.rm=T,from=0,to=1)->z
		nearest=which.min(abs(z$x-obs[i]))
		print(prob<-prob*z$y[nearest]/max(z$y))
		plot(z,xlab="",ylab="",main="",axes=F); legend("top",colnames(sf)[cols[i]],pch="",bty="n")
		if(z$x[nearest]<max(z$x)) points(z$x[nearest],z$y[nearest],pch="*",col="red",cex=4)
		else text(mean(z$x),mean(z$y),"[off]",col="blue")
	}
	2*params-2*log(prob)
}

mod.comp<-function(fs,cols,obs,p.acc=.01) {
#fs=_stat.csv files to compare, #cols=columns of summary stats to use, #obs=observed values, #p.acc=% of sims to accept in first file
sf=c()
for(i in 1:length(fs)) sf=rbind(sf,cbind(read.csv(fs[i]),i))
sds=apply(sf[,cols],2,sd)
euc=rep(0,nrow(sf))
for(i in 1:length(cols)) euc=euc+(sf[,cols[i]]-obs[i])^2/sds[i]
sort(euc)[nrow(sf)*p.acc]->d
#wt=ifelse(euc>d,0,1-euc/d)
wt=ifelse(euc>d,0,1)
res=c()
for(i in 1:length(fs)) res[i]=sum((sf[,ncol(sf)]==i)*wt)
names(res)<-fs
res/sum(res)
}

#An alternative method to convert priors into posteriors
SIR<-function(statfile,cols=NA,obs=NA) {
read.csv(statfile)->sf
if(is.na(cols[1])) {
print(colnames(sf))
readline("Which column/s (eg 4,23,27)? ")->cols
as.numeric(strsplit(cols,',')[[1]])->cols
}
if(is.na(obs[1])) {
readline("Observed values: ")->obs
as.numeric(strsplit(obs,',')[[1]])->obs
}
which(colnames(sf)=="PRIORS")->pr
pr=c((pr+1):ncol(sf))
euc=0
for(i in 1:length(cols)) 
	euc=euc+((sf[,cols[i]]-obs[i])/diff(range(sf[,cols[i]])))^2
sqrt(euc)->euc
kw=0; density(sf[,pr[i]],bw=diff(range(sf[,pr[i]]))/100)->gth
d=sd(euc)
for(i in 1:nrow(sf)) {
	prior.prob=1
	for(j in 1:length(cols)) prior.prob=prior.prob*gth$y[sum(gth$x<sf[i,cols[j]])]
	kw[i]=exp(-euc[i]/d)/prior.prob
}	
kw=kw/sum(kw)
sample(nrow(sf),5000,replace=T,prob=kw)->rs
par(mar=c(2,2,2,.1),mgp=c(3,1,0))
layout(matrix(1:length(pr),nr=1))
for(i in 1:length(pr)) plot(density(sf[rs,pr[i]]),main=paste("Prior",i))
#filled.contour(kde2d(sf[rs,pr[1]],sf[rs,pr[2]]),color.palette=heat.colors)
}

SSC.Like<-function(statfile,cols=NA,obs=NA,accept=.1) {
read.csv(statfile,na.strings="-1.#IND")->sf
if(is.na(cols[1])) {
print(colnames(sf))
readline("Which column/s (eg 4,23,27)? ")->cols
as.numeric(strsplit(cols,',')[[1]])->cols
}
if(is.na(obs[1])) {
readline("Observed values: ")->obs
as.numeric(strsplit(obs,',')[[1]])->obs
}
mean(sf[,cols],na.rm=T)->cx
euc=0
for(i in 1:length(cols)) 
	euc=euc+((sf[,cols[i]]-cx[i])/diff(range(sf[,cols[i]],na.rm=T)))^2
sqrt(euc)->euc
sf=sf[-which(euc>quantile(euc,.9,na.rm=T) | is.na(euc)),]
sqrt(length(cols))*accept->accept
euc=0
for(i in 1:length(cols)) 
	euc=euc+((sf[,cols[i]]-obs[i])/diff(range(sf[,cols[i]],na.rm=T)))^2
sqrt(euc)->euc
which(euc<accept & !is.na(euc))->ix
floor(sqrt(length(cols)))->flc
layout(matrix(1:(flc*ceiling(length(cols)/flc)),nc=flc))
par(mar=c(3,3,2,.5))
for(i in 1:length(cols)) {
	pretty(range(c(sf[,cols[i]],obs[i]),na.rm=T),n=30)->buckets
	hist(sf[,cols[i]],main=paste(i,": ",colnames(sf)[cols[i]],sep=""),breaks=buckets,xlab="",ylab="")
	hist(sf[ix,cols[i]],col="blue",breaks=buckets,add=T)
	points(obs[i],0,pch="*",col="red",cex=4)
}
return(length(ix)/nrow(sf))
}

approx.posterior<-function(x) {
#x should be a list of "accepted" values
fit.distrib<-function(p,t) {
	ex=switch(t,
		exp=dexp(x,p),
		unif=dunif(x,min(p),max(p)),
		norm=dnorm(x,p[1],p[2]),
		gam=dgamma(x,p[1],p[2]))
	-sum(log(ex))
}
edgeless.density<-function(x,refl=2) {
	#refl: # of bandwidths to extend data
	range(x)->z
	2*bw.nrd0(x)->b
	which(x<z[1]+b)->ix.l
	which(x>z[2]-b)->ix.h
	density(c(2*z[1]-x[ix.l],x,2*z[2]-x[ix.h]))->r
	which(r$x<z[1])->ix.l
	ix.l[2:length(ix.l)-1]->ix.l
	which(r$x>z[2])[-1]->ix.h
	r$x[-c(ix.l,ix.h)]->r$x
	r$y[-c(ix.l,ix.h)]->r$y
	r
}
edgeless.density(x)->dx
exp(log(diff(range(x,na.rm=T)))-8)->ddx
c(min(x)-ddx,dx$x,max(x)+ddx)->dx$x
c(0,dx$y,0)->dx$y
gof=c()
gof[[1]]=optimize(fit.distrib,t="exp",c(1/mean(x)*c(.1,10)))
gof[[2]]=optim(range(x),fit.distrib,t="unif")
gof[[3]]=optim(c(mean(x),sd(x)),fit.distrib,t="norm")
gof[[4]]=optim(c(1,1/mean(x)),fit.distrib,t="gam")
plot(dx,ylim=c(0,max(dx$y)*1.23),xlab="Parameter value",main="Approximate Fits")
polygon(c(min(x)-ddx,dx$x,max(x)+ddx),c(0,dx$y,0),col="plum")
lines(dx$x,dexp(dx$x,gof[[1]]$minimum),lwd=2,col="blue")
lines(dx$x,dunif(dx$x,gof[[2]]$par[1],gof[[2]]$par[2]),lwd=2,col="red")
lines(dx$x,dnorm(dx$x,gof[[3]]$par[1],gof[[3]]$par[2]),lwd=2,col="orange")
lines(dx$x,dgamma(dx$x,gof[[4]]$par[1],gof[[4]]$par[2]),lwd=2,col="black")
c(gof[[1]]$objective,gof[[2]]$value,gof[[3]]$value,gof[[4]]$value)->gofv
sort(gofv,index.return=T)$ix->ix
legend("topright",paste(signif(gofv,4),c("Exp(","Unif(","Norm(","Gamma("),
	c(signif(1/gof[[1]]$minimum,3),paste(signif(gof[[2]]$par,3),collapse=","),
	paste(signif(gof[[3]]$par,3),collapse=","),
	paste(signif(gof[[4]]$par,3),collapse=",")),")")[ix],text.col=c("blue","red","orange","black")[ix],
	cex=.8
	)
names(gof)<-c("Exp","Unif","Norm","Gamma")
gof[ix]
}


#Beaumont's method for estimating parameter values
makepd4 <- function(target,x,sumstat,tol,gwt,rejmethod=T,transf="none",bb=c(0,0)) {
# target is the set of target summary stats
# x is the parameter vector (long vector of numbers from the simulations) and is the dependent variable for the regression
# sumstat is an array of simulated summary stats (i.e. independent variables).
# NBB this function originally used lm() and assumed 4 summary stats, and I edited by hand for other numbers. 
# NBB I've now modified it using lsfit() (following Shola Ajayi) so that it will take an arbitrary number of summary stats.
# tol is the required proportion of points nearest the target values
# gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
# if rejmethod=T it doesn't bother with the regression, and just does rejection.


# If rejmethod=F it returns a list with the following components:-

# $x regression adjusted values
# $vals - unadjusted values in rejection region (i.e. normal rejection)
# $wt - the regression weight (i.e. the Epanechnikov weight)
# $ss - the sumstats corresponding to these points
# $predmean - estimate of the posterior mean
# $fv - the fitted value from the regression

if(sum(transf == c("none","log","logit")) == 0) stop("transf must be none, log, or logit")
if(transf=="logit") if(bb[1] >= bb[2]) stop("bounds wrong for logit")
if(missing(gwt))gwt <- rep(T,length(sumstat[,1]))
nss <- length(sumstat[1,])
# scale everything 
    scaled.sumstat <- sumstat  
    for(j in 1:nss) scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
    target.s <- target
    for(j in 1:nss) target.s[j] <- normalise(target[j],sumstat[,j][gwt])
    
# calc euclidean distance
    sum1 <- 0
    for(j in 1:nss) sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
   dst <- sqrt(sum1)
# includes the effect of gwt in the tolerance
    dst[!gwt] <- floor(max(dst[gwt])+10)
    

# wt1 defines the region we're interested in 
    abstol <- quantile(dst,tol)
    wt1 <- dst < abstol
    
    if(transf == "log"){
    	if(min(x) <= 0){
    		print("log transform: val out of bounds - correcting")
    		x.tmp <- ifelse(x <= 0,max(x),x)
    		x.tmp.min <- min(x.tmp)
    		x <- ifelse(x <= 0, x.tmp.min,x)
    	}
    	x <- log(x)
    }
    else if(transf == "logit"){
    	if(min(x) <= bb[1]){
    		x.tmp <- ifelse(x <= bb[1],max(x),x)
    		x.tmp.min <- min(x.tmp)
    		x <- ifelse(x <= bb[1], x.tmp.min,x)
    	}
    	if(max(x) >= bb[2]){
    		x.tmp <- ifelse(x >= bb[2],min(x),x)
    		x.tmp.max <- max(x.tmp)
    		x <- ifelse(x >= bb[2], x.tmp.max,x)
    	}
    	x <- (x-bb[1])/(bb[2]-bb[1])
    	x <- log(x/(1-x))
    }
    
    if(rejmethod) { l1 <- list(x=x[wt1],wt=wt1)
    } else {
        regwt <- 1-dst[wt1]^2/abstol^2        
        fit1 <- lsfit(scaled.sumstat[wt1,],x[wt1],wt=regwt)
        predmean <- fit1$coeff %*% c(1,target.s)
        l1 <- list(x=fit1$residuals+predmean,vals=x[wt1],wt=regwt,ss=sumstat[wt1,],predmean=predmean,fv = x[wt1]-fit1$residuals,transf=transf)
    }
    if(transf == "log"){
    	l1$x <- exp(l1$x)
    	l1$vals <- exp(l1$vals)
    } else if(transf == "logit"){
    	l1$x <- exp(l1$x)/(1+exp(l1$x))
    	l1$x <- l1$x*(bb[2]-bb[1])+bb[1]
    	l1$vals <- exp(l1$vals)/(1+exp(l1$vals))
    	l1$vals <- l1$vals*(bb[2]-bb[1])+bb[1]
    }
    
    l1
}


normalise <- function(x,y) {
if(var(y) == 0)return (x - mean(y))
(x-(mean(y)))/sqrt(var(y))
}


calmod <- function(target,x,sumstat,tol,gwt,rejmethod=T) {
# this function uses categorical regression to estimate the posterior probability of a particular model 
#      under the ABC framework P(Y=y | S = s)
#
# target is the set of target summary stats - i.e. s, what you measured from the data.
# x is the parameter vector, Y (long vector of numbers from the simulations) and is the dependent variable for the regression
#          This is a categorical variable, taking values 1, 2, .. n where there are n categories (models)
# sumstat is an array of simulated summary stats, S (i.e. independent variables).
# tol is the required proportion of points nearest the target values
# gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
# if rejmethod=T it doesn't bother with the regression, and just does rejection.


# If rejmethod=F it returns a list with the following components:-

# $x1 expected value on a logit scale, with standard errors of the estimate
# $x2 expected value on a natural scale - i.e. p(Y=y | S = s) This is what we would normally report.
# $vals - the Ys in the rejection region. The proportion of the different Ys (different models) is a Pritchard etal-style, rejection-based
#                                         estimate of the posterior probability of the different models. You might get some improvement by 
#                                          weighting the frequencies with $wt.
# $wt - the Epanechnikov weight. 
# $ss - the sumstats corresponding to the points in $vals. 



if(missing(gwt))gwt <- rep(T,length(sumstat[,1]))

nss <- length(sumstat[1,])


# scale everything 

    scaled.sumstat <- sumstat
    
    for(j in 1:nss){
    
    	scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
    }
    target.s <- target
    
    for(j in 1:nss){
    
    	target.s[j] <- normalise(target[j],sumstat[,j][gwt])
    }
    
# calc euclidean distance

    sum1 <- 0
    for(j in 1:nss){
    	sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
   }
   dst <- sqrt(sum1)
# includes the effect of gwt in the tolerance
    dst[!gwt] <- floor(max(dst[gwt])+10)
    

# wt1 defines the region we're interested in 
    abstol <- quantile(dst,tol)
    wt1 <- dst < abstol
    nit <- sum(wt1)
    
    if(rejmethod){
        l1 <- list(x=x[wt1],wt=0)
    }
    else{
        regwt <- 1-dst[wt1]^2/abstol^2
        
        catx <- as.numeric(names(table(x[wt1])))
        ncat <- length(catx)
        yy <- matrix(nrow=nit,ncol=ncat)
        for(j in 1:ncat){
        	yy[,j] <- as.numeric(x[wt1] == catx[j])
        }
        
        tr <- list()

        for(j in 1:nss){
        	tr[[j]] <- scaled.sumstat[wt1,j]
        }
        
        xvar.names <- paste("v",as.character(c(1:nss)),sep="")
        
        names(tr) <- xvar.names
        
        fmla <- as.formula(paste("yy ~ ", paste(xvar.names, collapse= "+")))
        
#        fit1 <- vglm(fmla,data=tr,multinomial) this is the version described in the 
#               manuscript,which did not use the Epanechnikov weights. 
        

        fit1 <- vglm(fmla,data=tr,multinomial,weights=regwt)
        
        target.df <- list()
        for(j in 1:nss){
        	target.df[[j]] <- target.s[j]
        }
        names(target.df) <- xvar.names
        
        prediction1 <- predict.vglm(fit1,target.df,se.fit=T)
        prediction2 <- predict.vglm(fit1,target.df,type="response")
        
        l1 <- list(x1=prediction1,x2=prediction2,vals=x[wt1],wt=regwt,ss=sumstat[wt1,])
        
    }
    l1
}

#Make pretty, stacked haplotype networks, v1.5
TempNet<-function(file,ftype="fasta",ages=NA,mut_size=1.5,nohap_size=.5,layernm=NA,
	invert=F,planes=F,confirm=T,color=NA,vcol=F,theta=0,phi=pi/6,hapid=F) { 

	in3d<-function(pts,theta,drawcmd,phi=pi/6,...) { #pts must be a matrix where columns are x,y,z coords; returns 2d coords
		rotm<-matrix(c(cos(theta),-sin(theta),0,0,-sin(phi),cos(phi)),nrow=3)
		if(drawcmd=="seg") {
			pts[,1:2]<-pts[,1:3]%*%rotm
			pts[,3:4]<-pts[,4:6]%*%rotm
			pts[,1:4]->pts
		} else as.matrix(pts)%*%rotm->pts
		if(drawcmd=="pl") plot(pts,...)
		if(drawcmd=="pt") points(pts,...)
		if(drawcmd=="poly") polygon(pts,...)
		if(drawcmd=="seg") segments(pts[,1],pts[,2],pts[,3],pts[,4],...)
		if(drawcmd=="line") lines(pts,...)
		if(drawcmd=="text") text(pts,...)
		pts
	}

	z0plane<-function(x) c((x[1]-x[2]*sin(theta)/sin(phi))/cos(theta),-x[2]/sin(phi))
	#ages must be specified as 1,2,3,...; 1 will be on the bottom
    library(pegas); library(ape)
    read.dna(file,format=ftype)->temp
	nseq=0
    if (is.na(ages[1])) {
		ag <- dimnames(temp)[[1]]
		if(is.null(ag)) { ag<-names(temp); nseq=length(temp) } else nseq=nrow(temp)
		ages <- as.integer(gsub("^.*[$]", "", ag)[attr(regexpr("^.*[$]",ag),"match.length")>0])
		if (length(ages)==0) ages <- rep(1,nseq)
	} else { 
		ages<-as.integer(factor(ages))
		length(ages)->nseq
	}
	if(sum(is.na(ages))) warning("Could not assign all samples to a layer.")
	
    if(nseq!=length(ages)) warning("Number of ages does not match number of sequences!")
    haplotype(temp)->s
    haploNet(s)->net
    if(is.na(layernm[1])) layernm=paste("Layer",1:length(unique(ages)))
	if(is.na(color[1])) color=rainbow(length(unique(ages)))
	if(invert) { ages=max(ages)-ages+1; layernm<-rev(layernm) }
        
   	matrix(0,nr=nrow(s),nc=length(unique(ages)))->lvl
    for(i in 1:nrow(s)) for(j in unique(ages)) lvl[i,j]=sum(ages[attr(s,"index")[[i]]]==j)
        
    table(c(net[,1:2]))->links
    elip=matrix(NA,nr=nrow(s),nc=3)
    elip[net[1,1],1:3]=0
    for(i in 1:nrow(net)) {
		while(sum(is.na(elip[net[i,1:2],1]))==2) #neither is yet part of the network
			net=rbind(net[-i,],net[i,])     
        net[i,1+is.na(elip[net[i,1:2],1])]->newl #[from,to] haplogroup
		which.max(apply(matrix(lvl[newl,],nrow=2),2,function(x) sum(sqrt(x/pi))))->bigcircs #find age group with max circle width
        d=sum(sqrt(lvl[newl,bigcircs]/pi))+net[i,3] #dist b/wn circle centers
		elip[newl[2],1:2]=elip[newl[1],1:2]+d*c(cos(elip[newl[1],3]),sin(elip[newl[1],3]))
        elip[newl[2],3]=elip[newl[1],3]+pi
        elip[newl,3]=elip[newl,3]+2*pi/links[newl]
    }
    #draw dat crazy graph
    par()$mar->oldmar
    clk=1; lbl.coords=NA
	windows(title="TempNet (press ESC when done)"); par(mar=c(0,0,0,0))
    while(clk!=0) {
		range(c(elip[,2]-sqrt(lvl/pi),elip[,2]+sqrt(lvl/pi)))->rng
		diff(rng)/2->h #z-distance between layers, needs to be *sin(phi) eventually...
		range(elip[,1])+c(-max(sqrt(lvl[which.min(elip[,1]),]/pi)),max(sqrt(lvl[which.max(elip[,1]),]/pi)))->xrng
        in3d(expand.grid(xrng*(1+planes/10),c(rng[1],rng[2]),h*(1:ncol(lvl)-1)),theta,"pl",phi,type="n")
		if(is.na(lbl.coords[1])) lbl.coords=c(xrng[1],rng[1])
		else {
			if(lbl.coords[1]<par('usr')[1]) lbl.coords[1]=par('usr')[1]
			if(lbl.coords[2]<par('usr')[3]) lbl.coords[2]=par('usr')[3]
		}
        seq(0,2*pi,len=31)->seg #can't draw real elipses, so draw a regular 30-gon instead
        for(i in 1:ncol(lvl)) {
			if(planes) in3d(cbind(c(xrng,rev(xrng)),rep(rng,e=2),h*(i-1)),theta,"poly",phi,col=rgb(.2,.2,.2,a=.4))
			for(j in rev(sort(elip[,2],ind=T)$ix)) { #elipses
				sqrt(lvl[j,i]/pi)->r; if(r==0) r=sqrt(nohap_size)/pi
				in3d(cbind(r*cos(seg)+elip[j,1],r*sin(seg)+elip[j,2],h*(i-1)),theta,"poly",phi,col=ifelse(lvl[j,i]==0,"white",color[i]))
				if(i<ncol(lvl)) if(lvl[j,i+1]>0 && lvl[j,i]>0) #exists at adjacent time points, vertical lines
					in3d(cbind(elip[j,1]+c(-r,r),rep(elip[j,2],2),h*(i-1),elip[j,1]+c(-1,1)*sqrt(lvl[j,i+1]/pi),rep(elip[j,2],2),h*i), theta,"seg",phi,col=ifelse(vcol,color[i+1],"black"))
			} 
            for(j in rev(sort(elip[,2],ind=T)$ix)) if(lvl[j,i]>0) in3d(cbind(elip[j,1],elip[j,2],h*(i-1)),theta,"text",phi,labels=ifelse(hapid,j,lvl[j,i])) #n in haplogroup
            for(j in 1:nrow(net)) { #lines; solid if both haplogroups exist at that time, dashed if not
				sqrt(lvl[net[j,1:2],i]/pi)->r; r[r==0]=sqrt(nohap_size)/pi
				th=atan2(diff(elip[net[j,1:2],2]),diff(elip[net[j,1:2],1]))
                x=elip[net[j,1:2],1]+c(1,-1)*r*cos(th)
                y=elip[net[j,1:2],2]+c(1,-1)*r*sin(th)
				in3d(cbind(x,y,h*(i-1)),theta,"line",phi,lwd=1+(sum(lvl[net[j,1:2],i]==0)==0),lty=c("dotted","solid")[1+(sum(lvl[net[j,1:2],i]==0)==0)])
                if(net[j,3]>1) in3d(cbind(seq(from=x[1],to=x[2],len=net[j,3]+1)[-c(1,net[j,3]+1)], seq(from=y[1],to=y[2],len=net[j,3]+1)[-c(1,net[j,3]+1)],h*(i-1)),theta,"pt",phi,pch=20,cex=mut_size)
            }
			in3d(matrix(c(lbl.coords,h*(i-1)),nrow=1),theta,"text",phi,label=layernm[i],col=color[i],pos=4,cex=1.5) #layer labels
		}
		if(confirm==T) {
			spot=z0plane(unlist(locator(1)))		
			apply(elip[,1:2]-rep(spot,e=nrow(elip)),1,function(x) sum(x^2))->dst
			if(min(dst)<1) {
				clk=which.min(dst)
				elip[clk,1:2]=z0plane(unlist(locator(1)))
			} else { 
				if(sum(spot-lbl.coords)<1) lbl.coords=z0plane(unlist(locator(1)))
				else cat("\a")
			}
		} else clk=0
	}
        
	par(mar=oldmar)
	rep(NA,nrow(temp))->x
	names(x)=rownames(temp)
	attr(s,"index")->y
	for(i in 1:length(y)) x[y[[i]]]=i
	list(x,lvl,net)
}

na.per.row<-function(x) apply(x,1,function(y) sum(is.na(y)))

viewlog<-function(path,logNe=F,demenames=NA,maxt=NA) {
	#make a spindle graph of Nes and Ks through time
	read.table(path,head=T)->z
	ncol(z)/2-1->npop
	if(logNe) for(i in 1:npop) z[,i*2+1]=log(z[,i*2+1])
	which.max(apply(z[,1:npop*2+1],1,sum))->tmax
	sum(z[tmax,1:npop*2+1])->NeScl  #largest census size
	xc=as.numeric(c(0,cumsum(as.numeric(z[tmax,1:npop*2+1]))[-npop])+z[tmax,1:npop*2+1]/2) #center lines
	
	plot(NA,xlim=c(0,NeScl),ylim=c(0,ifelse(is.na(maxt),max(z$Time),maxt)),xlab="Census (Ne)",ylab="Time (yr bp)",axes=F)
	axis(2)
	if(logNe) for(i in 1:npop) {
		floor(log10(exp(max(z[,i*2+1]))))->j
		axis(1,at=xc[i]+log(10^(0:j))/2,labels=c("1","10","100","1k","10k","100k","1M","10M")[0:j+1])
	} else for(i in 1:npop) {
		pretty(c(0,max(z[,i*2+1])))->tick
		axis(1,at=xc[i]+tick/2,labels=tick)
	}
	for(i in 1:npop) {
		data.frame(x=c(xc[i]-z[,i*2+1]/2,xc[i]+rev(z[,i*2+1]/2)),y=c(z[,2],rev(z[,2])))->lin
		polygon(lin[which(lin[,1]!=Inf & lin[,1]!=-Inf),],col=rainbow(npop)[i])
		ps=z[,i*2+2]/z[,i*2+1] #percent of deme sampled
		data.frame(x=c(xc[i]-z[,i*2+2]/(2*max(ps)),xc[i]+rev(z[,i*2+2]/(2*max(ps)))),y=c(z[,2],rev(z[,2])))->lin
		polygon(lin[which(c(ps,rev(ps))>0),],col="black")
		mtext(ifelse(is.na(demenames[1]),paste("Deme",i-1),demenames[i]),
		     at=xc[i],col=rainbow(npop)[i],cex=1.7,line=.5)
	}

}
