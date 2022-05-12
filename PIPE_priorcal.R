#PIPE Prior Calibration
library(pipe.design)
set.seed(200)

#values grid
low=c(0.025,0.05,0.075,0.10)
inc=c(0.025,0.05,0.075,0.10)
#function for prior median matrix
prior.m=function(low,inc){
  t1=low
  t2=low+inc
  t3=low+2*inc
  t4=low+3*inc
  t5=low+4*inc
  prior.med=matrix(nrow=3,ncol=3,c(t1,t2,t3,t2,t3,t4,t3,t4,t5))
  return(prior.med)
}
ss=c(1/72,1/36,1/18,1/9)

#scenarios to calibrate over
true1=matrix(nrow=3,ncol=3,c(0.05,0.1,0.15,0.1,0.15,0.2,0.15,0.2,0.3))
true2=matrix(nrow=3,ncol=3,c(0.05,0.1,0.3,0.1,0.2,0.45,0.2,0.3,0.55))
true3=matrix(nrow=3,ncol=3,c(0.15,0.3,0.45,0.3,0.45,0.55,0.45,0.55,0.65))
true4=matrix(nrow=3,ncol=3,c(0.3,0.45,0.5,0.45,0.5,0.55,0.5,0.55,0.6))


sc1=array(NA,dim=c(length(low),length(inc),length(ss)))
for(x1 in 1:length(low)){
  for(x2 in 1:length(inc)){
    for(x3 in 1:length(ss)){
      design = pipe.design(N=36,S=2000,c=3,theta=0.30,pi=true1, 
                           prior.med=prior.m(low[x1],inc[x2]), prior.ss=matrix(nrow=3,ncol=3,ss[x3]),
                           strategy="ss-random",admis="closest",constraint="neighbouring-nodiag",
                           uppertox.constraint=NULL,epsilon=1,mode="sim")
      dose = pickone(design)
      mtd = which(true1==0.3,arr.ind=T)
      same = matrix(NA,nrow=2000,ncol=dim(mtd)[1])
      for(q in 1:2000){
        for(z in 1:dim(mtd)[1]){
          same[q,z] = 1*(mtd[z,1]==dose[q,1] & mtd[z,2]==dose[q,2])
        }
      }
      success=apply(same,1,sum)
      sc1[x1,x2,x3] = 100*mean(success)
    }
  }
}
round(sc1,1)

sc2=array(NA,dim=c(length(low),length(inc),length(ss)))
for(x1 in 1:length(low)){
  for(x2 in 1:length(inc)){
    for(x3 in 1:length(ss)){
      design = pipe.design(N=36,S=2000,c=3,theta=0.3,pi=true2, 
                           prior.med=prior.m(low[x1],inc[x2]), prior.ss=matrix(nrow=3,ncol=3,ss[x3]),
                           strategy="ss-random",admis="closest",constraint="neighbouring-nodiag",
                           uppertox.constraint=NULL,epsilon=1,mode="sim")
      dose = pickone(design)
      mtd = which(true2==0.3,arr.ind=T)
      same = matrix(NA,nrow=2000,ncol=dim(mtd)[1])
      for(q in 1:2000){
        for(z in 1:dim(mtd)[1]){
          same[q,z] = (mtd[z,1]==dose[q,1] & mtd[z,2]==dose[q,2])
        }
      }
      success=apply(same,1,sum)
      sc2[x1,x2,x3] = 100*mean(success)
    }
  }
}
round(sc2,1)

sc3=array(NA,dim=c(length(low),length(inc),length(ss)))
for(x1 in 1:length(low)){
  for(x2 in 1:length(inc)){
    for(x3 in 1:length(ss)){
      design = pipe.design(N=36,S=2000,c=3,theta=0.3,pi=true3, 
                           prior.med=prior.m(low[x1],inc[x2]), prior.ss=matrix(nrow=3,ncol=3,ss[x3]),
                           strategy="ss-random",admis="closest",constraint="neighbouring-nodiag",
                           uppertox.constraint=NULL,epsilon=1,mode="sim")
      dose = pickone(design)
      mtd = which(true3==0.3,arr.ind=T)
      same = matrix(NA,nrow=2000,ncol=dim(mtd)[1])
      for(q in 1:2000){
        for(z in 1:dim(mtd)[1]){
          same[q,z] = (mtd[z,1]==dose[q,1] & mtd[z,2]==dose[q,2])
        }
      }
      success=apply(same,1,sum)
      sc3[x1,x2,x3] = 100*mean(success)
    }
  }
}
round(sc3,1)

sc4=array(NA,dim=c(length(low),length(inc),length(ss)))
for(x1 in 1:length(low)){
  for(x2 in 1:length(inc)){
    for(x3 in 1:length(ss)){
      design = pipe.design(N=36,S=2000,c=3,theta=0.3,pi=true4, 
                           prior.med=prior.m(low[x1],inc[x2]), prior.ss=matrix(nrow=3,ncol=3,ss[x3]),
                           strategy="ss-random",admis="closest",constraint="neighbouring-nodiag",
                           uppertox.constraint=NULL,epsilon=1,mode="sim")
      dose = pickone(design)
      mtd = which(true4==0.3,arr.ind=T)
      same = matrix(NA,nrow=2000,ncol=dim(mtd)[1])
      for(q in 1:2000){
        for(z in 1:dim(mtd)[1]){
          same[q,z] = (mtd[z,1]==dose[q,1] & mtd[z,2]==dose[q,2])
        }
      }
      success=apply(same,1,sum)
      sc4[x1,x2,x3] = 100*mean(success)
    }
  }
}
round(sc4,1)

sc1
sc2
sc3
sc4

########################################################################

## Average PCS across scenarios 

#geometric mean
gm=mm=array(NA,dim=c(length(low),length(inc),length(ss)))
for(i in 1:length(low)){
  for(j in 1:length(inc)){
    for(k in 1:length(ss)){
      gm[i,j,k]=(sc1[i,j,k]*sc2[i,j,k]*sc3[i,j,k]*sc4[i,j,k])^(0.25)
      mm[i,j,k]=min(c(sc1[i,j,k],sc2[i,j,k],sc3[i,j,k],sc4[i,j,k]))
    }
  }
}
round(gm,1)
round(mm,1)
min(gm)
which(gm==max(gm),arr.ind=T)
which(gm>38,arr.ind=T)
which(mm==max(mm),arr.ind=T)


#graphs
library(ggplot2)
x <- c("0.025","0.05","0.075","0.10")
y <- c("0.025","0.05","0.075","0.10")
data <- expand.grid(X=x, Y=y)
data$GM1 <- as.vector(gm[,,1])
data$GM2 <- as.vector(gm[,,2])
data$GM3 <- as.vector(gm[,,3])
data$GM4 <- as.vector(gm[,,4])

pdf(file="fig_PIPE_prior1.pdf",width=4,height=3)
ggplot(data,aes(x=data$X,y=data$Y,fill=(data$GM1))) + geom_tile() + scale_fill_gradient(low="white", high="darkblue",limits=c(30,40)) +
  labs(x=expression(rho),y=expression(delta),fill="Mean \nPCS (%)") + theme(legend.key.height=unit(1, "cm")) +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))
dev.off()

pdf(file="fig_PIPE_prior2.pdf",width=4,height=3)
ggplot(data,aes(x=data$X,y=data$Y,fill=(data$GM2))) + geom_tile() + scale_fill_gradient(low="white", high="darkblue",limits=c(30,40)) +
  labs(x=expression(rho),y=expression(delta),fill="Mean \nPCS (%)") + theme(legend.key.height=unit(1, "cm")) +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))
dev.off()

pdf(file="fig_PIPE_prior3.pdf",width=4,height=3)
ggplot(data,aes(x=data$X,y=data$Y,fill=(data$GM3))) + geom_tile() + scale_fill_gradient(low="white", high="darkblue",limits=c(30,40)) +
  labs(x=expression(rho),y=expression(delta),fill="Mean \nPCS (%)") + theme(legend.key.height=unit(1, "cm")) +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))
dev.off()

pdf(file="fig_PIPE_prior4.pdf",width=4,height=3)
ggplot(data,aes(x=data$X,y=data$Y,fill=(data$GM4))) + geom_tile() + scale_fill_gradient(low="white", high="darkblue",limits=c(30,40)) +
  labs(x=expression(rho),y=expression(delta),fill="Mean \nPCS (%)") + theme(legend.key.height=unit(1, "cm")) +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))
dev.off()

#########################################################################





## calibrating the overdosing rule
set.seed(800)

#take best values
low=0.05
inc=0.025
ss=1/18


epsilon=seq(from=1,to=0.4,by=-0.05)


s8=matrix(c(0.05,0.10,0.30,0.10,0.20,0.45,0.20,0.30,0.55),nrow=3,ncol=3)
s10=matrix(c(0.15,0.30,0.45,0.30,0.45,0.55,0.45,0.55,0.65),nrow=3,ncol=3)
s13=matrix(c(0.30,0.45,0.50,0.45,0.50,0.55,0.50,0.55,0.60),nrow=3,ncol=3)
s14=matrix(c(0.45,0.50,0.55,0.50,0.55,0.60,0.55,0.60,0.65),nrow=3,ncol=3)


tox_s8<-c()
overdoses8<-matrix(NA,nrow=length(epsilon),ncol=2000)
for(y1 in 1:length(epsilon)){
  simulation = pipe.design(N=36,S=2000,c=3,theta=0.3,pi=s8, 
                           prior.med=prior.m(low,inc), prior.ss=matrix(nrow=3,ncol=3,ss),
                           strategy="ss-random",admis="closest",constraint="neighbouring-nodiag",
                           uppertox.constraint=NULL,epsilon=epsilon[y1],mode="sim")
  overdoses8[y1,]<-sapply(c(1:2000), function(x) sum(simulation$n.sim[[x]]*as.numeric(simulation$pi>0.3)))
  dose = pickone(simulation)
  mtd = which(s8==0.3,arr.ind=T)
  same = matrix(NA,nrow=2000,ncol=dim(mtd)[1])
  for(q in 1:2000){
    for(z in 1:dim(mtd)[1]){
      same[q,z] = (mtd[z,1]==dose[q,1] & mtd[z,2]==dose[q,2])
    }
  }
  success=apply(same,1,sum)
  tox_s8[y1]=100*mean(success)
}
tox_s10<-c()
overdoses10<-matrix(NA,nrow=length(epsilon),ncol=2000)
for(y1 in 1:length(epsilon)){
  simulation = pipe.design(N=36,S=2000,c=3,theta=0.3,pi=s10, 
                           prior.med=prior.m(low,inc), prior.ss=matrix(nrow=3,ncol=3,ss),
                           strategy="ss-random",admis="closest",constraint="neighbouring-nodiag",
                           uppertox.constraint=NULL,epsilon=epsilon[y1],mode="sim")
  overdoses10[y1,]<-sapply(c(1:2000), function(x) sum(simulation$n.sim[[x]]*as.numeric(simulation$pi>0.3)))
  dose = pickone(simulation)
  mtd = which(s10==0.3,arr.ind=T)
  same = matrix(NA,nrow=2000,ncol=dim(mtd)[1])
  for(q in 1:2000){
    for(z in 1:dim(mtd)[1]){
      same[q,z] = (mtd[z,1]==dose[q,1] & mtd[z,2]==dose[q,2])
    }
  }
  success=apply(same,1,sum)
  tox_s10[y1]=100*mean(success)
}

tox_s13<-c()
overdoses13<-matrix(NA,nrow=length(epsilon),ncol=2000)
for(y1 in 1:length(epsilon)){
  simulation = pipe.design(N=36,S=2000,c=3,theta=0.3,pi=s13, 
                           prior.med=prior.m(low,inc), prior.ss=matrix(nrow=3,ncol=3,ss),
                           strategy="ss-random",admis="closest",constraint="neighbouring-nodiag",
                           uppertox.constraint=NULL,epsilon=epsilon[y1],mode="sim")
  overdoses13[y1,]<-sapply(c(1:2000), function(x) sum(simulation$n.sim[[x]]*as.numeric(simulation$pi>0.3)))
  dose = pickone(simulation)
  mtd = which(s13==0.3,arr.ind=T)
  same = matrix(NA,nrow=2000,ncol=dim(mtd)[1])
  for(q in 1:2000){
    for(z in 1:dim(mtd)[1]){
      same[q,z] = (mtd[z,1]==dose[q,1] & mtd[z,2]==dose[q,2])
    }
  }
  success=apply(same,1,sum)
  tox_s13[y1]=100*mean(success)
}

tox_s14<-c()
overdoses14<-matrix(NA,nrow=length(epsilon),ncol=2000)
for(y1 in 1:length(epsilon)){
  simulation = pipe.design(N=36,S=2000,c=3,theta=0.3,pi=s14, 
                           prior.med=prior.m(low,inc), prior.ss=matrix(nrow=3,ncol=3,ss),
                           strategy="ss-random",admis="closest",constraint="neighbouring-nodiag",
                           uppertox.constraint=NULL,epsilon=epsilon[y1],mode="sim")
  overdoses14[y1,]<-sapply(c(1:2000), function(x) sum(simulation$n.sim[[x]]*as.numeric(simulation$pi>0.3)))
  pstop = sum(simulation$n.rpII==0)/2000
  tox_s14[y1] = 100*pstop
}
#plots
pipe_plot8a<-data.frame(epsilon,rowMeans(overdoses8),"Sc 8")
names(pipe_plot8a)<-c("epsilon","value","Scenario")
pipe_plot10a<-data.frame(epsilon,rowMeans(overdoses10),"Sc 10")
names(pipe_plot10a)<-c("epsilon","value","Scenario")
pipe_plot13a<-data.frame(epsilon,rowMeans(overdoses13),"Sc 13")
names(pipe_plot13a)<-c("epsilon","value","Scenario")
pipe_plot14a<-data.frame(epsilon,rowMeans(overdoses14),"Sc 14")
names(pipe_plot14a)<-c("epsilon","value","Scenario")
pipe_plotGa<-data.frame(epsilon,(rowMeans(overdoses8)*rowMeans(overdoses10)*rowMeans(overdoses13)*rowMeans(overdoses14))^0.258,"Geo mean")
names(pipe_plotGa)<-c("epsilon","value","Scenario")

pipe_plota<-rbind(pipe_plot8a,pipe_plot10a,pipe_plot13a,pipe_plot14a,pipe_plotGa)
names(pipe_plota)<-c("epsilon","value","Scenario")

ggplot(pipe_plota) + geom_line(aes(x=epsilon, y=value, colour=Scenario)) +
  scale_colour_manual(values=c("red","green","blue","purple","black")) +
  ylim(0,36) +
  labs(x="Epsilon",y="Patients Treated on Over-toxic Dose (mean)",color=NULL) 

plot(epsilon,rowMeans(overdoses8),type="l",ylim=c(0,40))
lines(epsilon,rowMeans(overdoses10),col="red")
lines(epsilon,rowMeans(overdoses13),col="blue")
lines(epsilon,rowMeans(overdoses14),col="green")

plot(epsilon,tox_s8,type="l",ylim=c(0,100))
lines(epsilon,tox_s10,col="red")
lines(epsilon,tox_s13,col="blue")
lines(epsilon,tox_s14,col="green")



pipe_plot8<-data.frame(epsilon,tox_s8,"Sc 8")
names(pipe_plot8)<-c("epsilon","value","Scenario")
pipe_plot10<-data.frame(epsilon,tox_s10,"Sc 10")
names(pipe_plot10)<-c("epsilon","value","Scenario")
pipe_plot13<-data.frame(epsilon,tox_s13,"Sc 13")
names(pipe_plot13)<-c("epsilon","value","Scenario")
pipe_plot14<-data.frame(epsilon,tox_s14,"Sc 14")
names(pipe_plot14)<-c("epsilon","value","Scenario")
pipe_plotG<-data.frame(epsilon,(tox_s8*tox_s10*tox_s13*tox_s14)^0.258,"Geo mean")
names(pipe_plotG)<-c("epsilon","value","Scenario")
pipe_plot<-rbind(pipe_plot8,pipe_plot10,pipe_plot13,pipe_plot14,pipe_plotG)
names(pipe_plot)<-c("epsilon","value","Scenario")

ggplot(pipe_plot) + geom_line(aes(x=epsilon, y=value, colour=Scenario)) +
  scale_colour_manual(values=c("red","green","blue","purple","black")) +
  ylim(0,100) +
  labs(x="Epsilon",y="Correct Outcome (%)",color=NULL) 

#tox_sc8=vector()
overdoses8<-matrix(NA,nrow=length(epsilon),ncol=2000)
for(y1 in 1:length(epsilon)){
  simulation = pipe.design(N=36,S=2000,c=3,theta=0.3,pi=s8, 
                           prior.med=prior.m(low,inc), prior.ss=matrix(nrow=3,ncol=3,ss),
                           strategy="ss-random",admis="closest",constraint="neighbouring-nodiag",
                           uppertox.constraint=NULL,epsilon=epsilon[y1],mode="sim")
  overdoses8[y1,]<-sapply(c(1:2000), function(x) sum(simulation$n.sim[[x]]*as.numeric(simulation$pi>0.3)))
  #  pstop = sum(simulation$n.rpII==0)/2000
  # tox_sc1[y1] = 100*pstop
}
tox_sc1

tox_sc2a=vector()
for(y1 in 1:length(epsilon)){
  design = pipe.design(N=36,S=2000,c=3,theta=0.3,pi=trueSAFE, 
                       prior.med=prior.m(low,inc), prior.ss=matrix(nrow=3,ncol=3,ss),
                       strategy="ss-random",admis="closest",constraint="neighbouring-nodiag",
                       uppertox.constraint=NULL,epsilon=epsilon[y1],mode="sim")
  dose = pickone(design)
  mtd = which(trueSAFE==0.3,arr.ind=T)
  same = matrix(NA,nrow=2000,ncol=dim(mtd)[1])
  for(q in 1:2000){
    for(z in 1:dim(mtd)[1]){
      same[q,z] = (mtd[z,1]==dose[q,1] & mtd[z,2]==dose[q,2])
    }
  }
  success=apply(same,1,sum)
  tox_sc2a[y1]=100*mean(success)
}
tox_sc2a

pdf(file="fig_PIPE_safe.pdf",width=4.6,height=3.4)
ggplot() + geom_line(aes(x=epsilon,y=tox_sc1,col="red"),size=0.8) + geom_line(aes(x=epsilon,y=tox_sc2a,col="green"),size=0.8) + ylim(0,100) +
  labs(x=expression(epsilon),y="Correct Outcome (%)",color=NULL) +  scale_color_manual(values=c("green"="green","red"="red"), labels=c("Sc 13","Sc 14")) +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))
dev.off()
