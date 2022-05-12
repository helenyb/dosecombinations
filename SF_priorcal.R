#setwd("//lancs/homes/45/georgem1/My Documents/DISS/Surf")
library(rjags)

mono_prior=function(clinA,clinB,effSS){
  thetahat=tauhat=NULL
  ndoseA=length(clinA)
  ndoseB=length(clinB)
  thetahat[1]=tauhat[1]=1-clinA[1]-clinB[1]+clinA[1]*clinB[1]
  for(i in 2:(ndoseA)){
    thetahat[i]=(1-clinA[i])/(1-clinA[i-1])
  }
  for(j in 2:(ndoseB)){
    tauhat[j]=(1-clinB[j])/(1-clinB[j-1])
  }
  tauhat=tauhat[-1]
  
  c=rep(effSS,ndoseA)
  t=rep(effSS,ndoseB-1)
  a=c*thetahat
  b=c*(1-thetahat)
  e=t*tauhat
  f=t*(1-tauhat)
  
  return(list(thetahat=thetahat,tauhat=tauhat,a=a,b=b,e=e,f=f,
              ndoseA=ndoseA,ndoseB=ndoseB))
}

op_prior=function(mean.ratio,effSS,space){
  ndoseA=space[1];ndoseB=space[2]
  thetahat=rep(mean.ratio,space[1])
  tauhat=rep(mean.ratio,(space[2]-1))
  a=effSS*thetahat;b=effSS*(1-thetahat)
  e=effSS*tauhat;f=effSS*(1-tauhat)
  
  return(list(thetahat=thetahat,tauhat=tauhat,a=a,b=b,e=e,f=f))
}

MTD=function(theta,space,tt){
  tox=matrix(NA,nrow=space[1],ncol=space[2])
  for(i in 1:space[1]){
    for(j in 1:space[2]){
      if(j==1){
        tox[i,1] = 1 - prod(theta[1:i])
      }else{
        tox[i,j] = 1 - prod(theta[1:i]) * prod(theta[(space[1]+1):(space[1]+j-1)]) 
      }
    }
  }
  return(tox)
}

eli=function(ratio,space,tt){
  eliminate=matrix(NA,nrow=space[1],ncol=space[2])
  for(i in 1:space[1]){
    for(j in 1:space[2]){
      if(j==1){
        if(i==1){
          eliminate[i,1] = mean((1 - ratio[,1:i]) > tt)
        }else{
          eliminate[i,1] = mean((1 - apply(ratio[,1:i],1,prod)) > tt)
        }
      }else{
        if(j==2 & i==1){
          eliminate[i,j] = mean((1 - ratio[,1:i] * ratio[,(space[1]+1):(space[1]+j-1)]) > tt)
        }else if(j==2 & i>=2){
          eliminate[i,j] = mean((1 - apply(ratio[,1:i],1,prod) * ratio[,(space[1]+1):(space[1]+j-1)]) > tt)
        }else if(j>=3 & i==1){
          eliminate[i,j] = mean((1 - ratio[,1:i] * apply(ratio[,(space[1]+1):(space[1]+j-1)],1,prod)) > tt)
        }else{
          eliminate[i,j] = mean((1 - apply(ratio[,1:i],1,prod) * apply(ratio[,(space[1]+1):(space[1]+j-1)],1,prod)) > tt)
        }
      }
    }
  }
  return(eliminate)
}

cat("model{
    
    ### Define likelihood
    for (i in 1:m[1]){
    for(j in 1:m[2]){
    x[i,j] ~ dbin(p[i,j],n[i,j])
    p[i,j] = (1-step(j-1.5)) * (1 - prod(par[1:i])) +
    step(j-1.5) * (1 - prod(par[1:i]) * prod(par[(m[1]+step(j-1.5)):(m[1]+j-1)]) )
    }
    }
    ### Define the priors
    for (k in 1:(sum(m)-1)){
    par[k] ~ dbeta(hyper[k],hyper[k+sum(m)-1])T(0,0.999999)
    }
    
    }", file="model_sf.txt")

###################################################################################

SFdesign=function(clinA=NULL, clinB=NULL, space, prior.type, mean.ratio, effSS,
                  n.patients, n.dlts, no.diag=TRUE, safety.thresh, tt, currentdose){
  
  ### Prior knowledge and data
  if(prior.type=="operational"){
    prior=op_prior(mean.ratio,effSS,space)
  }
  if(prior.type=="monotherapy"){
    prior=mono_prior(clinA,clinB,effSS)
  }
  a=prior$a;b=prior$b
  e=prior$e;f=prior$f
  thetahat=prior$thetahat
  tauhat=prior$tauhat
  ndoseA=space[1];ndoseB=space[2]
  
  hyper=c(a,e,b,f)
  inits=list(par=c(thetahat,tauhat))
  params=rep(0,(ndoseA+ndoseB-1))
  for(i in 1:(ndoseA+ndoseB-1)){
    params[i]=paste("par[", i, "]", sep="")
  }
  n=matrix(nrow=ndoseA,ncol=ndoseB,byrow=F,n.patients)
  x=matrix(nrow=ndoseA,ncol=ndoseB,byrow=F,n.dlts)
  data=list("n" = n, "x" = x, "hyper" = hyper, "m" = c(ndoseA,ndoseB),
            "p" = matrix(nrow=ndoseA,ncol=ndoseB))
  
  ### Run JAGS model
  jags.m=jags.model(file="model_sf.txt",data=data,inits=inits,n.chains=1,n.adapt=500,quiet=TRUE)
  theta=rep(0,(ndoseA+ndoseB-1))
  ratio=matrix(NA,nrow=2000,ncol=(ndoseA+ndoseB-1))
  js=jags.samples(jags.m,params,n.iter=2000,type="trace",progress.bar="none")
  for(i in 1:(ndoseA+ndoseB-1)){
    t=paste("par[", i, "]", sep="")
    ratio[,i]=as.vector(js[[t]])
    theta[i]=mean(as.vector(js[[t]]))
  }
  
  ### Admissible doses
  tox=MTD(theta,space,tt)
  if(no.diag==TRUE){
    potentialdose=rbind(currentdose,currentdose+c(1,0),currentdose+c(0,1),currentdose+c(1,-1),
                        currentdose+c(-1,1),currentdose-c(1,0),currentdose-c(0,1),currentdose-c(1,1))
    potentialdose=potentialdose[which(potentialdose[,1]!=0),]
    potentialdose=potentialdose[which(potentialdose[,2]!=0),]
    potentialdose=potentialdose[which(potentialdose[,1]<=ndoseA),]
    potentialdose=potentialdose[which(potentialdose[,2]<=ndoseB),]
    potentialdose=matrix(potentialdose,nrow=dim(potentialdose)[1])
    
    ### Overdose control
    elim=eli(ratio=ratio,space=space,tt=tt)
    rej=rep(0,dim(potentialdose)[1])
    for(z in 1:dim(potentialdose)[1]){
      rej[z]=1*(elim[potentialdose[z,1],potentialdose[z,2]] >= safety.thresh)
    }
    potentialdose=matrix(potentialdose[which(rej==0),],ncol=2)
    
    if(dim(potentialdose)[1]==0){
      prob.overdose=1
      rule=c(0,0)
    }else{
      prob.overdose=0
      poss=rep(NA,dim(potentialdose)[1])
      for(i in 1:dim(potentialdose)[1]){
        a1=potentialdose[i,1];a2=potentialdose[i,2]
        poss[i]=tox[a1,a2]
      }
      b=which(abs(poss-tt)==min(abs(poss-tt)))
      rule=potentialdose[b,]
    }
  }
  
  ### MTD
  if(safety.thresh!=FALSE & prob.overdose==1){
    mtd=c(0,0)
  }else{
    mtd=which(abs(tox-tt)==min(abs(tox-tt)),arr.ind=TRUE)
  }
  
  return(list(prob.overdose=prob.overdose,n=n,x=x,tox=tox,rule=rule,mtd=mtd))
}

# SIMULATIONS -------------------------------------------------------------

sim=function(m, ss, true, start.dose, clinA=NULL, clinB=NULL, space, prior.type,
             mean.ratio, effSS, no.diag, safety.thresh, tt, n.sims){
  
  if(safety.thresh>1 | safety.thresh<=0) {
    stop("safety.thresh must be a number in the interval (0,1]")
  }
  if(ss/m != floor(ss/m)){
    stop("sample size ss must be divisible by cohort size m")
  }
  if(no.diag!=TRUE & no.diag!=FALSE){
    stop("no.diag must be set to TRUE or FALSE")
  }
  
  ndoseA=space[1];ndoseB=space[2]
  
  trueMAT=matrix(true,nrow=ndoseA,ncol=ndoseB)
  true.mtd=which(trueMAT==tt,arr.ind=T)
  true.mtds=which(trueMAT>0.16 & trueMAT<0.33,arr.ind=T)
  true.od=which(trueMAT>tt,arr.ind=T)
  
  n_track=list()
  y_track=list()
  mtd_track=list()
  
  PCS=vector()
  PAS=vector()
  PCE=vector()
  PAE=vector()
  OVER=vector()
  PEOVER=vector()
  ntox=vector()
  npat=vector()
  ZERO=vector()
  
  for(k in 1:n.sims){
    n=x=rep(0,ndoseA*ndoseB)
    rule=start.dose
    new=start.dose[1]+ndoseA*(start.dose[2]-1)
    abandon=0
    i=0
    while(abandon==0 & i<ss){
      i=i+m
      currentdose=rule
      n[new]=n[new]+m
      dlt=rbinom(1,size=m,prob=c(true[new],1-true[new]))
      x[new]=x[new]+dlt
      d=SFdesign(clinA=clinA, clinB=clinB, space=space, prior.type=prior.type, mean.ratio=mean.ratio, effSS=effSS,
                 n.patients=n, n.dlts=x, no.diag=no.diag, safety.thresh=safety.thresh, tt=tt, currentdose=currentdose)
      rule=d$rule
      new=rule[1]+ndoseA*(rule[2]-1)
      
      if(d$prob.overdose==1){
        abandon=1
      }
    }
    if(n.sims==1){
      return(d)
    }
    else{
      npat[k]=sum(n)
      ntox[k]=sum(x)
      n_track[[k]]=matrix(n,nrow=ndoseA,ncol=ndoseB)
      y_track[[k]]=matrix(x,nrow=ndoseA,ncol=ndoseB)
      mtd_track[[k]]=d$mtd
      
      if(dim(true.mtd)[1]!=0){
        correct=pce=vector()
        for(g in 1:dim(true.mtd)[1]){
          correct[g]=100*(d$mtd[1]==true.mtd[g,1] & d$mtd[2]==true.mtd[g,2])
          pce[g]=100*n_track[[k]][true.mtd[g,1],true.mtd[g,2]]/sum(n)
        }
        PCS[k]=sum(correct)
        PCE[k]=sum(pce)
      }else{
        PCS[k]=0
        PCE[k]=0
      }
      
      if(dim(true.mtds)[1]!=0){
        corrects=pae=vector()
        for(g in 1:dim(true.mtds)[1]){
          corrects[g]=100*(d$mtd[1]==true.mtds[g,1] & d$mtd[2]==true.mtds[g,2])
          pae[g]=100*n_track[[k]][true.mtds[g,1],true.mtds[g,2]]/sum(n)
        }
        PAS[k]=sum(corrects)
        PAE[k]=sum(pae)
      }else{
        PAS[k]=0
        PAE[k]=0
      }
      
      if(dim(true.od)[1]!=0){
        over=vector()
        peover=vector()
        for(g in 1:dim(true.od)[1]){
          over[g]=100*(d$mtd[1]==true.od[g,1] & d$mtd[2]==true.od[g,2])
          peover[g]=100*n_track[[k]][true.od[g,1],true.od[g,2]]/npat[k]
        }
        OVER[k]=sum(over)
        PEOVER[k]=sum(peover)
      }else{
        OVER[k]=0
        PEOVER[k]=0
      }
      
      ZERO[k]=(d$mtd[1]==0 & d$mtd[2]==0)
      
    }
  }
  return(list(n=n_track,y=y_track,mtd=mtd_track,PCS=mean(PCS),PAS=mean(PAS),PCE=mean(PCE),PAE=mean(PAE),
              OVER=mean(OVER),ntox=mean(ntox),npat=mean(npat),ZERO=100*mean(ZERO),PEOVER=mean(PEOVER)))
}

############################


true1=as.vector(matrix(nrow=3,ncol=3,c(0.05,0.1,0.15,0.1,0.15,0.2,0.15,0.2,0.3)))
true2=as.vector(matrix(nrow=3,ncol=3,c(0.05,0.1,0.3,0.1,0.2,0.45,0.2,0.3,0.55)))
true3=as.vector(matrix(nrow=3,ncol=3,c(0.15,0.3,0.45,0.3,0.45,0.55,0.45,0.55,0.65)))
true4=as.vector(matrix(nrow=3,ncol=3,c(0.3,0.45,0.5,0.45,0.5,0.55,0.5,0.55,0.6)))
mean.ratio=c(0.95,0.925,0.9,0.875,0.85)
effSS=c(1,2,3,4,5)
sss=1000

set.seed(1110)
sc1=matrix(NA,nrow=length(mean.ratio),ncol=length(effSS))
for(x1 in 1:length(mean.ratio)){
  for(x2 in 1:length(effSS)){
    simulation = sim(m=3, ss=36, true=true1, start.dose=c(1,1), clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational",
                     mean.ratio=mean.ratio[x1], effSS=effSS[x2], no.diag=TRUE, safety.thresh=1, tt=0.3, n.sims=sss)
    sc1[x1,x2] = simulation$PCS
  }
}

set.seed(2220)
sc2=matrix(NA,nrow=length(mean.ratio),ncol=length(effSS))
for(x1 in 1:length(mean.ratio)){
  for(x2 in 1:length(effSS)){
    simulation = sim(m=3, ss=36, true=true2, start.dose=c(1,1), clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational",  
                     mean.ratio=mean.ratio[x1], effSS=effSS[x2], no.diag=TRUE, safety.thresh=1, tt=0.3, n.sims=sss)
    sc2[x1,x2] = simulation$PCS
  }
}


set.seed(3330)
mean.ratio=c(0.95,0.925,0.9,0.875,0.85)
effSS=c(1,2,3,4,5)
sss=1000

sc3=matrix(NA,nrow=length(mean.ratio),ncol=length(effSS))
for(x1 in 1:length(mean.ratio)){
  for(x2 in 1:length(effSS)){
    simulation = sim(m=3, ss=36, true=true3, start.dose=c(1,1), clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational",  
                     mean.ratio=mean.ratio[x1], effSS=effSS[x2], no.diag=TRUE, safety.thresh=1, tt=0.3, n.sims=sss)
    sc3[x1,x2] = simulation$PCS
  }
}


set.seed(4440)
mean.ratio=c(0.95,0.925,0.9,0.875,0.85)
effSS=c(1,2,3,4,5)
sss=1000

sc4=matrix(NA,nrow=length(mean.ratio),ncol=length(effSS))
for(x1 in 1:length(mean.ratio)){
  for(x2 in 1:length(effSS)){
    simulation = sim(m=3, ss=36, true=true4, start.dose=c(1,1), clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational",  
                     mean.ratio=mean.ratio[x1], effSS=effSS[x2], no.diag=TRUE, safety.thresh=1, tt=0.3, n.sims=sss)
    sc4[x1,x2] = simulation$PCS
  }
}


##finding best parameter values

geom=function(x){
  n=length(x)
  mean=prod(x[1:n])^(1/n)
  return(mean)
}

gm=mm=arith=matrix(NA,nrow=5,ncol=5)
for(i in 1:5){
  for (j in 1:5){
    gm[i,j]=geom(c(sc1[i,j],sc2[i,j],sc3[i,j],sc4[i,j]))
    mm[i,j]=min(c(sc1[i,j],sc2[i,j],sc3[i,j],sc4[i,j]))
    arith[i,j]=mean(c(sc1[i,j],sc2[i,j],sc3[i,j],sc4[i,j]))
  }
}
round(gm,1)
round(mm,1)
which(g>58.5,arr.ind=T)
which(mm>52.5,arr.ind=T)

library(ggplot2)
X = c("0.95","0.925","0.90","0.875","0.85")
Y = c("1","2","3","4","5")
data <- expand.grid(X=X, Y=Y)
data$GM <- as.vector(gm)
data$MM <- as.vector(mm)

pdf(file="fig_SFD_prior_final.pdf",width=4,height=3)
ggplot(data,aes(x=X,y=Y,fill=GM)) + geom_tile() + scale_fill_gradient(low="white", high="darkblue",limits=c(54,60)) +
  labs(y="Prior effective sample size",x="Prior mean ratio",fill="Mean \nPCS (%)") + theme(legend.key.height=unit(1, "cm")) +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))
dev.off()


################
## Calibrating the safety constraint


set.seed(2222)

#take best values
mean.ratio=0.875
effSS=4
#epsilon=seq(from=1,to=0.60,by=-0.02)
sss=1000

#scenarios
s8=matrix(c(0.05,0.10,0.30,0.10,0.20,0.45,0.20,0.30,0.55),nrow=3,ncol=3)
s10=matrix(c(0.15,0.30,0.45,0.30,0.45,0.55,0.45,0.55,0.65),nrow=3,ncol=3)
s13=matrix(c(0.30,0.45,0.50,0.45,0.50,0.55,0.50,0.55,0.60),nrow=3,ncol=3)
s14=matrix(c(0.45,0.50,0.55,0.50,0.55,0.60,0.55,0.60,0.65),nrow=3,ncol=3)


#range of epsilon
epsilon=seq(from=1,to=0.4,by=-0.05)


SFtox_sc14=vector()
SFtox_sc14b=vector()
for(y1 in 1:length(epsilon)){
  simulation = sim(m=3, ss=36, true=s14, start.dose=c(1,1), clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational",  
                   mean.ratio=mean.ratio, effSS=effSS, no.diag=TRUE, safety.thresh=epsilon[y1], tt=0.3, n.sims=sss)
  SFtox_sc14[y1] = simulation$ZERO
  SFtox_sc14b[y1] = mean(sapply(c(1:sss),function(x) sum(simulation$n[[x]]*as.numeric(s14>0.3))))
}

SFtox_sc8=vector()
SFtox_sc8b=vector()
for(y1 in 1:length(epsilon)){
  simulation = sim(m=3, ss=36, true=s8, start.dose=c(1,1), clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational",
                   mean.ratio=mean.ratio, effSS=effSS, no.diag=TRUE, safety.thresh=epsilon[y1], tt=0.3, n.sims=sss)
  SFtox_sc8[y1] = simulation$PCS
  SFtox_sc8b[y1] = mean(sapply(c(1:sss),function(x) sum(simulation$n[[x]]*as.numeric(s8>0.3))))
}

SFtox_sc10=vector()
SFtox_sc10b=vector()
for(y1 in 1:length(epsilon)){
  simulation = sim(m=3, ss=36, true=s10, start.dose=c(1,1), clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational",
                   mean.ratio=mean.ratio, effSS=effSS, no.diag=TRUE, safety.thresh=epsilon[y1], tt=0.3, n.sims=sss)
  SFtox_sc10[y1] = simulation$PCS
  SFtox_sc10b[y1] = mean(sapply(c(1:sss),function(x) sum(simulation$n[[x]]*as.numeric(s10>0.3))))
}

SFtox_sc13=vector()
SFtox_sc13b=vector()
for(y1 in 1:length(epsilon)){
  simulation = sim(m=3, ss=36, true=s13, start.dose=c(1,1), clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational",
                   mean.ratio=mean.ratio, effSS=effSS, no.diag=TRUE, safety.thresh=epsilon[y1], tt=0.3, n.sims=sss)
  SFtox_sc13[y1] = simulation$PCS
  SFtox_sc13b[y1] = mean(sapply(c(1:sss),function(x) sum(simulation$n[[x]]*as.numeric(s13>0.3))))
}


#graphs
SF_plot8a<-data.frame(epsilon,SFtox_sc8b,"Sc 8")
names(SF_plot8a)<-c("epsilon","value","Scenario")
SF_plot10a<-data.frame(epsilon,SFtox_sc10b,"Sc 10")
names(SF_plot10a)<-c("epsilon","value","Scenario")
SF_plot13a<-data.frame(epsilon,SFtox_sc13b,"Sc 13")
names(SF_plot13a)<-c("epsilon","value","Scenario")
SF_plot14a<-data.frame(epsilon,SFtox_sc14b,"Sc 14")
names(SF_plot14a)<-c("epsilon","value","Scenario")
SF_plotGa<-data.frame(epsilon,(SFtox_sc8b*SFtox_sc10b*SFtox_sc13b*SFtox_sc14b)^0.25,"Geo mean")
names(SF_plotGa)<-c("epsilon","value","Scenario")

SF_plota<-rbind(SF_plot8a,SF_plot10a,SF_plot13a,SF_plot14a,SF_plotGa)
names(SF_plota)<-c("epsilon","value","Scenario")

ggplot(SF_plota) + geom_line(aes(x=epsilon, y=value, colour=Scenario)) +
  scale_colour_manual(values=c("red","green","blue","purple","black")) +
  ylim(0,36.01) +
  labs(x="Epsilon",y="Patients Treated on Over-toxic Dose (mean)",color=NULL) 




SF_plot8<-data.frame(epsilon,SFtox_sc8,"Sc 8")
names(SF_plot8)<-c("epsilon","value","Scenario")
SF_plot10<-data.frame(epsilon,SFtox_sc10,"Sc 10")
names(SF_plot10)<-c("epsilon","value","Scenario")
SF_plot13<-data.frame(epsilon,SFtox_sc13,"Sc 13")
names(SF_plot13)<-c("epsilon","value","Scenario")
SF_plot14<-data.frame(epsilon,SFtox_sc14,"Sc 14")
names(SF_plot14)<-c("epsilon","value","Scenario")
SF_plotG<-data.frame(epsilon,(SFtox_sc8*SFtox_sc10*SFtox_sc13*SFtox_sc14)^0.25,"Geo mean")
names(SF_plotG)<-c("epsilon","value","Scenario")
SF_plot<-rbind(SF_plot8,SF_plot10,SF_plot13,SF_plot14,SF_plotG)
names(SF_plot)<-c("epsilon","value","Scenario")

ggplot(SF_plot) + geom_line(aes(x=epsilon, y=value, colour=Scenario)) +
  scale_colour_manual(values=c("red","green","blue","purple","black")) +
  ylim(0,100) +
  labs(x="Epsilon",y="Correct Outcome (%)",color=NULL)



