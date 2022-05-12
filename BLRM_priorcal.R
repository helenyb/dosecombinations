
##BLRM PRIOR CALIBRATION 
library(rjags)
library(mvtnorm)


#the gibbsa sampler in rjags is computationally intensive so this code was run on the university's cluster
#24 cores are used in parallel for efficiency
library("doParallel")
registerDoParallel(cores=24)



##MODEL FILE FOR RJAGS
cat("model{
    
  ### Define likelihood
  for(i in 1:ndose[1]){
    A[i] = doseA[i]/dosestar[1]
    odds1[i] = alpha[1] * pow(A[i],beta[1])
  }
  for(j in 1:ndose[2]){
    B[j] = doseB[j]/dosestar[2]
    odds2[j] = alpha[2] * pow(B[j],beta[2])
  }
  for(i in 1:ndose[1]){ 
    for(j in 1:ndose[2]){
      odds12_null[i,j] = odds1[i] + odds2[j] + odds1[i]*odds2[j]
      odds12[i,j] = odds12_null[i,j] * exp(eta*A[i]*B[j])
  
      n[i,j] ~ dbin(pi[i,j],y[i,j])
      logit(pi[i,j]) = log(odds12[i,j])
    }
  }
  
  ### Define the priors
  parA[1:2] ~ dmnorm.vcov(mu0A[],S2A[,])
  parB[1:2] ~ dmnorm.vcov(mu0B[],S2B[,])
  log_alpha = c(parA[1],parB[1])
  alpha = c(exp(log_alpha[1]),exp(log_alpha[2]))
  log_beta = c(parA[2],parB[2])
  beta = c(exp(log_beta[1]),exp(log_beta[2]))
  eta ~ dnorm(mu0C,1/S2C)
  
  }", file="model_blrm_2.txt")

#FUNCTION TO ESTIMATE P(DLT) FROM MODEL PARAMETERS
pi_hat=function(currentdose,doseA,doseB,dosestar,
                meanp1,meanp2,meanp3,meanp4,meanp5){
  A=doseA[currentdose[1]]/dosestar[1]
  B=doseB[currentdose[2]]/dosestar[2]
  odds1_mod=exp(meanp2)*(A^exp(meanp4))
  odds2_mod=exp(meanp3)*(B^exp(meanp5))
  odds12_mod=(odds1_mod + odds2_mod + odds1_mod*odds2_mod) * exp(meanp1*A*B)
  prob_mod=odds12_mod/(1+odds12_mod)
  return(prob_mod)
}

####################################################################
####################################################################
##SMALL FUNCTIONS 
EWOC=function(x){
  mean(x>0.33)
}
interval=function(x){
  mean(x<0.33 & x>0.16)
}
## FUNCTION TO EXECUTE THE BLRM 
## FUNCTION TO EXECUTE THE BLRM 
BLRM=function(n,y,tt,currentdose,mean1,mean2,var1,var2,var3,safe){
  #dose values and reference dose
  doseA=c(100,200,300)
  doseB=c(100,200,300)
  ndoseA=ndoseB=3
  dosestar=c(max(doseA),max(doseB)) 
  ndose=c(length(doseA),length(doseB))
  #prior parameter and data setup for gibbs sampler 
  mu0A=c(mean1,mean2) #2 params
  mu0B=c(mean1,mean2) #2 params
  mu0C=0 #1 param
  corr=0 #1 param
  S2A=matrix(c(var1,corr,corr,var2),nrow=2) #2 params
  S2B=matrix(c(var1,corr,corr,var2),nrow=2) #2 params
  S2C=var3                                  #1 param
  TauA=matrix(c(var1^(-2),corr,corr,var2^(-2)),nrow=2,ncol=2)
  TauB=matrix(c(var1^(-2),corr,corr,var2^(-2)),nrow=2,ncol=2)
  TauC=var3^(-2)
  inits=list(parA=c(rmvnorm(1,mu0A,S2A)),parB=c(rmvnorm(1,mu0B,S2B)),eta=rnorm(1,mu0C,sqrt(S2C)))
  
  n.patient=matrix(n,nrow=3,byrow=F)
  n.dlt=matrix(y,nrow=3,byrow=F)
  data=list("y"=n.patient,"n"=n.dlt,"doseA"=doseA,"doseB"=doseB,"S2A"=S2A,"S2B"=S2B,"S2C"=S2C,
            "mu0A"=mu0A,"mu0B"=mu0B,"mu0C"=mu0C,"ndose"=ndose,"dosestar"=dosestar,
            "pi"=matrix(NA,nrow=ndose[1],ncol=ndose[2]))
  
  #gibbs sampler
  jags.m=jags.model(file="model_blrm_2.txt",data=data,inits=inits,n.chains=1,n.adapt=1000,quiet=TRUE)
  params=c("eta","log_alpha[1]","log_alpha[2]","log_beta[1]","log_beta[2]")
  js=jags.samples(jags.m,params,n.iter=4000,type="trace",progress.bar="none")
  #organize output
  eta=as.vector(js[[1]])
  log_alpha1=as.vector(js[[2]])
  log_alpha2=as.vector(js[[3]])
  log_beta1=as.vector(js[[4]])
  log_beta2=as.vector(js[[5]])
  
  #find admissible doses
  potentialdose=rbind(currentdose,currentdose-c(1,0),currentdose-c(0,1),currentdose+c(1,0),
                      currentdose+c(0,1),currentdose-c(1,1),currentdose+c(-1,1),currentdose-c(1,-1))
  potentialdose=potentialdose[which(potentialdose[,1]!=0),]
  potentialdose=potentialdose[which(potentialdose[,2]!=0),]
  potentialdose=potentialdose[which(potentialdose[,1]<=ndoseA),]
  potentialdose=potentialdose[which(potentialdose[,2]<=ndoseB),]
  potentialdose=matrix(potentialdose,nrow=dim(potentialdose)[1])
  
  #EWOC procedure to eliminate overly toxic doses
  new_prob=apply(potentialdose,1,pi_hat,doseA=doseA,doseB=doseB,dosestar=dosestar,
                 meanp1=eta,meanp2=log_alpha1,meanp3=log_alpha2,meanp4=log_beta1,meanp5=log_beta2)
  ewoc=apply(new_prob,2,EWOC)
  
  #choose the one that is most likely to be in the acceptable toxicity interval (0.16, 0.33), if none exist we stop early
  if(min(ewoc)<safe & !is.na(min(ewoc))){
    potentialdose=matrix(potentialdose[which(ewoc<safe),],ncol=2)
    new_prob=matrix(new_prob[,which(ewoc<safe)],nrow=4000)
    prob_int=apply(new_prob,2,interval)
    choice=which(prob_int==max(prob_int))
    rule=c(potentialdose[choice,1],potentialdose[choice,2])
    abandon=0
  }else{
    abandon=1
    rule=c(0,0)
  }
  mtd=c(0,0)
  
  #prob that toxicity of each dose lies in the interval (0.16, 0.33)
  if(is.na(sum(n)>=36 & abandon==0)){
    mtd=c(0,0)
    abandon=1
    rule=c(0,0)
  }else{
    if(sum(n)>=36 & abandon==0){
      int=matrix(0,nrow=ndoseA,ncol=ndoseB)
      for(i in 1:ndoseA){
        for(j in 1:ndoseB){
          if(n.patient[i,j]>=6){
            probs = pi_hat(c(i,j),doseA=doseA,doseB=doseB,dosestar=dosestar,
                           meanp1=eta,meanp2=log_alpha1,meanp3=log_alpha2,meanp4=log_beta1,meanp5=log_beta2)
            int[i,j] = interval(probs)
          }
        }
      }
      mtd=which(int==max(int),arr.ind=TRUE)
      if(dim(mtd)[1]>1){
        mtd=mtd[1,]
      }
    }
  }
  
  return(list(rule=rule,mtd=mtd,abandon=abandon))
}
################################################

# SIMULATION FUNCTION
SIMpar=function(k,ss=36,m=3,tt=0.3,true,mean1,mean2,var1,var2,var3,safe){
  set.seed(k)
  ndoseA=ndoseB=3
  start.dose=c(1,1)
  #define truths based on scenario spec
  trueMAT=matrix(true,nrow=ndoseA,ncol=ndoseB)
  true.mtd=which(trueMAT==tt,arr.ind=T)
  true.mtds=which(trueMAT>0.16 & trueMAT<0.33,arr.ind=T)
  true.od=which(trueMAT>tt,arr.ind=T)
  #true.ud=which(trueMAT<tt,arr.ind=T)
  
  #initializing
  n_track=list()
  y_track=list()
  mtd_track=list()
  
  PCS=vector()
  PAS=vector()
  PCE=vector()
  PAE=vector()
  OVER=vector()
  ntox=vector()
  npat=vector()
  ZERO=vector()
  
  
  n=y=rep(0,ndoseA*ndoseB)
  rule=start.dose
  new=start.dose[1]+ndoseA*(start.dose[2]-1)
  abandon=0
  i=0
  
  #update for each cohort
  while(abandon==0 & i<ss){
    i=i+m
    currentdose=rule
    n[new]=n[new]+m
    dlt=rbinom(1,size=m,prob=c(true[new]))
    y[new]=y[new]+dlt
    d=BLRM(n,y,tt,currentdose,mean1,mean2,var1,var2,var3,safe)
    abandon=d$abandon
    if(abandon==0){
      rule=d$rule
      new=rule[1]+ndoseA*(rule[2]-1) 
    }
  }
  
  
  #assigning values for output
  npat=sum(n)
  ntox=sum(y)
  n_track=matrix(n,nrow=ndoseA,ncol=ndoseB)
  y_track=matrix(y,nrow=ndoseA,ncol=ndoseB)
  mtd_track=d$mtd
  
  if(dim(true.mtd)[1]!=0){
    correct=pce=vector()
    for(g in 1:dim(true.mtd)[1]){
      correct[g]=100*(d$mtd[1]==true.mtd[g,1] & d$mtd[2]==true.mtd[g,2])
      pce[g]=100*n_track[true.mtd[g,1],true.mtd[g,2]]/sum(n)
    }
    PCS=sum(correct)
    PCE=sum(pce)
  }else{
    PCS=0
    PCE=0
  }
  
  if(dim(true.mtds)[1]!=0){
    corrects=pae=vector()
    for(g in 1:dim(true.mtds)[1]){
      corrects[g]=100*(d$mtd[1]==true.mtds[g,1] & d$mtd[2]==true.mtds[g,2])
      pae[g]=100*n_track[true.mtds[g,1],true.mtds[g,2]]/sum(n)
    }
    PAS=sum(corrects)
    PAE=sum(pae)
  }else{
    PAS=0
    PAE=0
  }
  
  if(dim(true.od)[1]!=0){
    nover= over=vector()
    for(g in 1:dim(true.od)[1]){
      over[g]=100*(d$mtd[1]==true.od[g,1] & d$mtd[2]==true.od[g,2])
      nover[g]=n_track[true.od[g,1],true.od[g,2]]
    }
    NOVER=sum(nover)
    OVER=sum(over)
  }else{
    OVER=0
    NOVER=0
  }
  
  ZERO=(d$mtd[1]==0 & d$mtd[2]==0)
  
  return(c(PCS,NOVER))
}

####################################################################################################################
#SCENARIOS FOR CALIBRATION
set.seed(1000)
true1=matrix(nrow=3,ncol=3,c(0.05,0.1,0.15,0.1,0.15,0.2,0.15,0.2,0.3))
true2=matrix(nrow=3,ncol=3,c(0.05,0.1,0.3,0.1,0.2,0.45,0.2,0.3,0.55))
true3=matrix(nrow=3,ncol=3,c(0.15,0.3,0.45,0.3,0.45,0.55,0.45,0.55,0.65))
true4=matrix(nrow=3,ncol=3,c(0.3,0.45,0.5,0.45,0.5,0.55,0.5,0.55,0.6))


#GRID FOR PARAMETERS (These were changed and rerun in calibrations in order to efficiently find best values)
c1=c(-1,-1.5,-2.5)
c2=0
v1=c(1,2,3)
v2=c(0.25,0.5,1)
v3=c(0.03)
nsims=500 #only 500 simulation for calibration to ensure it is done in reasonable time


#arrays to fill with results
sc1=array(NA,dim=c(length(c1),length(c2),length(v1),length(v2),length(v3)))
sc2=array(NA,dim=c(length(c1),length(c2),length(v1),length(v2),length(v3)))
sc3=array(NA,dim=c(length(c1),length(c2),length(v1),length(v2),length(v3)))
sc4=array(NA,dim=c(length(c1),length(c2),length(v1),length(v2),length(v3)))

sc1o=array(NA,dim=c(length(c1),length(c2),length(v1),length(v2),length(v3)))
sc2o=array(NA,dim=c(length(c1),length(c2),length(v1),length(v2),length(v3)))
sc3o=array(NA,dim=c(length(c1),length(c2),length(v1),length(v2),length(v3)))
sc4o=array(NA,dim=c(length(c1),length(c2),length(v1),length(v2),length(v3)))


#simulations for all combinations of parameters across the four scenarios 
for(x1 in 1:length(c1)){
  for(x2 in 1:length(c2)){
    for(x3 in 1:length(v1)){
      for(x4 in 1:length(v2)){
        for(x5 in 1:length(v3)){
          
          meanvec<-c(NA,NA)
          try(
            meanvec<-colMeans(foreach(i=1:nsims, .combine=rbind ) %dopar% {
              SIMpar(i,true=true1,mean1=c1[x1],mean2=c2[x2],var1=v1[x3],var2=v2[x4],var3=v3[x5],safe=0.25)
            })
          )
          
          sc1[x1,x2,x3,x4,x5]<-meanvec[1] 
          sc1o[x1,x2,x3,x4,x5]<-meanvec[2]
          
          save.image("BLRM_prior.RData")
          
          meanvec<-c(NA,NA)
          try(
            meanvec<-colMeans(
              foreach(i=1:nsims, .combine=rbind ) %dopar% {
                SIMpar(i,true=true2,mean1=c1[x1],mean2=c2[x2],var1=v1[x3],var2=v2[x4],var3=v3[x5],safe=0.25)
              }
            ))
          sc2[x1,x2,x3,x4,x5]<-meanvec[1]
          sc2o[x1,x2,x3,x4,x5]<-meanvec[2]
          
          save.image("BLRM_prior.RData")
          
          meanvec<-c(NA,NA)
          try(
            meanvec<-colMeans(foreach(i=1:nsims, .combine=rbind ) %dopar% {
              SIMpar(i,true=true3,mean1=c1[x1],mean2=c2[x2],var1=v1[x3],var2=v2[x4],var3=v3[x5],safe=0.25)
            }))
          sc3[x1,x2,x3,x4,x5]<-meanvec[1]
          sc3o[x1,x2,x3,x4,x5]<-meanvec[2]
          
          save.image("BLRM_prior.RData")
          
          meanvec<-c(NA,NA)
          try(
            meanvec<-colMeans(foreach(i=1:nsims, .combine=rbind ) %dopar% {
              SIMpar(i,true=true4,mean1=c1[x1],mean2=c2[x2],var1=v1[x3],var2=v2[x4],var3=v3[x5],safe=0.25)
            }))
          sc4[x1,x2,x3,x4,x5]<-meanvec[1]
          sc4o[x1,x2,x3,x4,x5]<-meanvec[2]
          
          save.image("BLRM_prior.RData")
          
          
        } 
      }
    }
  }
}
