##SF SIMS

library(doParallel)
library(rjags)

#the gibbsa sampler in rjags is computationally intensive so this code was run on the university's cluster
#24 cores are used in parallel for efficiency
registerDoParallel(cores=24)


#scenarios
s1=matrix(c(0.05,0.10,0.15,0.10,0.15,0.20,0.15,0.20,0.30),nrow=3,ncol=3)
s2=matrix(c(0.05,0.10,0.20,0.10,0.20,0.30,0.15,0.30,0.45),nrow=3,ncol=3)
s3=matrix(c(0.02,0.10,0.20,0.05,0.15,0.30,0.10,0.20,0.45),nrow=3,ncol=3)
s4=matrix(c(0.05,0.10,0.20,0.10,0.20,0.45,0.15,0.30,0.60),nrow=3,ncol=3)
s5=matrix(c(0.02,0.20,0.45,0.05,0.30,0.55,0.15,0.45,0.65),nrow=3,ncol=3)
s6=matrix(c(0.10,0.15,0.30,0.15,0.30,0.45,0.30,0.45,0.60),nrow=3,ncol=3)
s7=matrix(c(0.10,0.15,0.30,0.20,0.30,0.50,0.45,0.50,0.60),nrow=3,ncol=3)
s8=matrix(c(0.05,0.10,0.30,0.10,0.20,0.45,0.20,0.30,0.55),nrow=3,ncol=3)
s9=matrix(c(0.10,0.30,0.40,0.15,0.40,0.50,0.30,0.50,0.60),nrow=3,ncol=3)
s10=matrix(c(0.15,0.30,0.45,0.30,0.45,0.55,0.45,0.55,0.65),nrow=3,ncol=3)
s11=matrix(c(0.02,0.30,0.45,0.05,0.45,0.60,0.10,0.60,0.75),nrow=3,ncol=3)
s12=matrix(c(0.20,0.45,0.65,0.30,0.50,0.70,0.45,0.55,0.75),nrow=3,ncol=3)
s13=matrix(c(0.30,0.45,0.50,0.45,0.50,0.55,0.50,0.55,0.60),nrow=3,ncol=3)
s14=matrix(c(0.45,0.50,0.55,0.50,0.55,0.60,0.55,0.60,0.65),nrow=3,ncol=3)
s15=matrix(c(0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10),nrow=3,ncol=3)
scenario=list(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15)

#construct monotherapy prior
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
#construct operational prior given parameters
op_prior=function(mean.ratio,effSS,space){
  ndoseA=space[1];ndoseB=space[2]
  thetahat=rep(mean.ratio,space[1])
  tauhat=rep(mean.ratio,(space[2]-1))
  a=effSS*thetahat;b=effSS*(1-thetahat)
  e=effSS*tauhat;f=effSS*(1-tauhat)
  
  return(list(thetahat=thetahat,tauhat=tauhat,a=a,b=b,e=e,f=f))
}
#function for outputting p(dlt) using model parameter values
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
##function for elimination 
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


#model file for rjags
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


#function to execute the SF-design 
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
  js=jags.samples(jags.m,params,n.iter=2000,type="trace")
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
    #mtd=which(abs(tox-tt)==min(abs(tox-tt)),arr.ind=TRUE)
    mtd=which(abs(tox-tt)==min(abs(tox-tt)),arr.ind=TRUE)
  }
  
  return(list(prob.overdose=prob.overdose,n=n,x=x,tox=tox,rule=rule,mtd=mtd))
}

#function to conduct n.sims simulations
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
  
 # for(k in 1:n.sims){
    
    ind_trial<-function(k){
      set.seed(k)
      k<-1 # dummy variable for parallel sims
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
  
  return(list(n=n_track[[1]],y=y_track[[1]],mtd=mtd_track[[1]],PCS=mean(PCS),PAS=mean(PAS),PCE=mean(PCE),PAE=mean(PAE),
              OVER=mean(OVER),ntox=mean(ntox),npat=mean(npat),ZERO=100*mean(ZERO),PEOVER=mean(PEOVER)))
    }
    
    
    res_list<-foreach(i=1:n.sims, combine=list ) %dopar% {
    ind_trial(i)
    }
   
    return(list(n=lapply(res_list,"[[",'n'),y=lapply(res_list,"[[",'y'),mtd=lapply(res_list,"[[",'mtd'),
                PCS=mean(unlist(lapply(res_list,"[[",'PCS'))),PAS=mean(unlist(lapply(res_list,"[[",'PAS'))),
                PCE=mean(unlist(lapply(res_list,"[[",'PCE'))),PAE=mean(unlist(lapply(res_list,"[[",'PAE'))),
                OVER=mean(unlist(lapply(res_list,"[[",'OVER'))),ntox=mean(unlist(lapply(res_list,"[[",'ntox'))),
           npat=mean(unlist(lapply(res_list,"[[",'npat'))),ZERO=mean(unlist(lapply(res_list,"[[",'ZERO'))),
           PEOVER=mean(unlist(lapply(res_list,"[[",'PEOVER')))))
    
}
####################


# SIMULATIONS -------------------------------------------------------------

#conduct simulations for the 15 scenarios
nsims=2000

PCS=PAS=PCE=PAE=OVER=PEOVER=ntox=npat=ZERO=rep(NA,15)
EXTRA=SEL=EXP=list()
for(x1 in 1:15){
  trial=sim(m=3, ss=36, true=scenario[[x1]], start.dose=c(1,1), clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational",  
            mean.ratio=0.825, effSS=4, no.diag=TRUE, safety.thresh=0.65, tt=0.3, n.sims=nsims)
  
  #PCS, PAS, PCE
  PCS[x1]=trial$PCS
  PAS[x1]=trial$PAS
  PCE[x1]=trial$PCE
  PAE[x1]=trial$PAE
  
  #Overtoxic
  OVER[x1]=trial$OVER
  PEOVER[x1]=trial$PEOVER
  ZERO[x1]=trial$ZERO
  
  #Number of patients and DLTs
  ntox[x1]=trial$ntox
  npat[x1]=trial$npat
  
  #extra outputs 
  EXTRA[[x1]]=trial$mtd
  SELE=matrix(NA,nrow=3,ncol=3)
  for(i in 1:3){
    for(j in 1:3){
      sel=rep(NA,nsims)
      for(k in 1:nsims){
        if(sum(trial$mtd[[k]])==0){
          sel[k]=0
        }else{
          sel[k]=100*(trial$mtd[[k]][1,1]==i & trial$mtd[[k]][1,2]==j)
        }
      }
      SELE[i,j]=mean(sel)
    }
  }
  SEL[[x1]]=SELE
  
  EXPE=matrix(NA,nrow=3,ncol=3)
  for(i in 1:3){
    for(j in 1:3){
      exp=rep(NA,nsims)
      for(k in 1:nsims){
        exp[k]=(trial$n[[k]][i,j])
      }
      EXPE[i,j]=mean(exp)
    }
  }
  EXP[[x1]]=EXPE
}
#label output
PCS1=PCS
PAS1=PAS
PCE1=PCE
PAE1=PAE
OVER1=OVER
PEOVER1=PEOVER
ntox1=ntox
npat1=npat
ZERO1=ZERO
EXTRA1=EXTRA
SEL1=SEL
EXP1=EXP
NOVER1<-unlist(lapply(c(1:15), function(x) sum(EXP1[[x]][scenario[[x]]>0.33])))


save.image("SF_sims.RData")
