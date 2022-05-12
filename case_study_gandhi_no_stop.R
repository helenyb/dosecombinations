#case study example (gandhi)
##no stopping rules

#data from real trial
true_n<-matrix(c(4,5,4,4,5,6,8,2,NA),byrow = T,nrow=3)
true_y<-matrix(c(0,1,0,1,0,3,1,1,NA),byrow = T,nrow=3)

#sample/cohort size
tSS<-36
co_size<-3

#simulate all potential patient responses
patients<-array(NA,dim=c(3,3,tSS))
set.seed(100)
for(i in 1:3){
  for (j in 1:3){
    if(!is.na(true_n[i,j])){
      patients[i,j,c(1:true_n[i,j])]<-sample(c(rep(1,true_y[i,j]),rep(0,true_n[i,j]-true_y[i,j])),size=true_n[i,j],replace=FALSE)
      patients[i,j,c((true_n[i,j]+1):tSS)]<-generate_patients(npatients = tSS-true_n[i,j],Nobs=true_n[i,j],Xobs = true_y[i,j])
      
    }
  }
}

for(i in 1:3){
  for (j in 1:3){
    if(is.na(true_n[i,j])){
      patients[i,j,]<-generate_patients(npatients = tSS,Nobs=4,Xobs =2)
      
    }
  }
}

################################################################################
#boin
#initilaize
current_n<-matrix(c(co_size,rep(0,8)),byrow=T,nrow=3)
current_y<-matrix(rep(0,9),byrow=T,nrow=3)

#define parameters
target=0.3
a1=0.65
a2=1.4
max_cohort=tSS/co_size
cohort<-1

#empty results matrices
dose_esc<-matrix(NA,nrow=max_cohort,ncol=2)
dose_esc2<-matrix(NA,nrow=max_cohort,ncol=3)
dose_esc[1,]<-c(1,1)
dose_esc2[1,]<-c(1,1,0)
set.seed(2)
while(cohort<max_cohort){
  #choose next combination
  nextdose<-next.comb(target=target, npts=current_n, ntox=current_y, dose.curr=dose_esc[cohort,], n.earlystop=100,
                      p.saf=a1*target, p.tox=a2*target, cutoff.eli=0.84,
                      extrasafe=FALSE, offset=0.05)$next_dc
  dose_esc[cohort+1,]<-nextdose
  
  #update counts
  co_y<-sum(patients[nextdose[1],nextdose[2],c((current_n[nextdose[1],nextdose[2]]+1):(current_n[nextdose[1],nextdose[2]]+co_size))])
  current_n[nextdose[1],nextdose[2]]<-current_n[nextdose[1],nextdose[2]]+co_size
  current_y[nextdose[1],nextdose[2]]<-current_y[nextdose[1],nextdose[2]]+co_y
  dose_esc2[cohort+1,]<-c(nextdose,co_y)
  cohort<-cohort+1
  
  
  #recommended dose
  if(cohort==max_cohort){
    recdose<-select.mtd.comb(target=target, npts=current_n, ntox=current_y,
                             cutoff.eli=0.84,
                             extrasafe=FALSE, offset=0.05)$MTD
  }
  
}
#store
assign(paste(c("BOIN_rec_",tSS),collapse=""),recdose)
 assign(paste(c("BOIN_dose_esc_",tSS),collapse=""),dose_esc)
assign(paste(c("BOIN_dose_esc2_",tSS),collapse=""),dose_esc2)
 assign(paste(c("BOIN_nmat_",tSS),collapse=""),current_n)
assign(paste(c("BOIN_ymat_",tSS),collapse=""),current_y)

################################################################################
#pipe
#function to calculate prior median matrix for any size
prior.m=function(low,inc,rows,cols){
  t1=low
  t2=low+inc
  t3=low+2*inc
  t4=low+3*inc
  t5=low+4*inc
  prior.med=matrix(NA,nrow=3,ncol=3)
  
  
  for(i in 1:rows){
    for(j in 1:cols){
      prior.med[i,j]<-low+(i-1)*(j-1)*inc
    }
  }
  return(prior.med)
}

#function to reformat data into foramt for pipe.design
create_dataframe<-function(nmat,patients){
  data_out<-data.frame()
  pat_num<-1
  for (i in 1:3){
    for( j in 1:3){
      
      
      nij<-nmat[i,j]
      if((!is.na(nij))&(nij!=0)){
        patij<-patients[i,j,]
        #   browser()
        for(pat in 1:nij){
          pat_data<-c(pat_num,i,j,patij[pat])
          pat_num<-pat_num+1
          data_out<-rbind(data_out,pat_data)
        }
      }
    }
  }
  
  names(data_out)<-c("patient","doseA","doseB","tox")
  return(data_out)
}

#initilaize
current_n<-matrix(c(co_size,rep(0,8)),byrow=T,nrow=3)
current_y<-matrix(rep(0,9),byrow=T,nrow=3)

target=0.3
cohort<-1

#empty results matrices
dose_esc<-matrix(NA,nrow=max_cohort,ncol=2)
dose_esc2<-matrix(NA,nrow=max_cohort,ncol=3)
dose_esc[1,]<-c(1,1)
dose_esc2[1,]<-c(1,1,0)
set.seed(2)

while(cohort<max_cohort){
  #manipulate current data into format for pipe
  current_data<-create_dataframe(nmat=current_n,patients = patients)
  #choose next combination
  nextdoseall<-pipe.design(N=sum(current_n)+co_size, S=1, c=co_size,theta=target,
                           prior.med =prior.m(0.05, 0.025,3,3) ,prior.ss = matrix(nrow=3,ncol=3,1/18),
                           epsilon=0.5,data=current_data, strategy="ss-random", admis="closest", 
                           constraint="neighbouring-nodiag", uppertox.constraint=NULL, mode="sim")
  
  nextdose<- c(nextdoseall$rec.i.sim[,cohort+1,],nextdoseall$rec.j.sim[,cohort+1,])
  
  dose_esc[cohort+1,]<-nextdose
  
  #update counts
  co_y<-sum(patients[nextdose[1],nextdose[2],c((current_n[nextdose[1],nextdose[2]]+1):(current_n[nextdose[1],nextdose[2]]+co_size))])
  current_n[nextdose[1],nextdose[2]]<-current_n[nextdose[1],nextdose[2]]+co_size
  current_y[nextdose[1],nextdose[2]]<-current_y[nextdose[1],nextdose[2]]+co_y
  dose_esc2[cohort+1,]<-c(nextdose,co_y)
  cohort<-cohort+1
  
  
  
  #recommended dose
  if(cohort==max_cohort){
    current_data<-create_dataframe(nmat=current_n,patients = patients)
    nextdoseall<-pipe.design(N=sum(current_n), S=1, c=co_size,theta=target,
                             prior.med =prior.m(0.05, 0.025,3,3) ,prior.ss = matrix(nrow=3,ncol=3,1/18),
                             epsilon=0.5,data=current_data, strategy="ss-random", admis="closest", 
                             constraint="neighbouring-nodiag", uppertox.constraint=NULL, mode="sim")
    
    recdose<- which(nextdoseall$rec>0,arr.ind = T)
  }
}
#store

assign(paste(c("PIPE_rec_",tSS),collapse=""),recdose)
assign(paste(c("PIPE_dose_esc_",tSS),collapse=""),dose_esc)
assign(paste(c("PIPE_dose_esc2_",tSS),collapse=""),dose_esc2)
assign(paste(c("PIPE_nmat_",tSS),collapse=""),current_n)
assign(paste(c("PIPE_ymat_",tSS),collapse=""),current_y)


################################################################################
#keyboard
#initilaize
current_n<-matrix(c(co_size,rep(0,8)),byrow=T,nrow=3)
current_y<-matrix(rep(0,9),byrow=T,nrow=3)

target=0.3
cohort<-1

#empty results matrices
dose_esc<-matrix(NA,nrow=max_cohort,ncol=2)
dose_esc2<-matrix(NA,nrow=max_cohort,ncol=3)
dose_esc[1,]<-c(1,1)
dose_esc2[1,]<-c(1,1,0)
set.seed(2)
while(cohort<max_cohort){
  #choose next combination
  nextdose<-next.comb.kb.new(target=target, npts=current_n, ntox=current_y, dose.curr=dose_esc[cohort,], n.earlystop=100,
                             marginL = 0.05, 
                             marginR = 0.05, cutoff.eli=0.84,
                             extrasafe=FALSE, offset=0.05)$next_dc
  
  
  
  dose_esc[cohort+1,]<-nextdose
  
  #update counts
  co_y<-sum(patients[nextdose[1],nextdose[2],c((current_n[nextdose[1],nextdose[2]]+1):(current_n[nextdose[1],nextdose[2]]+co_size))])
  current_n[nextdose[1],nextdose[2]]<-current_n[nextdose[1],nextdose[2]]+co_size
  current_y[nextdose[1],nextdose[2]]<-current_y[nextdose[1],nextdose[2]]+co_y
  dose_esc2[cohort+1,]<-c(nextdose,co_y)
  cohort<-cohort+1
  
  
  #recommended dose
  if(cohort==max_cohort){
    recdose<-select.mtd.comb.kb(target=target, npts=current_n, ntox=current_y,
                                cutoff.eli=0.84,
                                extrasafe=FALSE, offset=0.05)$MTD
  }
  
}
#store
assign(paste(c("KEY_rec_",tSS),collapse=""),recdose)
assign(paste(c("KEY_dose_esc_",tSS),collapse=""),dose_esc)
assign(paste(c("KEY_dose_esc2_",tSS),collapse=""),dose_esc2)
assign(paste(c("KEY_nmat_",tSS),collapse=""),current_n)
assign(paste(c("KEY_ymat_",tSS),collapse=""),current_y)


################################################################################
##sf
library(rjags)

#functions for creating prior for input

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
  ratio=matrix(NA,nrow=5000,ncol=(ndoseA+ndoseB-1))
  js=jags.samples(jags.m,params,n.iter=5000,type="trace",progress.bar="none")
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



#initilaize
current_n<-matrix(c(co_size,rep(0,8)),byrow=T,nrow=3)
current_y<-matrix(rep(0,9),byrow=T,nrow=3)

target=0.3
cohort<-1

#empty results matrices
dose_esc<-matrix(NA,nrow=max_cohort,ncol=2)
dose_esc2<-matrix(NA,nrow=max_cohort,ncol=3)
dose_esc[1,]<-c(1,1)
dose_esc2[1,]<-c(1,1,0)
set.seed(2)

while(cohort<max_cohort){
  #choose next combination
  nextdose<- 
    SFdesign(clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational", mean.ratio=0.825, effSS=4,
             n.patients=current_n, n.dlts=current_y, no.diag=TRUE, safety.thresh=0.65, tt=0.3, currentdose=dose_esc[cohort,])$rule
  
  dose_esc[cohort+1,]<-nextdose
  
  #update counts
  co_y<-sum(patients[nextdose[1],nextdose[2],c((current_n[nextdose[1],nextdose[2]]+1):(current_n[nextdose[1],nextdose[2]]+co_size))])
  current_n[nextdose[1],nextdose[2]]<-current_n[nextdose[1],nextdose[2]]+co_size
  current_y[nextdose[1],nextdose[2]]<-current_y[nextdose[1],nextdose[2]]+co_y
  dose_esc2[cohort+1,]<-c(nextdose,co_y)
  cohort<-cohort+1
  
  
  #recommended dose
  if(cohort==max_cohort){
    recdose<-SFdesign(clinA=NULL, clinB=NULL, space=c(3,3), prior.type="operational", mean.ratio=0.825, effSS=4,
                      n.patients=current_n, n.dlts=current_y, no.diag=TRUE, safety.thresh=0.65, tt=0.3, currentdose=dose_esc[cohort,])$mtd
    
    
  }
  
}
#store
assign(paste(c("SF_rec_",tSS),collapse=""),recdose)
assign(paste(c("SF_dose_esc_",tSS),collapse=""),dose_esc)
assign(paste(c("SF_dose_esc2_",tSS),collapse=""),dose_esc2)
assign(paste(c("SF_nmat_",tSS),collapse=""),current_n)
assign(paste(c("SF_ymat_",tSS),collapse=""),current_y)



################################################################################
#blrm
#model file for rjags
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
##SMALL FUNCTIONS 
EWOC=function(x){
  mean(x>0.33)
}
interval=function(x){
  mean(x<0.33 & x>0.16)
}

## FUNCTION TO EXECUTE THE BLRM 
BLRM=function(n,y,tt,currentdose,mean1,mean2,var1,var2,var3,safe,SS){
  
  doseA=c(120,160,200)
  doseB=c(25,50,75)
  ndoseA=ndoseB=3
  dosestar=c(max(doseA),max(doseB)) 
  ndose=c(length(doseA),length(doseB))
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
  
  jags.m=jags.model(file="model_blrm_2.txt",data=data,inits=inits,n.chains=1,n.adapt=1000,quiet=TRUE)
  params=c("eta","log_alpha[1]","log_alpha[2]","log_beta[1]","log_beta[2]")
  js=jags.samples(jags.m,params,n.iter=4000,type="trace",progress.bar="none")
  
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
    choice=which(prob_int==max(prob_int))[1]
    rule=c(potentialdose[choice,1],potentialdose[choice,2])
    #browser()
    abandon=0
  }else{
    abandon=1
    rule=c(0,0)
  }
  mtd=c(0,0)
  
  #prob that toxicity of each dose lies in the interval (0.16, 0.33)
  if(is.na(sum(n)>=SS & abandon==0)){
    mtd=c(0,0)
    abandon=1
    rule=c(0,0)
  }else{
    if(sum(n)>=SS & abandon==0){
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
library(mvtnorm)
c1=-2.5
c2=0
v1=1
v2=0.5
v3=0.03
ewoc_val=0.25
#initilaize
current_n<-matrix(c(co_size,rep(0,8)),byrow=T,nrow=3)
current_y<-matrix(rep(0,9),byrow=T,nrow=3)
target=0.3
cohort<-1

#empty results matrices
dose_esc<-matrix(NA,nrow=max_cohort,ncol=2)
dose_esc2<-matrix(NA,nrow=max_cohort,ncol=3)
dose_esc[1,]<-c(1,1)
dose_esc2[1,]<-c(1,1,0)
set.seed(2)
while(cohort<max_cohort){
  #choose next combination
  nextdose<-BLRM(n=current_n,y=current_y,tt=target,currentdose=dose_esc[cohort,],mean1=c1,mean2=c2,var1=v1,var2=v2,var3=v3,safe=ewoc_val,SS=co_size*max_cohort)$rule
  
  dose_esc[cohort+1,]<-nextdose
  
  #update counts
  co_y<-sum(patients[nextdose[1],nextdose[2],c((current_n[nextdose[1],nextdose[2]]+1):(current_n[nextdose[1],nextdose[2]]+co_size))])
  current_n[nextdose[1],nextdose[2]]<-current_n[nextdose[1],nextdose[2]]+co_size
  current_y[nextdose[1],nextdose[2]]<-current_y[nextdose[1],nextdose[2]]+co_y
  dose_esc2[cohort+1,]<-c(nextdose,co_y)
  cohort<-cohort+1
  
  
  #recommended dose
  if(cohort==max_cohort){
    recdose<-BLRM(n=current_n,y=current_y,tt=target,currentdose=dose_esc[cohort,],mean1=c1,mean2=c2,var1=v1,var2=v2,var3=v3,safe=ewoc_val,SS=co_size*max_cohort)$mtd
    
  }
}




#store

assign(paste(c("BLRM_rec_",tSS),collapse=""),recdose)
assign(paste(c("BLRM_dose_esc_",tSS),collapse=""),dose_esc)
assign(paste(c("BLRM_dose_esc2_",tSS),collapse=""),dose_esc2)
assign(paste(c("BLRM_nmat_",tSS),collapse=""),current_n)
assign(paste(c("BLRM_ymat_",tSS),collapse=""),current_y)


#alternative prior
c1=log(0.25)
c2=0
v1=4
v2=1
v3=0.03
ewoc_val=0.25

#initilaize
current_n<-matrix(c(co_size,rep(0,8)),byrow=T,nrow=3)
current_y<-matrix(rep(0,9),byrow=T,nrow=3)
target=0.3
cohort<-1

#empty results matrices
dose_esc<-matrix(NA,nrow=max_cohort,ncol=2)
dose_esc2<-matrix(NA,nrow=max_cohort,ncol=3)
dose_esc[1,]<-c(1,1)
dose_esc2[1,]<-c(1,1,0)
set.seed(2)
while(cohort<max_cohort){
  #choose next combination
  nextdose<-BLRM(n=current_n,y=current_y,tt=target,currentdose=dose_esc[cohort,],mean1=c1,mean2=c2,var1=v1,var2=v2,var3=v3,safe=ewoc_val,SS=co_size*max_cohort)$rule
  
  dose_esc[cohort+1,]<-nextdose
  
  #update counts
  co_y<-sum(patients[nextdose[1],nextdose[2],c((current_n[nextdose[1],nextdose[2]]+1):(current_n[nextdose[1],nextdose[2]]+co_size))])
  current_n[nextdose[1],nextdose[2]]<-current_n[nextdose[1],nextdose[2]]+co_size
  current_y[nextdose[1],nextdose[2]]<-current_y[nextdose[1],nextdose[2]]+co_y
  dose_esc2[cohort+1,]<-c(nextdose,co_y)
  cohort<-cohort+1
  
  
  #recommended dose
  if(cohort==max_cohort){
    recdose<-BLRM(n=current_n,y=current_y,tt=target,currentdose=dose_esc[cohort,],mean1=c1,mean2=c2,var1=v1,var2=v2,var3=v3,safe=ewoc_val,SS=co_size*max_cohort)$mtd
    
  }
}


#store

assign(paste(c("BLRM_rec2_",tSS),collapse=""),recdose)
assign(paste(c("BLRM_dose_esc2_",tSS),collapse=""),dose_esc)
assign(paste(c("BLRM_dose_esc22_",tSS),collapse=""),dose_esc2)
assign(paste(c("BLRM_nmat2_",tSS),collapse=""),current_n)
assign(paste(c("BLRM_ymat2_",tSS),collapse=""),current_y)

