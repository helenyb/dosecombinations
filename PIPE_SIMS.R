#PIPE SIMS

library(pipe.design)

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

#function to output prior median matrix for a 3x3 combination space
prior.m=function(low,inc){
  t1=low
  t2=low+inc
  t3=low+2*inc
  t4=low+3*inc
  t5=low+4*inc
  prior.med=matrix(nrow=3,ncol=3,c(t1,t2,t3,t2,t3,t4,t3,t4,t5))
  return(prior.med)
}


#initialize empty lists
PCS=PAS=PCE=PAE=OVER=PEOVER=NOVER=ZERO=ntox=npat=rep(NA,15)
EXTRA=SEL=EXP=list()
#conduct simulations for all 15 scenarios
for(x1 in 1:15){
  trial = pipe.design(N=36,S=2000,c=3, theta=0.30, pi=scenario[[x1]], prior.med=prior.m(0.05, 0.025), prior.ss=matrix(nrow=3,ncol=3,1/18),
                      strategy="ss-random", admis="closest", constraint="neighbouring-nodiag", uppertox.constraint=NULL, epsilon=0.5, mode="sim")
  
  mtd=which(scenario[[x1]]==0.30,arr.ind=T)
  mtds=which((scenario[[x1]]>0.16 & scenario[[x1]]<0.33),arr.ind=T)
  over=which(scenario[[x1]]>0.33,arr.ind=T)
  dose = pickone(trial)
  
  #PCS, PAS, PCE
  if(dim(mtd)[1]!=0){
    pce=vector()
    pae=vector()
    
    same = matrix(NA,nrow=2000,ncol=dim(mtd)[1])
    for(q in 1:2000){
      for(z in 1:dim(mtd)[1]){
        same[q,z] = (mtd[z,1]==dose[q,1] & mtd[z,2]==dose[q,2])
      }
    }
    success=apply(same,1,sum)
    PCS[x1] = 100*mean(success)
    
    extra=vector()
    for(z in 1:dim(mtd)[1]){
      extra[z] = 100*mean(same[,z])
    }
    EXTRA[[x1]] = extra
    
    for(z in 1:dim(mtd)[1]){
      pce[z]=100*trial$exp[mtd[z,1],mtd[z,2]]
    }
    PCE[x1]=sum(pce)
    
    same2 = matrix(NA,nrow=2000,ncol=dim(mtds)[1])
    for(q in 1:2000){
      for(z in 1:dim(mtds)[1]){
        same2[q,z] = (mtds[z,1]==dose[q,1] & mtds[z,2]==dose[q,2])
      }
    }
    success2=apply(same2,1,sum)
    PAS[x1] = 100*mean(success2)
    
    for(z in 1:dim(mtds)[1]){
      pae[z]=100*trial$exp[mtds[z,1],mtds[z,2]]
    }
    PAE[x1]=sum(pae)
  }else{
    PCS[x1]=0
    PAS[x1]=0
    PCE[x1]=0
    PAE[x1]=0
  }
  
  #Overtoxic rec
  if(dim(over)[1]!=0){
    same3 = matrix(NA,nrow=2000,ncol=dim(over)[1])
    for(q in 1:2000){
      for(z in 1:dim(over)[1]){
        same3[q,z] = (over[z,1]==dose[q,1] & over[z,2]==dose[q,2])
      }
    }
    success3=apply(same3,1,sum)
    OVER[x1]=100*mean(success3)
    
    peover=vector()
    nover<-c()
    for(z in 1:dim(over)[1]){
      peover[z]=100*trial$exp[over[z,1],over[z,2]]
      nover[z]=mean(unlist(lapply(c(1:2000), function(x) trial$n.sim[[x]][over[z,1],over[z,2]])))
      
    }
    PEOVER[x1]=sum(peover)
    NOVER[x1]=sum(nover)
    
  }else{
    OVER[x1]=0
    PEOVER[x1]=0
  }
  
  #Number of patients and DLTs
  tox=pat=vector()
  for(z in 1:2000){
    tox[z]=sum(trial$r.sim[[z]])
    pat[z]=sum(trial$n.sim[[z]])
  }
  ntox[x1]=mean(tox)
  npat[x1]=mean(pat)
  
  SELE=matrix(NA,nrow=3,ncol=3)
  for(i in 1:3){
    for(j in 1:3){
      sel=rep(NA,2000)
      for(k in 1:2000){
        sel[k]=100*(dose[k,1]==i & dose[k,2]==j)
      }
      SELE[i,j]=mean(sel)
    }
  }
  SEL[[x1]]=SELE
  
  EXP[[x1]]=trial$exp * npat[x1]
}
#outputs are vectors or lists of length 15, corresponding to the scenario
p_PCS=PCS
p_PAS=PAS
p_OVER=OVER
p_NOVER=NOVER
p_PEOVER=PEOVER
p_PCE=PCE
p_PAE=PAE
p_ntox=ntox
p_npat=npat
p_EXTRA=EXTRA
p_SEL=SEL
p_EXP=EXP

save.image("PIPE_sims.RData")

