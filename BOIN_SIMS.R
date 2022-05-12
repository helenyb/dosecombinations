#BOIN SIMS

#Scenarios
library(BOIN)
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


#initialize empty lists
PCS=PAS=PCE=PAE=OVER=PEOVER=NOVER=ntox=npat=rep(NA,15)
EXTRA=SEL=EXP=list()
#conduct simulations for all 15 scenarios
for(x1 in 1:15){
  trial=get.oc.comb(target=0.3, p.true=scenario[[x1]], ncohort=12, cohortsize=3, n.earlystop = NULL, startdose = c(1,1), titration = FALSE,
                    p.saf = 0.3*0.65, p.tox = 0.3*1.4, cutoff.eli = 0.84, extrasafe = FALSE, offset = 0.05,
                    ntrial = 2000, mtd.contour = FALSE, seed = 6)
  mtd=which(scenario[[x1]]==0.30,arr.ind=T)
  mtds=which((scenario[[x1]]>0.16 & scenario[[x1]]<0.33),arr.ind=T)
  over=which(scenario[[x1]]>0.33,arr.ind=T)
  
  #PCS, PAS, PCE
  if(dim(mtd)[1]!=0){
    pcs=vector()
    pas=vector()
    pce=vector()
    pae=vector()
    for(z in 1:dim(mtd)[1]){
      pcs[z]=trial$selpercent[mtd[z,1],mtd[z,2]]
      pce[z]=100 * trial$npatients[mtd[z,1],mtd[z,2]]/trial$totaln
    }
    for(z in 1:dim(mtds)[1]){
      pas[z]=trial$selpercent[mtds[z,1],mtds[z,2]]
      pae[z]=100 * trial$npatients[mtds[z,1],mtds[z,2]]/trial$totaln
    }
    PCS[x1]=sum(pcs)
    PAS[x1]=sum(pas)
    PCE[x1]=sum(pce)
    PAE[x1]=sum(pae)
    EXTRA[[x1]]=pcs
  }else{
    PCS[x1]=0
    PAS[x1]=0
    PCE[x1]=0
    PAE[x1]=0
  }
  
  #Overtoxic rec
  if(dim(over)[1]!=0){
    otox=vector()
    peover=vector()
    nover<-c()
    for(z in 1:dim(over)[1]){
      otox[z]=trial$selpercent[over[z,1],over[z,2]]
      peover[z]=100*trial$npatients[over[z,1],over[z,2]]/trial$totaln
      nover[z]=trial$npatients[over[z,1],over[z,2]]
    }
    OVER[x1]=sum(otox)
    NOVER[x1]=sum(nover)
    PEOVER[x1]=sum(peover)
  }else{
    OVER[x1]=0
    PEOVER[x1]=0
  }
  
  #Number of patients and DLTs
  ntox[x1]=trial$totaltox
  npat[x1]=trial$totaln
  
  #sel and exp
  SEL[[x1]]=trial$selpercent
  EXP[[x1]]=trial$npatients
}
#outputs are vectors or lists of length 15, corresponding to the scenario
b_PCS=PCS
b_PAS=PAS
b_OVER=OVER
b_NOVER=NOVER
b_PEOVER=PEOVER
b_PEOVER[14]=100
b_PCE=PCE
b_PAE=PAE
b_ntox=ntox
b_npat=npat
b_EXTRA=EXTRA
b_SEL=SEL
b_EXP=EXP

save.image("BOIN_sims.RData")