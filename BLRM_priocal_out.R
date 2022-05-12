#graphs for BLRM prior calibration

load("BLRM_prior.RData")
ge0_mean=array(NA,dim=c(length(c1),length(c2),length(v1),length(v2),length(v3)))
for(x1 in 1:length(c1)){
  for(x2 in 1:length(c2)){
    for(x3 in 1:length(v1)){
      for(x4 in 1:length(v2)){
        for(x5 in 1:length(v3)){
          ge0_mean[x1,x2,x3,x4,x5]<- exp( mean(log(c(sc1[x1,x2,x3,x4,x5],sc4[x1,x2,x3,x4,x5],sc3[x1,x2,x3,x4,x5],sc4[x1,x2,x3,x4,x5]))))
        } 
      }
    }
  }
}

which(ge0_mean==max(ge0_mean,na.rm=T),arr.ind = T)


load("BLRM_priorcal_ewoc.RData")
BLRM_ewoc1<-sc_ewoc2[c(1:16),1]
BLRM_ewoc2<-sc_ewoc3[c(1:16),1]
BLRM_ewoc3<-sc_ewoc4[c(1:16),1]
BLRM_ewoc4<-100*sc_ewoc5[c(1:16),3]

BLRM_geomean_prior<-(BLRM_ewoc1*BLRM_ewoc2*BLRM_ewoc3*BLRM_ewoc4)^0.25

BLRM_ewoc_frame<-data.frame(rep(seq(0.25,0.1,-.01),times=5),
                            c(BLRM_ewoc1,BLRM_ewoc2,BLRM_ewoc3,BLRM_ewoc4,BLRM_geomean_prior),
                            
                            rep(c("Sc 8","Sc 10","Sc 13","Sc 14","Geo mean"),each=length(seq(0.25,0.1,-.01))))
names(BLRM_ewoc_frame)<-c("epsilon","value","Scenario")
BLRM_ewoc_frame$Scenario<-factor(BLRM_ewoc_frame$Scenario, levels=c("Sc 8","Sc 10","Sc 13","Sc 14","Geo mean"))

library(ggplot2)
ggplot(BLRM_ewoc_frame) + geom_line(aes(x=epsilon, y=value, colour=Scenario)) +
  scale_colour_manual(values=c("red","green","blue","purple","black")) +
  ylim(0,100) +
  labs(x="Epsilon",y="Correct Outcome (%)",color=NULL)


BLRM_ewoc1n<-sc_ewoc2[c(1:16),2]
BLRM_ewoc2n<-sc_ewoc3[c(1:16),2]
BLRM_ewoc3n<-sc_ewoc4[c(1:16),2]
BLRM_ewoc4n<-sc_ewoc5[c(1:16),2]

BLRM_geomean_priorn<-(BLRM_ewoc1n*BLRM_ewoc2n*BLRM_ewoc3n*BLRM_ewoc4n)^0.25

BLRM_ewoc_framen<-data.frame(rep(seq(0.25,0.1,-.01),times=5),
                             c(BLRM_ewoc1n,BLRM_ewoc2n,BLRM_ewoc3n,BLRM_ewoc4n,BLRM_geomean_priorn),
                             
                             rep(c("Sc 8","Sc 10","Sc 13","Sc 14","Geo mean"),each=length(seq(0.25,0.1,-.01))))
names(BLRM_ewoc_framen)<-c("epsilon","value","Scenario")
BLRM_ewoc_framen$Scenario<-factor(BLRM_ewoc_framen$Scenario, levels=c("Sc 8","Sc 10","Sc 13","Sc 14","Geo mean"))


ggplot(BLRM_ewoc_framen) + geom_line(aes(x=epsilon, y=value, colour=Scenario)) +
  scale_colour_manual(values=c("red","green","blue","purple","black")) +
  ylim(0,36.01) +
  labs(x="Epsilon",y="Patients Treated on Over-toxic Dose (mean)",color=NULL) 

