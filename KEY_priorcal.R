#KEY PRIOR CALIBRATION

library(Keyboard)
#source in the altered file to conduct KEYBOARD design as described by authors
source("KEY_SIM.R")
#scenarios
s8=matrix(c(0.05,0.10,0.30,0.10,0.20,0.45,0.20,0.30,0.55),nrow=3,ncol=3)
s10=matrix(c(0.15,0.30,0.45,0.30,0.45,0.55,0.45,0.55,0.65),nrow=3,ncol=3)
s13=matrix(c(0.30,0.45,0.50,0.45,0.50,0.55,0.50,0.55,0.60),nrow=3,ncol=3)
s14=matrix(c(0.45,0.50,0.55,0.50,0.55,0.60,0.55,0.60,0.65),nrow=3,ncol=3)


#values
epsilon=seq(from=1,to=0.4,by=-0.05)
a=0.65
b=1.4

KEYtox_sc14=vector()
KEYtox_sc14b=vector()
for(y1 in 1:length(epsilon)){
  simulation = 
  keySIM(
    target=0.3,
    p.true=s14,
    ncohort=12,
    cohortsize=3,
    n.earlystop = 1000,
    marginL = 0.09,
    marginR = 0.09,
    startdose = c(1, 1),
    cutoff.eli = epsilon[y1],
    extrasafe = FALSE,
    offset = 0.05,
    ntrial = 4000
  )
  
  
  pstop = 100-sum(simulation$selpercent)
  
  KEYtox_sc14[y1] = pstop
  KEYtox_sc14b[y1] =sum(simulation$nptsdose*as.numeric(simulation$p.true>0.3))
}
#tox_sc1

KEYtox_sc8=vector()
KEYtox_sc8a=vector()
KEYtox_sc8b=vector()
for(y1 in 1:length(epsilon)){
  simulation = 
    keySIM(
      target=0.3,
      p.true=s8,
      ncohort=12,
      cohortsize=3,
      n.earlystop = 1000,
      marginL = 0.09,
      marginR = 0.09,
      startdose = c(1, 1),
      cutoff.eli = epsilon[y1],
      extrasafe = FALSE,
      offset = 0.05,
      ntrial = 4000
    )
  pnostop = sum(simulation$selpercent)
  KEYtox_sc8b[y1]<- sum(simulation$nptsdose*as.numeric(simulation$p.true>0.3))
  pcs = simulation$pcs
  KEYtox_sc8[y1] = simulation$pcs
  KEYtox_sc8a[y1] = as.numeric(sub("%", "", pcs))
}



KEYtox_sc10=vector()
KEYtox_sc10a=vector()
KEYtox_sc10b=vector()
for(y1 in 1:length(epsilon)){
  simulation = 
    keySIM(
      target=0.3,
      p.true=s10,
      ncohort=12,
      cohortsize=3,
      n.earlystop = 1000,
      marginL = 0.09,
      marginR = 0.09,
      startdose = c(1, 1),
      cutoff.eli = epsilon[y1],
      extrasafe = FALSE,
      offset = 0.05,
      ntrial = 4000
    )
  pnostop = sum(simulation$selpercent)
  pcs = simulation$pcs
  KEYtox_sc10b[y1]<- sum(simulation$nptsdose*as.numeric(simulation$p.true>0.3))
  KEYtox_sc10[y1] = simulation$pcs
  KEYtox_sc10a[y1] = as.numeric(sub("%", "", pcs))
}


KEYtox_sc13=vector()
KEYtox_sc13a=vector()
KEYtox_sc13b=vector()
for(y1 in 1:length(epsilon)){
  simulation = 
    keySIM(
      target=0.3,
      p.true=s13,
      ncohort=12,
      cohortsize=3,
      n.earlystop = 1000,
      marginL = 0.09,
      marginR = 0.09,
      startdose = c(1, 1),
      cutoff.eli = epsilon[y1],
      extrasafe = FALSE,
      offset = 0.05,
      ntrial = 4000
    )
  pnostop = sum(simulation$selpercent)
  KEYtox_sc13b[y1]<- sum(simulation$nptsdose*as.numeric(simulation$p.true>0.3))
  pcs = simulation$pcs
  KEYtox_sc13[y1] = simulation$pcs
  KEYtox_sc13a[y1] = as.numeric(sub("%", "", pcs))
}

##plots

key_plot8a<-data.frame(epsilon,KEYtox_sc8b,"Sc 8")
names(key_plot8a)<-c("epsilon","value","Scenario")
key_plot10a<-data.frame(epsilon,KEYtox_sc10b,"Sc 10")
names(key_plot10a)<-c("epsilon","value","Scenario")
key_plot13a<-data.frame(epsilon,KEYtox_sc13b,"Sc 13")
names(key_plot13a)<-c("epsilon","value","Scenario")
key_plot14a<-data.frame(epsilon,KEYtox_sc14b,"Sc 14")
names(key_plot14a)<-c("epsilon","value","Scenario")
key_plotGa<-data.frame(epsilon,(KEYtox_sc8b*KEYtox_sc10b*KEYtox_sc13b*KEYtox_sc14b)^0.25,"Geo mean")
names(key_plotGa)<-c("epsilon","value","Scenario")

key_plota<-rbind(key_plot8a,key_plot10a,key_plot13a,key_plot14a,key_plotGa)
names(key_plota)<-c("epsilon","value","Scenario")

ggplot(key_plota) + geom_line(aes(x=epsilon, y=value, colour=Scenario)) +
  scale_colour_manual(values=c("red","green","blue","purple","black")) +
  ylim(0,36.01) +
  labs(x="Epsilon",y="Patients Treated on Over-toxic Dose (mean)",color=NULL) 




key_plot8<-data.frame(epsilon2,KEYtox_sc8a,"Sc 8")
names(key_plot8)<-c("epsilon","value","Scenario")
key_plot10<-data.frame(epsilon2,KEYtox_sc10a,"Sc 10")
names(key_plot10)<-c("epsilon","value","Scenario")
key_plot13<-data.frame(epsilon2,KEYtox_sc13a,"Sc 13")
names(key_plot13)<-c("epsilon","value","Scenario")
key_plot14<-data.frame(epsilon2,KEYtox_sc14,"Sc 14")
names(key_plot14)<-c("epsilon","value","Scenario")
key_plotG<-data.frame(epsilon2,(KEYtox_sc8a*KEYtox_sc10a*KEYtox_sc13a*KEYtox_sc14)^0.25,"Geo mean")
names(key_plotG)<-c("epsilon","value","Scenario")
key_plot<-rbind(key_plot8,key_plot10,key_plot13,key_plot14,key_plotG)
names(key_plot)<-c("epsilon","value","Scenario")

ggplot(key_plot) + geom_line(aes(x=epsilon, y=value, colour=Scenario)) +
  scale_colour_manual(values=c("red","green","blue","purple","black")) +
  ylim(0,100) +
  labs(x="Epsilon",y="Correct Outcome (%)",color=NULL)

