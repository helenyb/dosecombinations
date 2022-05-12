#BOIN prior calibration
library(BOIN)

#############################################################################   

## calibrating interval
set.seed(100)

#scenarios
true1=matrix(nrow=3,ncol=3,c(0.05,0.1,0.15,0.1,0.15,0.2,0.15,0.2,0.3))
true2=matrix(nrow=3,ncol=3,c(0.05,0.1,0.3,0.1,0.2,0.45,0.2,0.3,0.55))
true3=matrix(nrow=3,ncol=3,c(0.15,0.3,0.45,0.3,0.45,0.55,0.45,0.55,0.65))
true4=matrix(nrow=3,ncol=3,c(0.3,0.45,0.5,0.45,0.5,0.55,0.5,0.55,0.6))

a=seq(from=0.85,to=0.40,by=-0.05)
b=seq(from=1.15,to=1.60,by=0.05)

sc1=matrix(NA,nrow=length(a),ncol=length(b))
for(x1 in 1:length(a)){
  for(x2 in 1:length(b)){
    simulation=get.oc.comb(target=0.3, p.true=true1, ncohort=12, cohortsize=3, n.earlystop = NULL, startdose = c(1,1), titration = FALSE,
                           p.saf = 0.3*a[x1], p.tox = 0.3*b[x2], cutoff.eli = 1, extrasafe = FALSE, offset = 0.05,
                           ntrial = 4000, mtd.contour = FALSE, seed = sample(1:1000000,1,replace=F))
    sc1[x1,x2]=simulation$pcs
  }
}
sc1

sc2=matrix(NA,nrow=length(a),ncol=length(b))
for(x1 in 1:length(a)){
  for(x2 in 1:length(b)){
    simulation=get.oc.comb(target=0.3, p.true=true2, ncohort=12, cohortsize=3, n.earlystop = NULL, startdose = c(1,1), titration = FALSE,
                           p.saf = 0.3*a[x1], p.tox = 0.3*b[x2], cutoff.eli = 1, extrasafe = FALSE, offset = 0.05,
                           ntrial = 4000, mtd.contour = FALSE, seed = sample(1:1000000,1,replace=F))
    sc2[x1,x2]=simulation$pcs
  }
}
sc2

sc3=matrix(NA,nrow=length(a),ncol=length(b))
for(x1 in 1:length(a)){
  for(x2 in 1:length(b)){
    simulation=get.oc.comb(target=0.3, p.true=true3, ncohort=12, cohortsize=3, n.earlystop = NULL, startdose = c(1,1), titration = FALSE,
                           p.saf = 0.3*a[x1], p.tox = 0.3*b[x2], cutoff.eli = 1, extrasafe = FALSE, offset = 0.05,
                           ntrial = 4000, mtd.contour = FALSE, seed = sample(1:1000000,1,replace=F))
    sc3[x1,x2]=simulation$pcs
  }
}
sc3

sc4=matrix(NA,nrow=length(a),ncol=length(b))
for(x1 in 1:length(a)){
  for(x2 in 1:length(b)){
    simulation=get.oc.comb(target=0.3, p.true=true4, ncohort=12, cohortsize=3, n.earlystop = NULL, startdose = c(1,1), titration = FALSE,
                           p.saf = 0.3*a[x1], p.tox = 0.3*b[x2], cutoff.eli = 1, extrasafe = FALSE, offset = 0.05,
                           ntrial = 4000, mtd.contour = FALSE, seed = sample(1:1000000,1,replace=F))
    sc4[x1,x2]=simulation$pcs
  }
}
sc4

sc1=matrix(as.numeric(sub("%", "", sc1)),nrow=length(a),ncol=length(b))
sc2=matrix(as.numeric(sub("%", "", sc2)),nrow=length(a),ncol=length(b))
sc3=matrix(as.numeric(sub("%", "", sc3)),nrow=length(a),ncol=length(b))
sc4=matrix(as.numeric(sub("%", "", sc4)),nrow=length(a),ncol=length(b))

sc1
sc2
sc3
sc4

###############################################################################################

## Averaging PCS across scenarios

#geometric mean and maximin
gm=mm=matrix(NA,nrow=length(a),ncol=length(b))
for(i in 1:length(a)){
  for(j in 1:length(b)){
    gm[i,j]=(sc1[i,j]*sc2[i,j]*sc3[i,j]*sc4[i,j])^(0.25)
    mm[i,j]=min(c(sc1[i,j],sc2[i,j],sc3[i,j],sc4[i,j]))
  }
}
round(gm,1)
round(mm,1)
max(gm)
min(gm)
which(gm==max(gm),arr.ind=TRUE)

x = seq(0.85,0.40,by=-0.05)
y = seq(1.15,1.60,by=0.05)
data <- expand.grid(X=x, Y=y)
data$GM <- as.vector(gm)

library(ggplot2)
pdf(file="fig_BOIN_int.pdf",width=4,height=3)
ggplot(data,aes(x=X,y=Y,fill=GM)) + geom_tile() + scale_fill_gradient(low="white", high="darkblue",limits=c(53,63)) + 
  labs(x=expression(a[1]),y=expression(a[2]),fill="Mean \nPCS (%)") + theme(legend.key.height=unit(1, "cm")) +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))
dev.off()

ggplot(data,aes(x=data$X,y=data$Y,fill=(data$MM))) + geom_tile() + scale_fill_gradient(low="white", high="black")

#################################################################################################

## calibrate safety constraint
set.seed(400)


s8=matrix(c(0.05,0.10,0.30,0.10,0.20,0.45,0.20,0.30,0.55),nrow=3,ncol=3)
s10=matrix(c(0.15,0.30,0.45,0.30,0.45,0.55,0.45,0.55,0.65),nrow=3,ncol=3)
s13=matrix(c(0.30,0.45,0.50,0.45,0.50,0.55,0.50,0.55,0.60),nrow=3,ncol=3)
s14=matrix(c(0.45,0.50,0.55,0.50,0.55,0.60,0.55,0.60,0.65),nrow=3,ncol=3)



epsilon=seq(from=1,to=0.4,by=-0.05)
a=0.65
b=1.4

BOINtox_sc14=vector()
BOINtox_sc14b=vector()
for(y1 in 1:length(epsilon)){
  simulation = get.oc.comb(target=0.3, p.true=s14, ncohort=12, cohortsize=3, n.earlystop = NULL, startdose = c(1,1), titration = FALSE,
                           p.saf = 0.3*a, p.tox = 0.3*b, cutoff.eli = epsilon[y1], extrasafe = FALSE, offset = 0.05,
                           ntrial = 4000, mtd.contour = FALSE, seed = sample(1:1000000,1,replace=F))
  pstop = 100-sum(simulation$selpercent)
  
  BOINtox_sc14[y1] = pstop
  BOINtox_sc14b[y1] =sum(simulation$npatients*as.numeric(simulation$p.true>0.3))
}
#tox_sc1

BOINtox_sc8=vector()
BOINtox_sc8a=vector()
BOINtox_sc8b=vector()
for(y1 in 1:length(epsilon)){
  simulation = get.oc.comb(target=0.3, p.true=s8, ncohort=12, cohortsize=3, n.earlystop = NULL, startdose = c(1,1), titration = FALSE,
                           p.saf = 0.3*a, p.tox = 0.3*b, cutoff.eli = epsilon[y1], extrasafe = FALSE, offset = 0.05,
                           ntrial = 4000, mtd.contour = FALSE, seed = sample(1:1000000,1,replace=F))
  pnostop = sum(simulation$selpercent)
  BOINtox_sc8b[y1]<- sum(simulation$npatients*as.numeric(simulation$p.true>0.3))
  pcs = simulation$pcs
  BOINtox_sc8[y1] = simulation$pcs
  BOINtox_sc8a[y1] = as.numeric(sub("%", "", pcs))
}



BOINtox_sc10=vector()
BOINtox_sc10a=vector()
BOINtox_sc10b=vector()
for(y1 in 1:length(epsilon)){
  simulation = get.oc.comb(target=0.3, p.true=s10, ncohort=12, cohortsize=3, n.earlystop = NULL, startdose = c(1,1), titration = FALSE,
                           p.saf = 0.3*a, p.tox = 0.3*b, cutoff.eli = epsilon[y1], extrasafe = FALSE, offset = 0.05,
                           ntrial = 4000, mtd.contour = FALSE, seed = sample(1:1000000,1,replace=F))
  pnostop = sum(simulation$selpercent)
  pcs = simulation$pcs
  BOINtox_sc10b[y1]<- sum(simulation$npatients*as.numeric(simulation$p.true>0.3))
  BOINtox_sc10[y1] = simulation$pcs
  BOINtox_sc10a[y1] = as.numeric(sub("%", "", pcs))
}


BOINtox_sc13=vector()
BOINtox_sc13a=vector()
BOINtox_sc13b=vector()
for(y1 in 1:length(epsilon)){
  simulation = get.oc.comb(target=0.3, p.true=s13, ncohort=12, cohortsize=3, n.earlystop = NULL, startdose = c(1,1), titration = FALSE,
                           p.saf = 0.3*a, p.tox = 0.3*b, cutoff.eli = epsilon[y1], extrasafe = FALSE, offset = 0.05,
                           ntrial = 4000, mtd.contour = FALSE, seed = sample(1:1000000,1,replace=F))
  pnostop = sum(simulation$selpercent)
  BOINtox_sc13b[y1]<- sum(simulation$npatients*as.numeric(simulation$p.true>0.3))
  pcs = simulation$pcs
  BOINtox_sc13[y1] = simulation$pcs
  BOINtox_sc13a[y1] = as.numeric(sub("%", "", pcs))
}



boin_plot8a<-data.frame(epsilon,BOINtox_sc8b,"Sc 8")
names(boin_plot8a)<-c("epsilon","value","Scenario")
boin_plot10a<-data.frame(epsilon,BOINtox_sc10b,"Sc 10")
names(boin_plot10a)<-c("epsilon","value","Scenario")
boin_plot13a<-data.frame(epsilon,BOINtox_sc13b,"Sc 13")
names(boin_plot13a)<-c("epsilon","value","Scenario")
boin_plot14a<-data.frame(epsilon,BOINtox_sc14b,"Sc 14")
names(boin_plot14a)<-c("epsilon","value","Scenario")
boin_plotGa<-data.frame(epsilon,(BOINtox_sc8b*BOINtox_sc10b*BOINtox_sc13b*BOINtox_sc14b)^0.25,"Geo mean")
names(boin_plotGa)<-c("epsilon","value","Scenario")

boin_plota<-rbind(boin_plot8a,boin_plot10a,boin_plot13a,boin_plot14a,boin_plotGa)
names(boin_plota)<-c("epsilon","value","Scenario")

ggplot(boin_plota) + geom_line(aes(x=epsilon, y=value, colour=Scenario)) +
  scale_colour_manual(values=c("red","green","blue","purple","black")) +
  ylim(0,36.01) +
  labs(x="Epsilon",y="Patients Treated on Over-toxic Dose (mean)",color=NULL) 




boin_plot8<-data.frame(epsilon,BOINtox_sc8a,"Sc 8")
names(boin_plot8)<-c("epsilon","value","Scenario")
boin_plot10<-data.frame(epsilon,BOINtox_sc10a,"Sc 10")
names(boin_plot10)<-c("epsilon","value","Scenario")
boin_plot13<-data.frame(epsilon,BOINtox_sc13a,"Sc 13")
names(boin_plot13)<-c("epsilon","value","Scenario")
boin_plot14<-data.frame(epsilon,BOINtox_sc14,"Sc 14")
names(boin_plot14)<-c("epsilon","value","Scenario")
boin_plotG<-data.frame(epsilon,(BOINtox_sc8a*BOINtox_sc10a*BOINtox_sc13a*BOINtox_sc14)^0.25,"Geo mean")
names(boin_plotG)<-c("epsilon","value","Scenario")
boin_plot<-rbind(boin_plot8,boin_plot10,boin_plot13,boin_plot14,boin_plotG)
names(boin_plot)<-c("epsilon","value","Scenario")

ggplot(boin_plot) + geom_line(aes(x=epsilon, y=value, colour=Scenario)) +
  scale_colour_manual(values=c("red","green","blue","purple","black")) +
  ylim(0,100) +
  labs(x="Epsilon",y="Correct Outcome (%)",color=NULL)






library(ggplot2)
pdf(file="fig_BOIN_safe.pdf",width=4.6,height=3.4)
ggplot() + geom_line(aes(x=epsilon,y=tox_sc1,col="red"),size=0.8) + geom_line(aes(x=epsilon,y=tox_sc2a,col="green"),size=0.8) + ylim(0,100) +
  labs(x=expression(epsilon),y="Correct Outcome (%)",color=NULL) +  scale_color_manual(values=c("green"="green","red"="red"), labels=c("Sc 13","Sc 14")) +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))
dev.off()
