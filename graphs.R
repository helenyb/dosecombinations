##plots for results
library(ggplot2)


#loading environments
load("BOIN_sims.RData")
load("KEY_sims.RData")
load("PIPE_sims.RData")
load("BLRM_sims.RData")
load("SF_sims.RData")

#append mean of 1:13
mean_1415<-function(vector){
  return(c(vector,mean(vector[-c(14,15)])))
}


##PCS and PCS in data frame format for ggplot
PCS_data<-data.frame(c(mean_1415(p_PCS2),mean_1415(k_PCS),mean_1415(b_PCS),mean_1415(PCS1),mean_1415(blrm_PCS)),
                     rep(factor(c(1:15,"Mean"),levels=unique(c(1:15,"Mean"))),times=5),
                     rep(c("PIPE","KEYBOARD","BOIN","S-F","BLRM"),each=16))
names(PCS_data)<-c("Result","Scenario","Method")

PAS_data<-data.frame(c(mean_1415(p_PAS2),mean_1415(k_PAS),mean_1415(b_PAS),mean_1415(PAS1),mean_1415(blrm_PAS)),
                     rep(factor(c(1:15,"Mean"),levels=unique(c(1:15,"Mean"))),times=5),
                     rep(c("PIPE","KEYBOARD","BOIN","S-F","BLRM"),each=16))
names(PAS_data)<-c("Result","Scenario","Method")




ggplot()+
  geom_col(data=PCS_data, aes(x=Scenario, y=Result, fill=Method, alpha=factor(1)), width=0.75, position="dodge") +
  
  # geom_col(data=PCS_data, aes(x=Scenario, y=Result, fill=Method)) +
  geom_col(data=PAS_data, aes(x=Scenario, y=Result, fill=Method, alpha=factor(0.9)), width=0.75, position="dodge")+ 
  # geom_bar(stat="identity", position=position_dodge()) +
  scale_alpha_discrete(name="Selection", labels=c("Acceptable","Correct"), range=c(0.17,1)) +
  
  scale_fill_manual(name="Design", values=c("red","black","cyan","green","purple")) +
  ylim(c(0,100))+
  labs(y="% Recommendation", title= "Percentage of Correct & Acceptable Selections (New)") +
  scale_x_discrete(name="Scenario", breaks=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","Mean"), labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","Mean"))

#append mean of 2:14
mean_115<-function(vector){
  return(c(vector,mean(vector[-c(1,15)])))
}



#number of overdose patients in data frame
NOVER_data<-data.frame(c(mean_115(p_NOVER2),mean_115(k_NOVER),mean_115(b_NOVER),mean_115(NOVER1),mean_115(blrm_NOVER)),
                       rep(factor(c(1:15,"Mean"),levels=unique(c(1:15,"Mean"))),times=5),
                       rep(c("PIPE","KEYBOARD","BOIN","S-F","BLRM"),each=16))
names(NOVER_data)<-c("Result","Scenario","Method")



ggplot()+
  geom_col(data=NOVER_data, aes(x=Scenario, y=Result, fill=Method), width=0.75, position="dodge") +
  
  # geom_col(data=PCS_data, aes(x=Scenario, y=Result, fill=Method)) +
  # geom_col(data=PAS_data, aes(x=Scenario, y=Result, fill=Method, alpha=factor(0.9)), width=0.75, position="dodge")+ 
  # geom_bar(stat="identity", position=position_dodge()) +
  # scale_alpha_discrete(name="Selection", labels=c("Acceptable","Correct"), range=c(0.17,1)) +
  
  scale_fill_manual(name="Design", values=c("red","black","cyan","green","purple")) +
  ylim(c(0,36))+
  labs(y="Mean number of Patients", title= "Mean number of Patients Treated on Overly Toxic Doses (New)") +
  scale_x_discrete(name="Scenario", breaks=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","Mean"), labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","Mean"))


#overly toxic recomendations
OVER_data<-data.frame(c(mean_115(p_OVER2),mean_115(k_OVER),mean_115(b_OVER),mean_115(OVER1),mean_115(blrm_OVER)),
                      rep(factor(c(1:15,"Mean"),levels=unique(c(1:15,"Mean"))),times=5),
                      rep(c("PIPE","KEYBOARD","BOIN","S-F","BLRM"),each=16))
names(OVER_data)<-c("Result","Scenario","Method")



ggplot()+
  geom_col(data=OVER_data, aes(x=Scenario, y=Result, fill=Method), width=0.75, position="dodge") +
  
  # geom_col(data=PCS_data, aes(x=Scenario, y=Result, fill=Method)) +
  # geom_col(data=PAS_data, aes(x=Scenario, y=Result, fill=Method, alpha=factor(0.9)), width=0.75, position="dodge")+ 
  # geom_bar(stat="identity", position=position_dodge()) +
  # scale_alpha_discrete(name="Selection", labels=c("Acceptable","Correct"), range=c(0.17,1)) +
  
  scale_fill_manual(name="Design", values=c("red","black","cyan","green","purple")) +
  ylim(c(0,100))+
  labs(y="Percentage of Over-toxic Selections", title= "Mean Percentage of Over-toxic Selections (New)") +
  scale_x_discrete(name="Scenario", breaks=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","Mean"), labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","Mean"))



