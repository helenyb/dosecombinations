#binomial with prob as 1st argument
binom_prob<-function(prob,N,X){
  dbinom(x=X,size=N,prob=prob)
}

#posterior
cdf_binom1<-function(prob,N,X){
  binom_prob(prob,N,X)/as.numeric(integrate(binom_prob,0,1,N,X)[1])
}

#function to solve for uniroot
cdf_binom3<-function(c,y,N,X){
  as.numeric( integrate(cdf_binom1,0,c,N,X)[1] )-y
}

#probability integral transform
cdf_binom4<-function(rand,N,X){
  as.numeric(uniroot(cdf_binom3,c(0,1),rand,N,X)$root)
  
}

#generate patients
generate_patients<-function(npatients,Nobs,Xobs){
  ex_samps<-c()
  for (i in 1:npatients){
    prob<-cdf_binom4(runif(1),Nobs,Xobs)
    ex_samps[i]<-rbinom(1,1,prob)
  }
  return(ex_samps)
}