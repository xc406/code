
####
## clip read counts to avoid extreme values
#####
file <- '~/Documents/google-python-exercises/gold_standard/MAX100_DGF_K562_htseq_output'
count <- read.delim(file,header=F)
count <- count[-(nrow(count) - 0:4),]
c <- count[,2]
c<-c[which(c>20)]

#file <- '~/Documents/google-python-exercises/aupr_curves/Irf4100_FAIRE_htseq_output'
#faire <- read.delim(file,header=F)
#faire <- faire[-(nrow(faire) - 0:4),]
#f <- faire[,2]
#f<-f[which(f>20)]

##initialize parameters
initNegBinParams <- function(R,Ez,DampNegBin=0.0){
  Ez <- Ez*(1-DampNegBin)+0.5*(DampNegBin);
  myM<-sum(R*Ez)/sum(Ez)
  myV<-sum(R^2 * Ez)/sum(Ez) - myM^2
  B1<-myM/myV;
      if(B1 > 0.999999){
      B1 <- 0.999999   
      }      
  A1<-myM*B1/(1-B1);
  list(logA=log(A1),logitB=qlogis(B1))
}

rc <- range(c)

q90 <- quantile(c,0.9)
c90 <- c[which(c>q90)]
c <- count[,2]
library(MASS)

cf90 <- fitdistr(c90, "Negative Binomial") 

alpha1 <- cf90$estimate[1]
mean1 <- cf90$estimate[2]

#q10 <- quantile(c,0.1)
c10 <- c[which(c<q90)]

cf10 <- fitdistr(c10, "Negative Binomial")

alpha0 <- cf10$estimate[1]
mean0 <- cf10$estimate[2]

d1<-dnbinom(c,size=alpha1,mu=mean1)
d0<-dnbinom(c,size=alpha0,mu=mean0)
plot(x=c,y=d1,log='xy',col='blue',xlim=c(1,10000),xlab='read count',ylab='probabilities')
points(x=c,y=d0,col='red')
legend('bottomleft',col=c('blue','red'),legend=c('90 percentile read count-open chromatin','90 percentile read count-closed chromatin'),bty='o',pch=c(0,0))
#pval <- pnbinom(c, size=3.79, mu=578 , lower.tail = FALSE, log.p = FALSE)

## E-step
estep <- function(c,pi,alpha1,mean1,alpha0,mean0){
	pi_estep <- (pi*dnbinom(c,size=alpha1,mu=mean1)) / ( pi*dnbinom(c,size=alpha1,mu=mean1) + (1-pi)*dnbinom(c,size=alpha0,mu=mean0) )
    return(pi_estep)
}

tau1func <- function(alpha1_temp,e.step,obs){
	pi_temp <- mean(e.step)
	tau1_temp <- alpha1_temp*sum(e.step[which(e.step>pi_temp)])/(alpha1_temp*sum(e.step[which(e.step>pi_temp)])+sum(obs[which(e.step>pi_temp)]*e.step[which(e.step>pi_temp)]))
	return(tau1_temp)
}

tau0func <- function(alpha0_temp,e.step,obs){
	pi_temp <- mean(e.step)
	tau0_temp <- alpha0_temp*(max(obs)-sum(e.step[which(e.step>pi_temp)]))/(sum(obs[which(e.step>pi_temp)])+alpha0_temp*(max(obs)-sum(e.step[which(e.step>pi_temp)]))-sum(obs[which(e.step>pi_temp)]*e.step[which(e.step>pi_temp)]))
	return(tau0_temp)
}

qfunc <- function(alpha1,alpha0,e.step,tau1,tau0){
	sum(e.step%*%(log10(0.5)+dnbinom(c,size=alpha1,mu=alpha1/tau1-alpha1,log=TRUE)))+sum((1-e.step)%*%(log10(0.5)+dnbinom(c,size=alpha0,mu=alpha0/tau0-alpha0,log=TRUE)))
}

## Maximization Step
mstep <- function(obs,e.step){
  
  # estimate pi
  pi_temp <- mean(e.step)
  
  # estimate alpha,mu
  
  alpha1_temp <- optim(c(alpha1,alpha0,e.step,tau1func(alpha1,e.step,obs),tau0func(alpha0,e.step,obs)),qfunc,gr=NULL,method="BFGS")
  alpha0_temp <- sum(obs*e.step) / sum(e.step)
  tau0_temp <- alpha0_temp*(max(obs)-sum(e.step[which(e.step>pi_temp)]))
		/(sum(obs[which(e.step>pi_temp)])+alpha0_temp*(max(obs)-sum(e.step[which(e.step>pi_temp)]))-sum(obs[which(e.step>pi_temp)]*e.step[which(e.step>pi_temp)]))
  tau1_temp <- alpha1_temp*sum(e.step[which(e.step>pi_temp)])/(alpha1_temp*sum(e.step[which(e.step>pi_temp)])+sum(obs[which(e.step>pi_temp)]*e.step[which(e.step>pi_temp)]))
  
  mean1_temp <- (alpha1_temp/tau1_temp)-alpha1_temp  
  mean0_temp <- (alpha0_temp/tau0_temp)-alpha0_temp
  
  list(pi_temp,alpha1_temp,mean1_temp,alpha0_temp,mean0_temp)   
}


## EM Algorithm
em.algo <- function(obs,pi_inits,alpha1_inits,alpha0_inits,maxit=1000,tol=1e-6){
  # Initial parameter estimates
  flag <- 0
  pi_cur <- pi_inits; p_cur <- p_inits; q_cur <- q_inits
   
  # Iterate between expectation and maximization steps
  for(i in 1:maxit){
    cur <- c(pi_cur,p_cur,q_cur)
    new <- mstep(obs,estep(obs, pi_cur, p_cur, q_cur))
    pi_new <- new[[1]]; p_new <- new[[2]]; q_new <- new[[3]]
    new_step <- c(pi_new,p_new,q_new)
    
    # Stop iteration if the difference between the current and new estimates is less than a tolerance level
    if( all(abs(cur - new_step) < tol) ){ flag <- 1; break}


    # Otherwise continue iteration
    pi_cur <- pi_new; p_cur <- p_new; q_cur <- q_new
  }
  if(!flag) warning("Didn’t converge\n")
  
  list(pi_cur, p_cur, q_cur)
}



## Calculate Information matrix
Info.Mat.function <- function(obs, pi.est, p.est, q.est){
  estep.est <- estep(obs,pi.est, p.est, q.est)
  info.mat <- matrix(rep(0,9),3,3)
  info.mat[1,1] <- - sum(estep.est)/(pi.est^2)  - sum((1-estep.est))/((1-pi.est)^2) 
  info.mat[2,2] <- - sum(estep.est*obs)/(p.est^2) - sum(estep.est*(1-obs))/((1-p.est)^2)
  info.mat[3,3] <- - sum((1-estep.est)*obs)/(q.est^2) - sum((1-estep.est)*(1-obs))/((1-q.est)^2)
  return(-info.mat)
}


## Generate sample data
n <- 5000 
pi_true <- 0.90 # prob of using first coin
p_true <-  0.60 # the first coin has P(heads) = 0.60
q_true <-  0.50 # the second coin has P(heads) = 0.50
true <- c(pi_true,p_true,q_true)
u <- ifelse(runif(n)<pi_true, rbinom(n,1,p_true),rbinom(n,1,q_true))


## Set parameter estimates
pi_init = 0.70; p_init = 0.70; q_init = 0.60


## Run EM Algorithm
output <- em.algo(u,pi_init,p_init,q_init)


## Calculate Confidence Intervals
sd.out <- sqrt(diag(solve(Info.Mat.function(u,output[[1]],output[[2]],output[[3]]))))
data.frame("Truth" = true, "EM Estimate" = unlist(output), "Lower CI" = unlist(output) - qnorm(.975)*sd.out, "Upper CI" = unlist(output) + qnorm(.975)*sd.out)