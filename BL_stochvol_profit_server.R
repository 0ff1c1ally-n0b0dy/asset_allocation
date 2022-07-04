#Profit sensitivity to confidence levels
setwd("/home-b/mandrei/Desktop/Research")

err <- file("FB error stochvol parallel 10^-7 10^-8 10^-7 10^-8.txt", open="wt")
sink(err, type="message")

writeLines(c(""), "test.txt")


################################
###############################
#BL time changing volatility 
###############################
################################

library(data.table)
library(psych)
library(GPArotation)
library(stochvol)
library(MCMCpack)
library(expm)
library(mvtnorm)
library(armspp)
library(invgamma)
library(RGeode)
library(Rcpp)
library(coda)
library(quadprog)
#set.seed(12345)

########################
#Functions 
########################
#Need to create a function for irrational powers of symmetric matrix. 
#This is because we need Q_t^(d/2), where d has a posterior dist that is cont. 
#The %^% was tested and it does not work for the interval (-1,1) and 
#matrix.power from matrixcalc doesn't work for non-integers
power.matrix<-function(M,power)
{
  diagonalize<-eigen(M)
  if(isSymmetric(M)==TRUE)
  {
    if(sum(diagonalize$values<0)==0)
    {
      res<-diagonalize$vectors%*%diag(diagonalize$values^power)%*%t(diagonalize$vectors)
      return(res)
    }
    else
    {
      return(print("Matrix is not positive definite"))
    }
  }
  else
  {
    #return(print("Matrix is not symmetric"))
    if(sum(diagonalize$values<0)==0)
    {
      res<-diagonalize$vectors%*%diag(diagonalize$values^power)%*%solve(diagonalize$vectors)
      return(res) 
    }
    else
    {
      return(print("Matrix is not positive definite"))
    }
  }
}

#Those functions coded like thie give us very close to 0 for numerator and very large for denominator
#(we go above the machine memory allocation for variables maybe). Also, I do not like having 
# workarounds such as bigq(), which is working with fractional approximations instead of the numbers
#One idea is to code the whole ratio as a function by grouping the terms similar to each other also.

ratio.mh.t<-function(k,d,Qt.cur,Qt.cand,A.inv,Qtp1.cur.inv,epsilont)
{
  Sigma.epsilont.cur<-diag(diag(Qt.cur)^(-1/2))%*%Qt.cur%*%diag(diag(Qt.cur)^(-1/2))
  Sigma.epsilont.cand<-diag(diag(Qt.cand)^(-1/2))%*%Qt.cand%*%diag(diag(Qt.cand)^(-1/2))
  det1<-(det(Sigma.epsilont.cand)/det(Sigma.epsilont.cur))^(-0.5)
  det2<-(det(Qt.cand)/det(Qt.cur))^((k*d+1)/2)
  expo.ratio1<-exp(-0.5*k*sum(diag((power.matrix(Qt.cand,d/2)%*%A.inv%*%power.matrix(Qt.cand,d/2)-
                                      power.matrix(Qt.cur,d/2)%*%A.inv%*%power.matrix(Qt.cur,d/2))%*%Qtp1.cur.inv)))
  expo.ratio2<-exp(-0.5*sum(diag((solve(Sigma.epsilont.cand)-solve(Sigma.epsilont.cur)-
                                    solve(Qt.cand)+solve(Qt.cur))%*%epsilont%*%t(epsilont))))
  res<-det1*det2*expo.ratio1*expo.ratio2
  return(res)
}

ratio.mh.T<-function(QT.cur,QT.cand,epsilonT)
{
  Sigma.epsilonT.cur<-diag(diag(QT.cur)^(-1/2))%*%QT.cur%*%diag(diag(QT.cur)^(-1/2))
  Sigma.epsilonT.cand<-diag(diag(QT.cand)^(-1/2))%*%QT.cand%*%diag(diag(QT.cand)^(-1/2))
  res<-exp(-0.5*(sum(diag((solve(Sigma.epsilonT.cand)-solve(Sigma.epsilonT.cur)-
                             solve(QT.cand)+solve(QT.cur))%*%epsilonT%*%t(epsilonT)))))
  return(res)
}

#optimized for speed versions of the above 2 functions; we do not need to compute again Sigma.epsilon.cur.inv and Q.cur.inv
#since we have them already computed in the big Gibbs Sampler loop

ratio.mh.t.optim.solve<-function(k,d,Qt.cur, Qt.cur.inv, Qt.cand,A.inv,Qtp1.cur.inv,epsilont, Sigma.epsilont.cur.inv)
{
  Sigma.epsilont.cand<-diag(diag(Qt.cand)^(-1/2))%*%Qt.cand%*%diag(diag(Qt.cand)^(-1/2))
  det1<-(det(Sigma.epsilont.cand)^(-0.5))/(det(Sigma.epsilont.cur.inv)^(0.5))
  det2<-(det(Qt.cand)/det(Qt.cur))^((k*d+1)/2)
  expo.ratio1<-exp(-0.5*k*sum(diag((power.matrix(Qt.cand,d/2)%*%A.inv%*%power.matrix(Qt.cand,d/2)-
                                      power.matrix(Qt.cur,d/2)%*%A.inv%*%power.matrix(Qt.cur,d/2))%*%Qtp1.cur.inv)))
  expo.ratio2<-exp(-0.5*sum(diag((solve(Sigma.epsilont.cand)-Sigma.epsilont.cur.inv-
                                    solve(Qt.cand)+Qt.cur.inv)%*%epsilont%*%t(epsilont))))
  res<-det1*det2*expo.ratio1*expo.ratio2
  return(res)
}

ratio.mh.T.optim.solve<-function(QT.cur.inv,QT.cand,epsilonT, Sigma.epsilonT.cur.inv)
{
  Sigma.epsilonT.cand<-diag(diag(QT.cand)^(-1/2))%*%QT.cand%*%diag(diag(QT.cand)^(-1/2))
  res<-exp(-0.5*(sum(diag((solve(Sigma.epsilonT.cand)-Sigma.epsilonT.cur.inv-
                             solve(QT.cand)+QT.cur.inv)%*%epsilonT%*%t(epsilonT)))))
  return(res)
}
#Use chol2inv(chol()) instead of solve(), since it is faster especially for large matrices
ratio.mh.t.optim<-function(k,d,Qt.cur, Qt.cur.inv, Qt.cand,A.inv,Qtp1.cur.inv,epsilont, Sigma.epsilont.cur.inv)
{
  Sigma.epsilont.cand<-diag(diag(Qt.cand)^(-1/2))%*%Qt.cand%*%diag(diag(Qt.cand)^(-1/2))
  det1<-(det(Sigma.epsilont.cand)^(-0.5))/(det(Sigma.epsilont.cur.inv)^(0.5))
  det2<-(det(Qt.cand)/det(Qt.cur))^((k*d+1)/2)
  expo.ratio1<-exp(-0.5*k*sum(diag((power.matrix(Qt.cand,d/2)%*%A.inv%*%power.matrix(Qt.cand,d/2)-
                                      power.matrix(Qt.cur,d/2)%*%A.inv%*%power.matrix(Qt.cur,d/2))%*%Qtp1.cur.inv)))
  expo.ratio2<-exp(-0.5*sum(diag((chol2inv(chol(Sigma.epsilont.cand))-Sigma.epsilont.cur.inv-
                                    chol2inv(chol(Qt.cand))+Qt.cur.inv)%*%epsilont%*%t(epsilont))))
  res<-det1*det2*expo.ratio1*expo.ratio2
  return(res)
}

ratio.mh.T.optim<-function(QT.cur.inv,QT.cand,epsilonT, Sigma.epsilonT.cur.inv)
{
  Sigma.epsilonT.cand<-diag(diag(QT.cand)^(-1/2))%*%QT.cand%*%diag(diag(QT.cand)^(-1/2))
  res<-exp(-0.5*(sum(diag((chol2inv(chol(Sigma.epsilonT.cand))-Sigma.epsilonT.cur.inv-
                             chol2inv(chol(QT.cand))+QT.cur.inv)%*%epsilonT%*%t(epsilonT)))))
  return(res)
}

#Take logs and use -2*log(u) is chi square 2 degrees of freedom (or just sample a uniform and plug it into -2*log(u))
ratio.mh.t.optim.log<-function(k,d,Qt.cur, Qt.cur.inv, Qt.cand,A.inv,Qtp1.cur.inv,epsilont, Sigma.epsilont.cur.inv)
{
  Sigma.epsilont.cand<-diag(diag(Qt.cand)^(-1/2))%*%Qt.cand%*%diag(diag(Qt.cand)^(-1/2))
  log.det1<-(-0.5)*log(det(Sigma.epsilont.cand))-(0.5)*log(det(Sigma.epsilont.cur.inv))
  log.det2<-((k*d+1)/2)*(log(det(Qt.cand))-log(det(Qt.cur)))
  log.expo.ratio1<-(-0.5)*k*sum(diag((power.matrix(Qt.cand,d/2)%*%A.inv%*%power.matrix(Qt.cand,d/2)-
                                        power.matrix(Qt.cur,d/2)%*%A.inv%*%power.matrix(Qt.cur,d/2))%*%Qtp1.cur.inv))
  log.expo.ratio2<-(-0.5)*sum(diag((chol2inv(chol(Sigma.epsilont.cand))-Sigma.epsilont.cur.inv-
                                      chol2inv(chol(Qt.cand))+Qt.cur.inv)%*%epsilont%*%t(epsilont)))
  log.res<-log.det1+log.det2+log.expo.ratio1+log.expo.ratio2
  return(log.res)
}

ratio.mh.T.optim.log<-function(QT.cur.inv,QT.cand,epsilonT, Sigma.epsilonT.cur.inv)
{
  Sigma.epsilonT.cand<-diag(diag(QT.cand)^(-1/2))%*%QT.cand%*%diag(diag(QT.cand)^(-1/2))
  log.res<-(-0.5)*(sum(diag((chol2inv(chol(Sigma.epsilonT.cand))-Sigma.epsilonT.cur.inv-
                               chol2inv(chol(QT.cand))+QT.cur.inv)%*%epsilonT%*%t(epsilonT))))
  return(log.res)
}


########################################
#Import data and set working directory
########################################
#setwd("D:/School USB/WPI/1st Year Postdoc/Research/Professor Zou/BL time changing volatility")

#r<-read.table("D:/School USB/WPI/1st Year Postdoc/Research/Professor Zou/BL time changing volatility
#             /data/S&P 500 Returns taking out.csv",sep=",", header = TRUE)

r<-read.table("S&P 500 Returns taking out.csv",sep=",", header = TRUE)

r<-r[,-1]

#####################################
#Initial values for the Gibbs-Sampler
#####################################
burn<-10

tickers<-colnames(r)
p<-length(tickers)
T.max<-nrow(r)

prices.test<-read.table("Prices Jan 2018.csv",header=T,fill=T,sep=",")
ld.test<-nrow(prices.test)
# returns.test<-matrix(0,nrow=1,ncol=lt)
# colnames(returns.test)<-tickers
# returns.test<-as.data.frame(returns.test)
drop<-c()

for (i in 1:p)
{
  if(is.na(sum(prices.test[,i])))
  {
    drop<-c(drop,i)
  }
  if(all(prices.test[,i]==rep(0,times=ld.test)))
  {
    drop<-c(drop,i)
  }
}

r<-r[,-drop]
#lt<-length(tickers)
p<-ncol(r)
tickers<-colnames(r)
n<-length(tickers)
prices.test<-prices.test[,-drop]
returns.test<-(prices.test[ld.test,]-prices.test[1,])/prices.test[1,]



###############################
#Input investor personal views  
###############################
mu<-apply(r,2,mean)
q0<-matrix(c(0.02,0.05),ncol=1)
Omega0<-diag(c(10^(-4),10^(-4)))
Omega0.inv<-diag(c(10^(4),10^(4)))

tickers.select<-c("AAPL","FB","GOOG","MSFT")
#Find the column index where the returns are.
tickers.select<-tickers.select[tickers.select %in% tickers  ]
lt.select<-length(tickers.select)
tickers.select.index<-matrix(0,nrow=lt.select,ncol=1)
for(i in 1:lt.select)
{
  tickers.select.index[i]<-which(tickers==tickers.select[i])
}
#Fill the matrix of views according to the columns where the selected tickers are.
P<-matrix(0,nrow=2,ncol=p)
colnames(P)<-tickers
P[1,tickers.select.index[1]]<--1
P[1,tickers.select.index[2]]<-1
P[2,tickers.select.index[3]]<-1
P[2,tickers.select.index[4]]<--1

q<-15

fit.fa<-fa(r,nfactors = q,rotate = "varimax",residuals = TRUE)
fit.scores<-factor.scores(x=r,f=fit.fa)

#############
# $weights is Lambda
# a row in $scores is f_t 
##############
F.mat<-t(fit.fa$scores)
Omega<-matrix(0,nrow=p,ncol=p)
#diag(Omega)<-diag(cov(r))         #choose a good starting point for Omega based from data
nu.0<-10
s.0<-0.01
alpha.0<-nu.0/2
beta.0<-nu.0*s.0/2
Omega<-diag(rinvgamma(n=p,shape = alpha.0,scale = beta.0))
Omega.inv<-matrix(0,nrow=p,ncol=p)
diag(Omega.inv)<-1/diag(Omega)
B<-fit.fa$weights     #use the weights from the initial FA fit (might increase convergence)
c.0<-5      

###############################################################
#Estimate initial h_t from the common factors using stochvol. 
#Using those h_t, estimate V_t using observation (i) page 185 in Asai's paper (our V_t is their W_t).
##################################################################
mu.h<-rep(0,times=q)
sigma2.h<-rep(5,times=q)
phi.h.c1<-rep(10,times=q)
phi.h.c2<-rep(1.5,times=q)
eta<-rep(0.5^2,times=q)
epsilon<-matrix(0,nrow=q,ncol=T.max)
V<-matrix(0,nrow=q,ncol=T.max)

#The svsample function needs draw at least 2, otherwie we get an error in ar(), which is inside the svsample. 
#So if we choose draw=2, we will keep the last draw (because inside the Gibbs Sampler h_t is just like the 
#other variables Q_t,k,d,etc).
#If the draw>=3, we will take the average, just like in Asai on page 185 part (i).
n.draws.sv<-2
burn.sv<-10^(2)

for (i in 1:q)
{
  fit.dummy<-svsample(F.mat[i,], priormu = c(mu.h[i], sigma2.h[i]), priorphi = c(phi.h.c1[i],phi.h.c2[i]),
                      quiet = TRUE, draws = n.draws.sv, burnin = burn.sv)
  if(n.draws.sv>=3)
  {
    V[i,]<-apply(exp(fit.dummy$latent),2,mean)
    epsilon[i,]<-diag(apply(exp(-0.5*fit.dummy$latent),2,mean))%*%F.mat[i,]       #see OBS 1 in the notes 
  }
  if(n.draws.sv==2)
  {
    V[i,]<-exp(fit.dummy$latent[n.draws.sv,])
    epsilon[i,]<-diag(exp(-0.5*fit.dummy$latent[n.draws.sv,]))%*%F.mat[i,]
  }
  if(n.draws.sv<=1)
  {
    print("Please input more draws for the function svsample")
    break
  }
}

#Initial parameters using Asai's paper
gamma.0<-q
lambda.0<-5
C.0.mat<-diag(rep(1/gamma.0,times=q)) 
C.0.mat.inv<-chol2inv(chol(C.0.mat))
#A<-diag(rep(1,times=q))
A<-riwish(gamma.0,C.0.mat.inv)
#k<-25
k<-rexptr(n=1,lambda=lambda.0,range = c(q,Inf))
(d<-runif(n=1,min=-1,max=1))
#d<-0.75
A.inv<-chol2inv(chol(A))

#Initially Q_t is identity.
Q<-array(0,dim=c(q,q,T.max))
Q.inv<-array(0,dim=c(q,q,T.max))
Sigma.epsilon<-array(0,dim=c(q,q,T.max))
Sigma.epsilon.inv<-array(0,dim=c(q,q,T.max))
for(t in 1:T.max)
{
  Q[,,t]<-diag(rep(1,times=q))
  Q.inv[,,t]<-chol2inv(chol(Q[,,t]))
  Sigma.epsilon[,,t]<-diag(diag(Q[,,t])^(-1/2))%*%Q[,,t]%*%diag(diag(Q[,,t])^(-1/2))
  Sigma.epsilon.inv[,,t]<-chol2inv(chol(Sigma.epsilon[,,t]))
}


###################################
#Start the Gibbs Sampler
#We will repeatedly use the second observation made in the notes: sample from inverse wishart cause the 
#parameter for the matrix contains only the inverse of A, which can be computed outside of the for loop only once.
#Otherwise we would have to compute the inverse of Q at each time point (so T.max=836 times).
###################################

lambda.stock<-2.5
capital.init<-10^5

library(foreach)
library(doParallel)
no_cores <- 6
cl<-makeCluster(no_cores)
registerDoParallel(cl)

o1.seq<-seq(from=10^(-7),to=10^(-8),length.out = 2)
o2.seq<-seq(from=10^(-7),to=10^(-8),length.out = 3)
l.o1<-length(o1.seq)
l.o2<-length(o2.seq)

results.vec<-foreach(o1=o1.seq,.combine = "rbind")%:%
  foreach(o2=o2.seq,.combine="rbind",.packages = c("MCMCpack","mvtnorm","stochvol","armspp","quadprog"))%dopar%
  {
    # o1<-10^(-8)
    # o2<-10^(-8)
    mu.gibbs<-matrix(0,nrow=burn,ncol=p)
    w<-matrix(0,nrow=burn,ncol=p)
    Sigma.gibbs<-matrix(0,nrow=p,ncol=p)
    
    for(j in 1:burn)
    {
      #If draws paramter in svsample() is set to even =1, we get an error in the dependent function ar().
      if(n.draws.sv<=1)
      {
        print("Please input more draws for the function svsample")
        break
      }
      #1. Sample the common factors - This for loop takes the longest time - approx 20-22s
      for (t in 1:T.max)
      {
        Sigma.F<-chol2inv(chol(t(B)%*%Omega.inv%*%B+diag(V[,t]^(-0.5))%*%Sigma.epsilon.inv[,,t]%*%diag(V[,t]^(-0.5))))
        mu.F<-Sigma.F%*%t(B)%*%Omega.inv%*%t(r[t,]-mu)
        F.mat[,t]<-rmvnorm(n=1,mean=mu.F,sigma=Sigma.F)
      }
      
      #2. Sample the weights B
      for (i in 1:p)
      {
        Sigma.B<-chol2inv(chol((1/Omega[i,i])*F.mat%*%t(F.mat)+(1/(c.0^2))*diag(rep(1,times=q))))
        mu.B<-(1/Omega[i,i])*Sigma.B%*%F.mat%*%(r[,i]-rep(mu[i],times=T.max))
        B[i,]<-rmvnorm(n=1,mean = mu.B, sigma = Sigma.B)
      }
      
      #3. Sample the mean of the returns mu
      r.star<-r-B%*%F.mat
      r.star.bar<-apply(r.star,2,mean)
      Sigma.mu<-chol2inv(chol(T.max*Omega.inv +t(P)%*%Omega0.inv%*%P))
      mean.mu<-Sigma.mu%*%(T.max*Omega.inv%*%r.star.bar +t(P)%*%Omega0.inv%*%q0)
      mu<-rmvnorm(n=1,mean=mean.mu,sigma=Sigma.mu)
      
      #4. Sample the sigma^2_i for i={1,...,p} 
      e<-r-mu-B%*%F.mat
      for(i in 1:p)
      {
        Omega[i,i]<-rinvgamma(n=1,shape = alpha.0+(T.max/2), scale = beta.0+0.5*sum(e[,i]^2))
        Omega.inv[i,i]<-1/Omega[i,i]
      }
      
      #5. Stochastic volatility to find V_t and \epsilon_t
      for (i in 1:q)
      {
        fit.dummy<-svsample(F.mat[i,], priormu = c(mu.h[i], sigma2.h[i]), priorphi = c(phi.h.c1[i],phi.h.c2[i]),
                            quiet = TRUE, draws = n.draws.sv, burnin = burn.sv)
        if(n.draws.sv>=3)
        {
          V[i,]<-apply(exp(fit.dummy$latent),2,mean)
          epsilon[i,]<-diag(apply(exp(-0.5*fit.dummy$latent),2,mean))%*%F.mat[i,]       #see OBS 1 in the notes 
        }
        if(n.draws.sv==2)
        {
          V[i,]<-exp(fit.dummy$latent[n.draws.sv,])
          epsilon[i,]<-diag(exp(-0.5*fit.dummy$latent[n.draws.sv,]))%*%F.mat[i,]
        }
      }
      
      #6. Sample Q_t using Metropolis-Hastings: sample a candidate Q from Inverse Wishart (see the notes for the reduced ratio)
      #Also, make the loops more efficient. Now we only compute Q.inv and Sigma.epsilon and Sigma.epsilon.inv when a candidate 
      #in Metropolis Hastings step is accepted
      count.mh<-0
      for (t in 2:(T.max-1))
      {
        Q.cand<-riwish(k+1,k*power.matrix(Q[,,t-1],d/2)%*%A.inv%*%power.matrix(Q[,,t-1],d/2)+epsilon[,t]%*%t(epsilon[,t]))
        rho<-ratio.mh.t.optim.log(k=k,d=d,Qt.cur = Q[,,t],Qt.cur.inv = Q.inv[,,t],A.inv=A.inv, Qtp1.cur.inv = Q.inv[,,t+1], 
                                  Qt.cand = Q.cand ,epsilont = epsilon[,t], Sigma.epsilont.cur.inv = Sigma.epsilon.inv[,,t])
        u<-runif(n=1,min=0,max=1)
        if (log(u)<rho)
        {
          Q[,,t]<-Q.cand
          Q.inv[,,t]<-chol2inv(chol(Q.cand))
          Sigma.epsilon[,,t]<-diag(diag(Q[,,t])^(-1/2))%*%Q[,,t]%*%diag(diag(Q[,,t])^(-1/2))
          Sigma.epsilon.inv[,,t]<-diag(diag(Q[,,t])^(1/2))%*%Q.inv[,,t]%*%diag(diag(Q[,,t])^(1/2))
          count.mh<-count.mh+1
        }
      }
      Q.cand<-riwish(k+1,k*power.matrix(Q[,,T.max-1],d/2)%*%A.inv%*%power.matrix(Q[,,T.max-1],d/2)+epsilon[,T.max]%*%t(epsilon[,T.max]))
      rho.Tmax<-ratio.mh.T.optim.log(QT.cur.inv = Q.inv[,,T.max],QT.cand = Q.cand, epsilonT = epsilon[,T.max], 
                                     Sigma.epsilonT.cur.inv = Sigma.epsilon.inv[,,T.max])
      u<-runif(n=1,min=0,max=1)
      if(log(u)<rho.Tmax)
      {
        Q[,,T.max]<-Q.cand
        Q.inv[,,T.max]<-chol2inv(chol(Q.cand))
        Sigma.epsilon[,,T.max]<-diag(diag(Q[,,T.max])^(-1/2))%*%Q[,,T.max]%*%diag(diag(Q[,,T.max])^(-1/2))
        Sigma.epsilon.inv[,,T.max]<-diag(diag(Q[,,T.max])^(1/2))%*%Q.inv[,,T.max]%*%diag(diag(Q[,,T.max])^(1/2))
        count.mh<-count.mh+1
      } 
      
      # 7. Sample A
      Sigma.A<-C.0.mat.inv
      for(t in 2:T.max)
      {
        Sigma.A<-Sigma.A+k*power.matrix(Q[,,t-1],d/2)%*%Q.inv[,,t]%*%power.matrix(Q[,,t-1],d/2)
      }
      A<-riwish(gamma.0+k*T.max, Sigma.A)
      A.inv<-chol2inv(chol(A))
      
      #memorize the \sum_{t=1}^T Q_{t-1}^{d/2}Q_t^{-1}Q_{t-1}^{d/2} because it shows in 2 other distributions 
      #and we do not want to waste time computing it multiple times.
      sum.quad.Q.d<-(Sigma.A-C.0.mat.inv)/k 
      
      #9. Sample d. We need arms for log concave sampling 
      
      d<-arms(n_samples = 1,
              log_pdf = function(d)
              {
                sum.dummy<-0
                for(t in 1:(T.max-1))
                {
                  sum.dummy<-sum.dummy+log(det(Q[,,t]))
                }
                res<-(k*d/2)*sum.dummy-(1/2)*sum(diag(k*A.inv%*%sum.quad.Q.d))
                return(res)
              },
              lower=-1,
              upper=1)
      
      #9. Sample k. We need arms again - this also takes quite a long time - 7-9s
      
      k<-arms(n_samples = 1,
              log_pdf = function(k)
              {
                sum.dummy1<-0
                for (i in 1:q)
                {
                  sum.dummy1<-sum.dummy1+log(gamma((k+1-i)/2))
                }
                sum.dummy1<-T.max*sum.dummy1
                res1<--lambda.0*k+((T.max*k*q)/2)*log(k/2)-(T.max*k/2)*log(det(A))-sum.dummy1
                sum.dummy2<-0
                for(t in 2:T.max)
                {
                  sum.dummy2<-sum.dummy2+log(det(power.matrix(Q[,,t-1],d/2)%*%Q.inv[,,t]%*%power.matrix(Q[,,t-1],d/2)))
                }
                res2<-(k/2)*sum.dummy2-(k/2)*sum(diag(A.inv%*%sum.quad.Q.d))
                return(res1+res2)
              },
              lower = q-1,
              upper = T.max)
      
      mu.gibbs[j,]<-mu
      Sigma.gibbs<-Sigma.gibbs+B%*%diag(V[,T.max]^0.5)%*%diag(diag(Q[,,T.max])^(-0.5))%*%Q[,,T.max]%*%
        diag(diag(Q[,,T.max])^(-0.5))%*%diag(V[,T.max]^0.5)%*%t(B)+Omega
      #qp<-solve.QP(Dmat = Sigma.gibbs, dvec = mu, Amat = matrix(rep(1,times=p),nrow=p,ncol=1), bvec = 1, meq = 1)
      #w[j,]<-qp$solution
      
      # print(paste0("Nr.MH accepted: ", count.mh, " out of ", T.max-1 ," ,MH at T: ", rho.Tmax,
      #              " burn period completion: ", 100*round(j/burn,digits = 4), "% at time ", Sys.time()))
      cat(paste("\n","Completed o1=",o1,"and o2=", o2," total: " ,round(j/burn,digits=4)*100, "% and nr.MH accepted: ", count.mh, " out of ", T.max-1 
                ," at ",Sys.time(),"\n"),file = "FB error parallel industries 10^-7 10^-8 10^-7 10^-8.txt",append = T)
      
    }
    
    # w.alt<-solve.QP(Dmat = sigma.model.post.mean, dvec = mu.post.mean, Amat = matrix(rep(1,times=n),nrow=n,ncol=1), bvec = 1, meq = 1)$solution
    # position.alt<-floor(w.alt*capital.init/prices.test[1,])
    # capital.stock.init.alt<-position.alt*prices.test[1,]
    # capital.stock.fin.alt<-position.alt*prices.test[ld.test,]
    
    colnames(mu.gibbs)<-tickers
    Sigma.post<-Sigma.gibbs/burn
    mu.post<-apply(mu.gibbs,2,mean)
    w.org<-(1/lambda.stock)*solve(Sigma.post)%*%mu.post
    w.org<-w.org/abs(sum(w.org))
    position<-floor(w.org*capital.init/prices.test[1,])
    capital.stock.init<-position*prices.test[1,]
    capital.stock.fin<-position*prices.test[ld.test,]
    
    distance<-norm(P%*%mu.post-q0,type="2")
    
    c(o1,o2,sum(capital.stock.fin-capital.stock.init),distance)
    
  }



colnames(results.vec)<-c("o1","o2","profit","distance")

write.table(results.vec,"FB error parallel industries 10^-7 10^-8 10^-7 10^-8.csv",sep=",",row.names = FALSE)

sink(type="message")
close(err)
stopCluster(cl)






