#######################################################
#         Data description                            #
#######################################################
j<-SN # unique individual identifier
n<-length(j) # Total number of individuals
i<-TwinID  # Cluster identifier
G<-length(unique(i) ) # Number of clusters/groups.
delta<-CENSORING  
time<-as.matrix(Time) # This is individual time to follow up
X<- as.matrix(cbind(covariates) ) # Covariate matrix
 m<-as.numeric(ncol(X))
#########################################################
# Marginal Likelihood estimation of parameters assuming #
# gamma frailty and splines for baseline hazard         #
#########################################################
#Number of  failures in each cluster
Di<-aggregate(delta, by=list ( i ) ,FUN=sum) #Total  number of observed failures in each cluster
e<-colSums(Di[ 2 : ( r+1) ] ) #e gives the sum of uncensored observation in all the clusters per dataset 
di<-Di [ 2 : ( r+1) ] #1st column of Di is pair identifier so needs to be omitted
#========================================================
Theta_Weib<-vector( "numeric" , length(1:r))
Beta_Weib<-vector( "numeric" , length(1:r))
#============================
#df = knots(K)+degree(D)+1 
alpha<-vector( "numeric" , df) #B Spline coefficients
#Let m be the number of covariates in dataset.
# p[1:m]<-beta ; p[m+1]<-theta ; p[m+2:m+df+1] <- B-spline coefficients
#Example, suppose we have 4 covariates, 7 knots for a cubic spline
# then we will have p[1:4]<-beta ; p[5]<-theta ; p[6:16] <- B-spline coefficients
likelihood.Weibull2<-function(p)
{
  baselinehaz <- bs(time, df, intercept=TRUE)   %*% 
    as.matrix(c( exp(p[m+2:m+df+1]) ) )
    cumbaselinehaz <- ibs(time, df, intercept=TRUE)   %*% 
    as.matrix(c( exp(p[m+2:m+df+1]) ) )
    cumhaz<-    cumbaselinehaz * exp(X%*%matrix(c(p[1:m])) )
  cumhaz<-as.numeric(sapply(split(cumhaz,i),sum))# Sum of j from 1 to ni. i.e cluster sum
    lnhaz<-delta*(  X%*%matrix(c(p[1:m]))  +log(baselinehaz))
  lnhaz<-as.numeric(sapply(split(lnhaz,i),sum))
 
  lik<-e*log( (p[m+1])  )- #1st term
    sum((di+1/((p[m+1])) )*log(1+cumhaz*( (p[m+1])  ) ))+ #4th term
    sum(lnhaz)- ##5th term
    G*log(gamma(1/( (p[m+1])  ) ) )+ #2nd term
    sum(log(gamma(di+1/( (p[m+1])) ) )  )# Third term
  -lik
}
  #Obtaining the MLE's
  initial2<-c(rep(0.5,m+1), rep(0.1,df)   ) # Initial values of beta,theta,B-spline coefficients.
  t2<- nlminb(initial2,likelihood.Weibull2,hessian=T, control = list(rel.tol = 1e-6) )
  #Parameter estimates
  std3 <- sqrt(diag(solve(hessian(likelihood.Weibull2, t3$par))))
  Beta.SE<-data.frame(t2$par[1:m] , std3[1:m]) # Beta output+standard error
  theta.se<-data.frame(t2$par[m+1] , std3[m+1]) #Theta output+standard error
  alpha =exp( t2$par[m+2:m+df+1]) 
#=================================================
#xx<-as.numeric(seq(0, 3, by = 0.1))
#Predicting cumulative baseline hazard at time tt for puroses of computation of
#survival probability
Int2<-predict(ibs(time, df=11, intercept = TRUE), tt)  
H0tt<-alpha%*%t(Int2) #Cumulative Baseline hazard
#Predicting survival probability at time tt for individual with covariate x_ij
#Let X_ij= X[1,]
Sp1WeibX<-(1+theta.se[1]*H0tt*exp( X[1,]%*% as.matrix(cbind( Beta.SE[1] ) )  ) )^(-(1/ theta.se[1]));Sp1WeibX



