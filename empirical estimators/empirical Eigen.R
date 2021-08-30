rm(list = ls())
packages <- c("ggplot2","tidyr","dplyr","latex2exp","gridExtra","pushoverr","Rcpp","scales","stargazer","selectiveInference" )
library(cccp)
lapply(packages, require, character.only = TRUE)
source("run_sim.R")
Rcpp::sourceCpp("zero_terms_functions.cpp")
start_time <- Sys.time()


sim_fn <- function(setting_param_list){
  Rcpp::sourceCpp("zero_terms_functions.cpp")
  
  EigenPrism <- function(y,X,invsqrtSig=NULL,alpha=0.05,target='beta2',zero.ind=c(),diagnostics=T){
    # Author: Lucas Janson (statweb.stanford.edu/~ljanson)
    # Runs EigenPrism procedure for estimating and generating confidence
    #  intervals for variance components in high-dimensional linear model:
    #       y = X%*%beta + e,   rows of X iid~ N(0,Sig),   e iid~ N(0,sigma^2)
    #  Requires cccp package for solving second order cone optimization.
    #  Note confidence interval endpoints may lie outside parameter domain, so it may be appropriate
    #   to clip them after the fact.
    # 
    # Inputs:
    #  y: response vector of length n (will automatically be centered)
    #  X: n by p design matrix; columns will automatically be centered and scaled to variance 1;
    #      should not contain intercept column, since both y and X will be centered
    #  invsqrtSig: if columns of X not independent, p by p positive definite matrix which is the square-root
    #               of the inverse of Sig, where Sig is the *correlation* matrix of the X (default is identity)
    #  alpha: significance level for confidence interval (default = 0.05)
    #  target: target of estimation/inference
    #		  'beta2' (default) is the squared 2-norm of the coefficient vector: sum(beta^2)
    #           'sigma2' is the noise variance sigma^2
    #           'heritability' is the fraction of variance of y explained by X%*%beta: t(beta)%*%Sig%*%beta/var(y)
    #  zero.ind: vector of which indices of the weight vector w to constrain to zero (default is none)
    #  diagnostics: boolean (default = T) for whether to generate diagnostic plots for the V_i, lambda_i, and w_i
    #  
    # Outputs:
    #  estimate: unbiased estimate of the target (for heritability, only approximately unbiased)
    #  CI: 100*(1-alpha)% confidence interval for target
    
    # Get dimensionality of problem
    n = nrow(X)
    p = ncol(X)
    
    # Transform y and X to proper form
    y = y-mean(y)
    X = scale(X,T,T)*n/(n-1)
    if(!is.null(invsqrtSig)) X = X%*%invsqrtSig
    
    # Take singular value decomposition and rescale singular values
    svd = svd(X)
    lambda = svd$d^2/p
    
    # Defined cone-constrained linear problem to optimize weights; [v; w] is vector of optimization variables
    q = c(1,rep(0,n)) #coefficient vector in objective function
    A = rbind(c(0,rep(1,n)),c(0,lambda)) #matrix for linear constraints
    b = c(0,1) #vector for linear constraints
    if(target=='sigma2') b = c(1,0) #switch constraints if target is sigma^2
    # Constrain some weights to be zero if desired
    if(!is.null(zero.ind)){
      A = rbind(A,cbind(rep(0,length(zero.ind)),diag(rep(1,n))[zero.ind,]))
      b = c(b,rep(0,length(zero.ind)))
    }
    # Define second-order cone constraints
    soc1 = socc(diag(c(1/4,rep(1,n))),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
    soc2 = socc(diag(c(1/4,lambda)),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
    prob = dlp(as.vector(q),A,as.vector(b),list(soc1,soc2))
    
    # Solve optimization problem and extract variables
    opt = cps(prob,ctrl(trace=F))
    v = getx(opt)[1]
    w = getx(opt)[-1]
    
    # Compute estimate and y's variance
    est = sum(w*(t(svd$u)%*%y)^2)
    yvar = sum(y^2)/n
    
    # Compute confidence interval
    CI = est + yvar*sqrt(v)*qnorm(1-alpha/2)*c(-1,1)
    if(target=='heritability'){
      est = est/yvar
      CI = CI/yvar
    }
    
    # Generate list with results
    result=list()
    result$estimate = est
    result$CI = CI
    
    # Generate diagnostic plots
    if(diagnostics){
      par(mfrow=c(1,3))
      
      # Check that eigenvectors are approximately Gaussian
      nV = floor(log10(n))
      srtV = svd$v[,10^(0:nV)]
      labs = c()
      for(i in 1:(nV+1)){
        srtV[,i] = sort(srtV[,i])
        ind = 10^(i-1)
        labs = c(labs,bquote(V[.(ind)]))
      }
      matplot(qnorm((1:p)/(p+1)),srtV,type="l",lwd=2,
              ylab="Quantiles of Eigenvectors",xlab="Gaussian Quantiles",
              main=expression(paste("Check Gaussianity of Eigenvectors ",V[i])))
      legend("topleft",as.expression(labs),col=1:(nV+1),lty=1:(nV+1),lwd=2)
      
      # Check that there are no outliers in the eigenvalues
      hist(lambda,main=expression(paste("Histogram of Normalized Eigenvalues ",lambda[i])),
           xlab=expression(lambda[i]))
      
      # Check that the weights are not dominated by just a few values
      srtw = sort(abs(w),T)
      plot(1:n,cumsum(srtw)/sum(srtw),type="l",lwd=2,
           main=expression(paste("Fraction of Total Weight in Largest k ",w[i])),
           xlab="k",ylab="Fraction of Total Weight")
    }
    
    return(result)
  }
  
  
  seed <-  setting_param_list$seed
  set.seed(seed)
  prop_large <-  setting_param_list$prop_large
  n <- setting_param_list$n
  p = n
  
  tau2_B=prop_large*tau2
  large = sqrt(tau2_B/c)
  small = sqrt((tau2-tau2_B)/(p-c) )  
  beta_range<-c(rep(large, c), rep(small,p-c))
  
  
  ################################################################
  ###### semi-supervised data of X
  ##########################################################
  X_pop<-matrix(
    #rnorm(n_pop*p)
    1-rexp(n_pop *p,1)
    ,nrow = n_pop, ncol = p)
  conditional_mean_p <- as.matrix(X_pop)%*%beta_range
  mehane <- var(double_dist_sum(X_pop))
  
  
  # sample data
  X_sample<-matrix(   
    #rnorm(n*p)   
    1-rexp(n *p,1) 
    ,nrow = n, ncol = p)
  conditional_mean <- as.matrix(X_sample)%*%beta_range 
  y_sample<- conditional_mean + rnorm(n,mean=0,sd=sqrt(sigma2))  # basic model
  
  Y=y_sample
  X=X_sample
  
  Eigen_tau2 <- EigenPrism(y=Y,X=X,target='beta2',diagnostics=F)$estimate  
  
  X<- as.matrix(X)
  g_int <- mean(double_dist_sum(X))
  ###   Bootstraping step 
  Eigen_Pri_b = g_int_b  = NA
  for (b in 1:B) {    
    # step 3.1: Resample with replacement
    indices <- sample(1:nrow(X),size = nrow(X) ,replace=TRUE) 
    X_b <- X[indices,]
    Y_b <- Y[indices]
    # step 3.2: calculate initials estimator:  EigenPrism
    Eigen_Pri_b[b] <- EigenPrism(y=Y_b,X=X_b,target='beta2',diagnostics=F)$estimate
    g_int_b[b] <- mean(double_dist_sum(X_b))
  }# end of bootstrapping loop
  
  cov_emp =cov(Eigen_Pri_b, g_int_b)
  c_star_mean <- n*cov_emp/mehane 
  
  
  return(tibble( 
    Eigenprism = Eigen_tau2
    ,"Empirical Eigen" = Eigen_tau2 - c_star_mean*g_int
    
  ))
}

######################################################################################

########################################
#### initial parameters settings
######################################
#availableCores()
n_cores =72  # my laptop has 8 cores.  Shimrit server has 72. 
set.seed(1)
n_simulation = 100
B=100
n_pop<-3000
tau2=1
sigma2=1
c = 5     #number of large betas

###############################################################################

setting_param_list <-  list(seed = 1:n_simulation, prop_large = c(0.33,0.66, 0.99),n=400)

# availableCores()
result <-run_sim(sim_fn, setting_param_list, nworkers = n_cores)
arrange(result, scenario)
raw_storage<-gather(result, method, value, 5:6) 
end_time <- Sys.time()
total_time  <- round(end_time - start_time)


###############################################
#    Boxplot 
##############################################
g_boxplot<-raw_storage %>%
  mutate(method=factor(method, levels=c(
    "Eigenprism"
    ,"Empirical Eigen"
    
  ) )) %>%ggplot( aes(x=as.factor(prop_large), y=value
                      ,fill=method)) + geom_boxplot() + theme_bw() +
  ylab("Value")+ labs(fill = "Estimator") + xlab(TeX("$\\tau^2_{B}$"))+
  geom_hline(yintercept=tau2, linetype="dashed", color = "red")+
  scale_fill_brewer(palette="Pastel1") 

###############################################
#    Summary Table
##############################################

stat<- raw_storage %>% group_by(method,prop_large,n) %>%
  summarise(Mean =round(mean(value),2)
            ,SD = round(sd(value),3)
            ,MSE=mean((value-tau2)^2) 
            ,SD_RMSE = round(sqrt(var((value-tau2)^2)/(4*MSE*n_simulation)),3)*1000
  ) %>% mutate(RMSE = round(sqrt(MSE),3), bias = round(tau2-Mean,2), tau2=round(tau2,1), n=n
  )%>% arrange(desc(prop_large))


stat <- stat %>%   dplyr::select(prop_large,n = n, Estimator = method,Mean,bias
                                 ,SD
                                 ,RMSE
                                 , SD_RMSE
)%>%
  arrange(prop_large) %>%
  mutate(tau2=round(tau2,1), sigma2=round(sigma2,1),
         n_sim=n_simulation
         , "# large betas" =c
         ,prop_large =  label_percent()(prop_large),
         total_time=total_time,  n_pop = n_pop )


g_boxplot
stat
