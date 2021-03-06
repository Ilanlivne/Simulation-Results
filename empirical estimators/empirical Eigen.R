rm(list=ls())
packages <- c("cccp","ggplot2","tidyr","dplyr","latex2exp","purrr","furrr",
              "gridExtra","grid","DT","pushoverr","Rcpp","expss","scales","cumstats","stargazer","xtable","selectiveInference" )


lapply(packages, require, character.only = TRUE)

start_time <- Sys.time()
library(tidyr)
library(tictoc)



norm_square <- function(x) (sum(x^2))
#set.seed(60)
raw_storage <- data.frame(method = NA,value =NA, n = NA)

for (n in c(100,300,500)) {
  # for (n in c(40,50)) {
  
  
  ### set parameters    
  n_simulation<-1000
  p = n
  B=100  # number of bootstrap samples
  tau2=2
  prop_large = 0.7   # proportion of "large betas" from tau2
  tau2_B=prop_large*tau2
  frac = 0.5  #fraction of observations to use for selection
  c = 4 # number of large betas
  large = sqrt(tau2_B/c)
  small = sqrt((tau2-tau2_B)/(p-c) )  
  beta_range<-c(rep(large,c),rep(small,p-c))
  sigma2<-1
  
  ###################################
  ###  coefficients' combination table
  
  df <- data.frame(j_1=rep(1:p,each=p), j_2=rep(1:p,times=p)) %>% 
    mutate( ind = ifelse(j_1==j_2,1,0)
    )%>%group_by(j_1,j_2) 
  
  ###################################################  
  
  
  
  ##########################################################################
  #                 define functions
  ###########################################################################
  
  ###  empirical estimator
  emp_estimator_f<-function (X) {
    bind_df<-bind_df %>% 
      filter(j_1 %in% estimated_indexes & j_2 %in% estimated_indexes)%>%
      group_by(j_1,j_2)%>% mutate(ind = ifelse(j_1==j_2,1,0))%>% 
      mutate( g = ifelse(ind == 1, mean(X[,j_1]^2-1), mean(X[,j_1]*X[,j_2]) ),
              a_g = c_est*g)
    return(sum(bind_df$a_g))
  }
  ###############################################################################
  
      
  ###     EigenPrism Estimator
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
  ######################################################################################################
  
    ### covariate selection algorithm (Algorithm 3)
  selection_algoritm <- function(X,Y){
    # step 1: calculate \hat\beta_j^2 for j=1,...,p:  
          W_3 = as.data.frame(X)*Y
          mean_squared_W<- colMeans(W_3^2)  #calculate first element of beta_square_hat
          var_W<-sapply(W_3,var) # #calculate the second element of beta
          beta_square_hat<- mean_squared_W- var_W # calculate vector of beta2_hat
    # step2:  calculate the differences lamda_j for j=2,...,p:  
          dt <- data.frame(
            beta = beta_range, beta_type = c(rep("big",c),rep("small",p-c)),
            j = 1:p,  beta_square_hat   ) %>% arrange(beta_square_hat) %>%
            mutate( index = 1:n(),  lag_1 = lag(beta_square_hat), 
                    diff = beta_square_hat - lag_1
            )%>%dplyr::select(-lag_1) %>% filter(diff != "NA")  
    #step 3: Select the covariates of B_gamma.
        #calculate  j_star:
          index_star = dt %>%  mutate(max_diff = max(diff, na.rm = T))%>% 
            filter(diff == max_diff) %>% dplyr::select(index) %>% unlist()
          if (p-index_star<2) {
            index_star = index_star  %>% unlist()-1
          }
      dt<- dt %>% mutate(pred =if_else(index >= index_star,"big","small") )
    return(dt)
  }
  
  ####################################
  # Algorithm 4
  ###################################
  
  # Input: a dataset (X,Y) from which I will resample the bootstrap samples.
  X<- matrix(rnorm(n*p),nrow = n, ncol = p)
  Y <-  X %*% beta_range + rnorm(n,mean=0,sd=sqrt(sigma2))  # basic model
  
  ###   Bootstraping step 
  boot_dt <-list()  # create empty list to store results from bootstrap samples
  for (b in 1:B) {    
    # step 3.1: Resample with replacement
    indices <- sample(1:nrow(X),size = nrow(X) ,replace=TRUE) 
    X_b <- X[indices,]
    Y_b <- Y[indices]
    # step 3.2: calculate initials estimator:  EigenPrism
    Eigen_tau2 <- EigenPrism(y=Y_b,X=X_b,target='beta2',diagnostics=F)$estimate
    
    ##################################################
    # step 3.3: For each bootstrap sample, calculate the zero-estimators g_{jj'} and store the result in
    # a list
    boot_dt[[b]]<- df %>% as_tibble() %>% group_by(j_1,j_2) %>%
      mutate( g = ifelse(ind == 1, mean(X_b[,j_1]^2-1), mean(X_b[,j_1]*X_b[,j_2]) ),
              Eigen =  Eigen_tau2   )    
  }# end of bootstrapping loop
  ####################################################################################
  ### step 4: approximate the coefficients c{jj'}:
  
  bind_df<- bind_rows(boot_dt)%>% 
    mutate(cov_emp = cov(Eigen,g), c_est = -cov_emp*n/2 )  %>%
    dplyr::select(j_1,j_2,c_est) %>%
    summarise( c_est = max(c_est) )
  
  
  ###############################################
  #                Simulations   
  ################################################
  # Now we simulate independent datasets (X,Y), calculate the estimators and store the results
  ncol_ = 2
  sim_results<- as.data.frame(matrix( NA, nrow = n_simulation, ncol = ncol_)) # create an empty df to store simultaions
  colnames(sim_results) <-  c("Empirical Eigen","Eigen")
  for (i in 1:n_simulation) {
    X<-matrix( rnorm(n*p) ,nrow = n, ncol = p)        
    Y<- X%*%beta_range + rnorm(n,mean=0,sd=sqrt(sigma2))  # basic model
    
    # Calculate the EigenPrism estimator
     Eigen_tau2 <- EigenPrism(y=Y,X=X,target='beta2',diagnostics=F)$estimate  
    
    ##############################################
    ### covariate selection step:
    dt<-selection_algoritm(X,Y)
    estimated_indexes <- filter(dt,pred == "big") %>% dplyr::select(j)
    estimated_indexes<- estimated_indexes[,1]  # extract thecovariates from the covariates-selection algorith 
    
    ### Store all estimators in a dataframe:
    sim_results$Eigen[i] <- Eigen_tau2
    sim_results$"Empirical Eigen"[i] <- Eigen_tau2 + emp_estimator_f(X)
  }# end of simulation for-loop.
  
  
  ### reshape and bind data
  compare<-gather(sim_results,"method","value",1:ncol_)
  raw_storage <- bind_rows(raw_storage,cbind(compare,n))   # save raw data for boxplot 
  
  
  
}  #end of for-loop over different values of n.


raw_storage = raw_storage[-1,]%>% mutate(tau2=tau2)  # remove 1st row which is NA by construction



###############################################
#    Summary Table
##############################################

stat<- raw_storage %>% group_by(method,n) %>%
  summarise(Mean =round(mean(value),2)
            ,SE = round(sd(value),2)
            ,MSE=mean((value-tau2)^2) 
            ,SE_RMSE = round(sqrt(var((value-tau2)^2)/(4*MSE*n_simulation)),3)*1000
  ) %>% mutate(RMSE = round(sqrt(MSE),2), bias = round(tau2-Mean,2), tau2=round(tau2,1),
               method=factor(method, levels=c("Eigen","Empirical Eigen")))%>%
  arrange(n,method) %>% dplyr::select(tau2,n,Estimator = method,Mean,bias
                                      ,SE, RMSE, SE_RMSE )

#########################################################################

# Output of simulation: summary-statistics table 
stat




# Barplot
max_n <-max(stat$n)
stat_max_n<-filter(stat,n==max_n)
ggplot(stat_max_n,aes(x=Estimator,y=SE,fill=as.factor(Estimator)))+
  geom_bar(stat="identity",position = "dodge")+
  theme(aspect.ratio = 0.6)

