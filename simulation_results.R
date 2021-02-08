rm(list=ls())
packages <- c("cccp","ggplot2","tidyr","dplyr","latex2exp","purrr","furrr","future",
              "gridExtra","grid","DT","pushoverr","Rcpp","expss","scales","cumstats","stargazer","xtable","selectiveInference" )

lapply(packages, require, character.only = TRUE)


Rcpp::sourceCpp("zero_terms_function.cpp")

raw_storage <- data.frame(method = NA,value =NA, n = NA)
#for (n in c(100)) {
for (n in c(100,300,500)) {

  ### set parameters
  n_simulation <- 100
  p = n
  tau2=1
  prop_large = 0.6
  tau2_B=prop_large*tau2
  frac = 0.5  #fraction of observations to use for selection
  smp_size <- floor(frac * n) 
  c = 4       #number of large betas
  K=10       #number of splits to average over
  large = sqrt(tau2_B/c)
  small = sqrt((tau2-tau2_B)/(p-c) )  
  beta_range<-c(rep(large, c), rep(small,p-c))
  sigma2<-1
  s=sum(beta_range^2)+sigma2
  
  ###################################
  ### optimal coefficients table
  df <- data.frame(j_1=rep(1:p,each=p), j_2=rep(1:p,times=p))%>%
    group_by(j_1,j_2)%>%  mutate(beta_jj  = beta_range[j_1]*beta_range[j_2])
  
  
  ##########################################################################
  #                 define functions
  ###########################################################################
  
  
  ###       Oracle Estimator
  T_oracle_f<-function (X) { 
           df<- df %>% mutate(mean_X_jj =  mean(X[,j_1]*X[,j_2]), 
           g_jj = ifelse(j_1==j_2,mean_X_jj - 1,mean_X_jj), psi = beta_jj * g_jj  )
    return( sum(df$psi)) 
  }
  ######################################
  ###      Proposed Estimator
  proposed_estimator <- function(X,Y) { 
    dt<- selection_algoritm(X,Y)   
    estimated_indexes <- filter(dt,pred == "big") %>% dplyr::select(j)
    df <- df %>%    filter(j_1 %in% estimated_indexes$j & j_2 %in% estimated_indexes$j)%>%
      group_by(j_1,j_2)%>%
      mutate(   psi_hat =  ifelse(j_1==j_2,func(W[,j_1],W[,j_2], X[,j_1]*X[,j_2]-1),
                                  func(W[,j_1],W[,j_2],X[,j_1]*X[,j_2])
      )  ) 
    return(sum(df$psi_hat))    
  }  
  ######################################
    ###  Modified Estimator
  modified_estimator <- function(X,Y) {
          
          ind <- sample(seq_len(n), size = smp_size) # randomly split the sample into two equal subsets:
        # covariate-selection subset and a subset for estimation of zero-estimators
          Y_1 <- Y[ind, ] 
          X_1 <- X[ind, ]
          X_2 <- X[-ind, ]  
          W_2 <- W[-ind, ]
          dt<- selection_algoritm(X_1,Y_1) #apply selection algorithm
          estimated_indexes <- filter(dt,pred == "big") %>% dplyr::select(j)
          df <- df %>%    filter(j_1 %in% estimated_indexes$j & j_2 %in% estimated_indexes$j)%>%
          group_by(j_1,j_2)%>%
          mutate(   psi_hat =  ifelse(j_1==j_2,func(W_2[,j_1],W_2[,j_2], X_2[,j_1]*X_2[,j_2]-1),
                                  func(W_2[,j_1],W_2[,j_2],X_2[,j_1]*X_2[,j_2])
      )  ) 
    return(sum(df$psi_hat))    
  }  
  ######################################
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
    ###############################################
  #                Simulations   
   ################################################
  # Now we simulate independent datasets (X,Y), calculate the estimators and store the results
  # into a dataframe.
  selected_covariates_c=NA
  ncol_ = 5
  sim_results<- as.data.frame(matrix( NA, nrow = n_simulation, ncol = ncol_)) # create an empty df to store simultaions
  colnames(sim_results) <- c("Naive", "Oracle","Proposed","Modified","PSI")
    for (i in 1:n_simulation) {
    X_sample<-matrix(   rnorm(n*p)    ,nrow = n, ncol = p)
    y_sample<- X_sample%*%beta_range + rnorm(n,mean=0,sd=sqrt(sigma2))  # basic model
    
    ### calculate naive estimator
        Y=y_sample
        X=X_sample
        W = as.data.frame(X)*Y
        mean_squared_W<- colMeans(W^2)  #calculate first element of beta_square_hat
        var_W<-sapply(W,var) # #calculate the second element of beta
        beta_square_hat<- mean_squared_W - var_W # calculate vector of beta2_hat
        naive<-sum(beta_square_hat)
   
    ### calculate PSI estimator
        PSI <-   estimateSigma(X, Y, intercept=TRUE, standardize=TRUE)   # estimate estimateSigma from selectiveInference package
       
    ### calculate Oracle estimator    
       T_oracle <- T_oracle_f(X)
    ### calculate Oracle estimator    
      T_proposed_estimator <- proposed_estimator(X,Y)
  
    ########################################
    ### calculate the Modified estimator 
    ######################################## 
    #refitting step: calculate psi_{jj'} elements
      phi=NA
      for (k in 1:K) {
          phi[k] <-  modified_estimator(X,Y)
      }
      Modified <- naive - 2*sum(phi)/K
    ################################################  
     ### Store all estimators in a dataframe:    
      sim_results$Naive[i]      <- naive
      sim_results$Oracle[i]   <- naive     -2 * T_oracle
      sim_results$Proposed[i] <- naive - 2 * T_proposed_estimator
      sim_results$Modified [i]   <- Modified  
      sim_results$PSI[i]<-  var(as.vector(Y)) - PSI$sigmahat^2
  } # this is the end of simulation for-loop.


    ### reshape and bind data
    compare<-gather(sim_results,"method","value",1:ncol_)
    raw_storage <- bind_rows(raw_storage,cbind(compare,n))   # save raw data for boxplot later on
}  #end of for-loop over different values of n.

raw_storage = raw_storage[-1,] %>% mutate(tau2=tau2) # remove 1st row which is NA by construction

###############################################
#    Boxplot 
##############################################
g_boxplot<-raw_storage %>%
  mutate(method=factor(method, levels=c(
    "Naive"
    ,"Modified"
    ,"Proposed"
    ,"Oracle"
    ,"PSI"
    ) )) %>%ggplot( aes(x=as.factor(n), y=value
            ,fill=method)) + geom_boxplot() + theme_bw() +
             ylab("Value")+ labs(fill = "Estimator") + xlab("Sample Size")+
             geom_hline(yintercept=tau2, linetype="dashed", color = "red")+
             scale_fill_brewer(palette="Pastel1") +ggtitle("Estimators")

###############################################
#    Summary Table
##############################################
stat<- raw_storage %>% group_by(method,n) %>%
       summarise(Mean =round(mean(value),2)
            ,SE = round(sd(value),2)
            ,MSE=mean((value-tau2)^2) 
            ,SE_RMSE = round(sqrt(var((value-tau2)^2)/(4*MSE*n_simulation)),3)*1000
              ) %>% mutate(RMSE = round(sqrt(MSE),2), bias = round(tau2-Mean,2), tau2=round(tau2,1),
            method=factor(method, levels=c("Naive","Modified","Proposed","Oracle","PSI")))%>%
            arrange(n,method) %>% dplyr::select(tau2,n,Estimator = method,Mean,bias
             ,SE, RMSE, SE_RMSE )

#########################################################################

#outputs of simulation: boxplot and summary-statistics table.
g_boxplot
stat
