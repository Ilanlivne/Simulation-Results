rm(list = ls())
packages <- c("ggplot2","tidyr","dplyr","latex2exp","gridExtra","pushoverr","Rcpp","scales","stargazer","selectiveInference" )
lapply(packages, require, character.only = TRUE)
start_time <- Sys.time()
source("run_sim.R")
Rcpp::sourceCpp("zero_terms_functions.cpp")



sim_fn <- function(setting_param_list){
  Rcpp::sourceCpp("zero_terms_functions.cpp")
  
  
  seed <-  setting_param_list$seed
  set.seed(seed)
  prop_large <-  setting_param_list$prop_large
  p = n
  
  tau2_B=prop_large*tau2
  large = sqrt(tau2_B/c)
  small = sqrt((tau2-tau2_B)/(p-c) )  
  beta_range<-c(rep(large, c), rep(small,p-c))
  
  ############################################
  ### covariate selection algorithm 
  ##############################################
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
  ###############################
  #Some more functions
  ###############################
  ###       Oracle Estimator
  T_oracle_f<-function (X) { 
    df<- df %>% dplyr::mutate(mean_X_jj =  mean(X[,j_1]*X[,j_2]), 
                              g_jj = ifelse(j_1==j_2,mean_X_jj - 1,mean_X_jj), psi = beta_jj * g_jj  )
    return( sum(df$psi)) 
  }
  #############################################################
  ###      Selection Estimator
  Selection_estimator <- function(X,Y) { 
    dt<- selection_algoritm(X,Y)   
    estimated_indexes <- filter(dt,pred == "big") %>% dplyr::select(j)
    df <- df %>%    filter(j_1 %in% estimated_indexes$j & j_2 %in% estimated_indexes$j)%>%
      group_by(j_1,j_2)%>%
      mutate(   psi_hat =  ifelse(j_1==j_2,func(W[,j_1],W[,j_2], X[,j_1]*X[,j_2]-1),
                                  func(W[,j_1],W[,j_2],X[,j_1]*X[,j_2])
      )  ) 
    return(sum(df$psi_hat))    
  }  
  ################################################################
  ###### semi-supervised data of X
  ##########################################################
  X_pop<-matrix(
    #rnorm(n_pop*p)
    1-rexp(n_pop *p,1)
    ,nrow = n_pop, ncol = p)
  conditional_mean_p <- as.matrix(X_pop)%*%beta_range
  prod (round(t(X_pop)%*%X_pop/n_pop) == diag(p))  # if prod == 1 then SIGMA = I
  mehane <- var(double_dist_sum(X_pop))
  
  ### optimal coefficients table
  df <- data.frame(j_1=rep(1:p,each=p), j_2=rep(1:p,times=p))%>%
    group_by(j_1,j_2)   %>%  dplyr::mutate(beta_jj  = beta_range[j_1]*beta_range[j_2])
  
  
  # sample data
  X_sample<-matrix(   
    #rnorm(n*p)   
    1-rexp(n *p,1) 
    ,nrow = n, ncol = p)
  conditional_mean <- as.matrix(X_sample)%*%beta_range 
  y_sample<- conditional_mean + rnorm(n,mean=0,sd=sqrt(sigma2))  # basic model
  ##############################################################################################333 
  ###                 naive estimator
  ##################################################################################3
  Y=y_sample
  X=X_sample
  W = as.data.frame(X)*Y
  mean_squared_W<- colMeans(W^2)  #calculate first element of beta_square_hat
  var_W<-sapply(W,var) # #calculate the second element of beta
  beta_square_hat<- mean_squared_W - var_W # calculate vector of beta2_hat
  naive<-sum(beta_square_hat)
  
  #########################################
  ###           Oracle estimator    
  OOE <- T_oracle_f(X)
  ##################################################3
  ###          Selection estimator (oracle)   
  T_Selection_estimator <- Selection_estimator(X,Y)
  #####################################################
  
  ###             Single coefficient estimator    
  X<- as.matrix(X)
  g <- mean(double_dist_sum(X))
  mone_s<-func22(X,Y)*4/(n*(n-1))
  (a_star_hat <- mone_s /mehane)
  single_coeff_estimator <-  naive - a_star_hat*g 
  ###################################################
  ###         PSI estimator   #####
  
  PSI <-   estimateSigma(X, Y, intercept=TRUE, standardize=TRUE)   # estimate estimateSigma from selectiveInference package
  ######################################################################
  
  
  #  create final table
  Naive <- naive
  OOE    <-  naive     -2 * OOE 
  Selection <- naive - 2 * T_Selection_estimator
  PSI <- var(as.vector(Y)) - PSI$sigmahat^2
  Single <- single_coeff_estimator
  
  return(tibble( Naive = naive
                 , OOE = OOE
                 , Selection = Selection
                 , PSI = PSI
                 , Single = Single
  ))
}


########################################
#### initial parameters settings
######################################
#availableCores()
n_cores =72  
set.seed(1)
n_simulation = 100
n_pop<-3000
tau2=1
sigma2=1
n=400
c = 5     #number of large betas

###########################################################################


setting_param_list <-  list(seed = 1:n_simulation, prop_large = c(0.33,0.66,0.99))
result <-run_sim(sim_fn, setting_param_list, nworkers = n_cores)
arrange(result, scenario)
raw_storage<-gather(result, method, value, 4:8) 
end_time <- Sys.time()
total_time  <- round(end_time - start_time)


###############################################
#    Boxplot 
##############################################
g_boxplot<-raw_storage %>%
  mutate(method=factor(method, levels=c(
    "Naive"
    ,"Single"
    ,"Selection"
    ,"OOE" 
    ,"PSI"
     ) )) %>%ggplot( aes(x=as.factor(prop_large), y=value
                      ,fill=method)) + geom_boxplot() + theme_bw() +
  ylab("Value")+ labs(fill = "Estimator") + xlab(TeX("$\\tau^2_{B}$"))+
  geom_hline(yintercept=tau2, linetype="dashed", color = "red")+
  scale_fill_brewer(palette="Pastel1") 




###############################################
#    Summary Table
##############################################

stat<- raw_storage %>% group_by(method,prop_large) %>%
  summarise(Mean =round(mean(value),2)
            ,SD = round(sd(value),3)
            ,MSE=mean((value-tau2)^2) 
            ,SD_RMSE = round(sqrt(var((value-tau2)^2)/(4*MSE*n_simulation)),3)*1000
  ) %>% mutate(RMSE = round(sqrt(MSE),3), bias = round(tau2-Mean,2), tau2=round(tau2,1), n=n
  )%>% arrange(desc(prop_large))


stat <- stat %>%   dplyr::select(prop_large,Estimator = method,Mean,bias
                                 ,SD
                                 ,RMSE
                                 , SD_RMSE
)%>%
  arrange(prop_large) %>%
  mutate(tau2=round(tau2,1), sigma2=round(sigma2,1),
         n_sim=n_simulation
         , "# large betas" =c
         ,prop_large =  label_percent()(prop_large),
         total_time=total_time, n=n, n_pop = n_pop )



g_boxplot
stat
