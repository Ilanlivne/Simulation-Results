rm(list = ls())
packages <- c("ggplot2","tidyr","dplyr","latex2exp","gridExtra","pushoverr","Rcpp","scales","stargazer","selectiveInference" )
library("selectiveInference")
lapply(packages, require, character.only = TRUE)
source("run_sim.R")
Rcpp::sourceCpp("zero_terms_functions.cpp")


start_time <- Sys.time()


sim_fn <- function(setting_param_list){
Rcpp::sourceCpp("zero_terms_functions.cpp")

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
  ##############################################################################################333 
  ###                 naive estimator
  ##################################################################################3
  Y=y_sample
  X=X_sample
  
  PSI <-   estimateSigma(X, Y, intercept=TRUE, standardize=TRUE)   
  PSI <- var(as.vector(Y)) - PSI$sigmahat^2

  X<- as.matrix(X)
  g_int <- mean(double_dist_sum(X))
  ###   Bootstraping step 
  PSI_Pri_b = g_int_b  = NA
  for (b in 1:B) {
    # step 3.1: Resample with replacement
    indices <- sample(1:nrow(X),size = nrow(X) ,replace=TRUE)
    X_b <- X[indices,]
    Y_b <- Y[indices]
    # step 3.2: calculate initials estimator:  
    
    PSI_b <-   estimateSigma(X_b, Y_b, intercept=TRUE, standardize=TRUE)   
    PSI_b <- var(as.vector(Y_b)) - PSI_b$sigmahat^2
    
    PSI_Pri_b[b] <-PSI_b
    g_int_b[b] <- mean(double_dist_sum(X_b))
  }# end of bootstrapping loop
  
  cov_emp =cov(PSI_Pri_b, g_int_b)
  c_star_mean <- n*cov_emp/mehane
  
  
  return(tibble( 
       PSI = PSI
      ,"Empirical PSI" = PSI - c_star_mean*g_int
    
  ))
}

######################################################################################

########################################
#### initial parameters settings
######################################
#availableCores()
n_cores =72  
#set.seed(1)
n_simulation = 100
B=100
n_pop<-3000
tau2=1
sigma2=1
c = 5     #number of large betas
#######################################################################################

setting_param_list <-  list(seed = 1:(n_simulation), prop_large = c(0.33,0.66,0.99),n=400)

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
    "PSI"
    ,"Empirical PSI"
    
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
