library(future)
library(lubridate)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)




#' An example for a sim_fn function
#' @param setting_param_list parameters that changes at each run
#' @param other_param_list static parameters for all the runs
#' @param nworkers number of cpus to use
#' @return a tibble with the results and also save an RDS of the results



run_sim <- function(sim_fn,setting_param_list,nworkers = 70){
  
  
  future::plan("multiprocess",workers =  min(nworkers,availableCores()-1))
  # the function safely used to avoid program failing.
  sim_fn <- purrr::safely(sim_fn)
  
  
  results <- do.call(tidyr::expand_grid, setting_param_list) %>%
    dplyr::mutate(scenario = 1:n())
  
  results <- results %>%
    # put the different parameters for each simulation in a tibble
    dplyr::group_by(scenario) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    ##Randomise the order of scenarios - helps share the load across cores
    dplyr::sample_frac(size = 1, replace = FALSE) %>%
    dplyr::mutate(res = 
                    furrr::future_map(data,
                                      ~sim_fn(.)
                                      ,      .progress = TRUE
                    )) %>%
    tidyr::unnest(cols = "data") 
  
  results <- dplyr::mutate(results,res = purrr::map(res,function(resi){return(resi$result)})) %>% 
    tidyr::unnest(cols = "res") 
  
  #name = paste0("res",now() %>% stringr::str_sub(1,16),".RDS")
  
  
 # saveRDS(object = results,file = name,version = 2)
  
  future::plan("sequential",.cleanup = TRUE)
  
  return(results)
  
}
