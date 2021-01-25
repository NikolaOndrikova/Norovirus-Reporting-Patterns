# source('./src/required_packages.R')

expandRegionalPreditor <- function(my_vector, my_Matrix, n){
  # expand regional vectors to the shape of hhh4 data matrix
  # n = # of age groups
  
  my_df = data.frame(t(my_vector))
  new_matrix = data.matrix(my_df[rep(seq_len(nrow(my_df)), each=nrow(my_Matrix)),])
  my_result = matrix(rep(t(new_matrix), n), ncol=ncol(my_Matrix))
  return(my_result)
}


makehhh4Control <- function(end_option, epi_option){
  # make controls for hhh4 models from the study
  
  control <- list(optimizer = list(stop = list(niter=500) #tol=1e-5,
  ),
  end = list( f = end_option, offset = shared_data$pop),
  ne = list( f = epi_option,
             weights = W_powerlaw(maxlag = 4, log = TRUE,
                                  normalize = FALSE, from0 = TRUE),
             scale = expandC(uk_con_norm, NREGIONS),
             normalize = TRUE),
  family = "NegBin1",
  data = c(group_indicators, shared_data, regional_predictors),
  subset = TRAIN,
  keep.terms = TRUE)
  return(control)   
}


updatehhh4 <- function(basic_model, end_option, epi_option, SEED = 1245){
  # update basic_model with different end_option and epi_option
  
  set.seed(SEED)
  ri_model <- update(basic_model,
                     family = "NegBin1",
                     optimizer = list(stop = list(tol=1e-5, niter=250),
                                      regression = list(method="nlminb"),
                                      variance = list(method="Nelder-Mead")),
                     end = list(f = end_option),
                     ne = list(f = epi_option),
                     subset = TRAIN,
                     keep.terms = TRUE
  )
  
  return(ri_model)   
  
}

selectBestModel <- function(owas, metrics = c("rps", "logs"), verbose = TRUE, plot = FALSE){
  # returns names of the best models per metric
  pvals = c()
  owas_mean_scores = c()
  perm_test = c()
  best_models = c()
  for (i in seq(1,length(metrics))){
    pvals[[i]] = unlist(sapply(owas, calibrationTest, which = metrics[i])["p.value",])
    calibrated_models = names(which(pvals[[i]]  > 0.050))
    
    owas_scores = lapply(owas[calibrated_models], scores, which = metrics, individual = TRUE, reverse = FALSE)
    owas_mean_scores = t(sapply(owas_scores, colMeans, dims = 2))
    ordered_error_index = order(owas_mean_scores[,i])
    models_to_compare = row.names(owas_mean_scores[ordered_error_index[1:2],])
    perm_test[[i]] = permutationTest(colMeans(owas_scores[[models_to_compare[1]]])[,i],
                                     colMeans(owas_scores[[models_to_compare[2]]])[,i], plot = plot)
    best_models[[metrics[i]]] = models_to_compare[1]
    
  }
  
  if (verbose == TRUE){
    options(digits = 3)
    cat("------------------------------------","\n")
    cat("Mean scores:","\n")
    print(owas_mean_scores)
    for (i in seq(1,length(metrics))){
      cat("////////////////////////////////////","\n")
      cat("------------------------------------","\n")
      cat(paste0(metrics[i], ": p-values"),"\n")
      print(data.frame(t(pvals[[i]])))
      cat("------------------------------------","\n")
      cat(paste0(metrics[i], ": permutation test"),"\n")
      print(data.frame(t(unlist(perm_test[[i]]))))
      cat("------------------------------------","\n")
      cat(paste0(metrics[i], ": best performing model (name)"),"\n")
      cat(best_models[[i]],"\n")
    }
  }
  
  return(best_models)
}
