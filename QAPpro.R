# ========= Initiation ===========
library(dplyr)
library(ggplot2)
library(car)
if (!require(clubSandwich)){
  install.packages("clubSandwich")
  library(clubSandwich)
}


# ========= Super Population Tests ============

#' DyadicLM: a function that runs linear regression on dyadic data
#'
#' @param form 
#' @param data_list 
#' @param var_scheme 
#'
#' @return
#' @export
#'
#' @examples
#' 
dyadicLM = function(form, data_list, var_scheme = "CR0"){
  ## dimension of regressors
  p = length(data_list) - 1
  n = nrow(data_list[[1]])
  
  ## from list of matrices to a vectorized data frame
  data = data.frame(
    sapply(
      1:length(data_list),
      function(x){
        c(data_list[[x]])
      }
    ))
  names(data) = names(data_list)
  
  ## run regression and obtain coefficients - must go with an intercept
  dyadic_fit = lm(formula = as.formula(form), data = data)
  
  ## obtain variance estimation
  ## require package: clubSandwich 
  ## several types of var_schemes: "CR0", "CR1", "CR1p", "CR1S", "CR2", or "CR3"
  cluster0 = factor(rep(1:n, each = n))
  var_mat = vcovCR(
    dyadic_fit,
    cluster = cluster0,
    type = var_scheme
  )
  
  list(
    coefs = dyadic_fit$coefficients,
    var_mat = 4 * var_mat,
    dyadic_fit = dyadic_fit
  )
  
}


## ============== Permutation tests ================

QAPpro = function(form, 
                  data_list, 
                  mode = "permute-Y", 
                  num_perms = 1000, 
                  var_scheme = "CR0", 
                  plot.flag = F, 
                  nameX.target = all.vars(formula(form))[-1],
                  user.seed = Sys.time()){
  
  ## dimension of regressors
  p = length(data_list) - 1
  n = nrow(data_list[[1]])
  nameY = all.vars(formula(form))[1]
  nameX = all.vars(formula(form))[-1]
  if (any(nameX == ".")){
    nameX = setdiff(names(data_list), nameY)
    nameX.target = nameX
  }

  
  ## compute observed value
  res = dyadicLM(form, data_list, var_scheme = var_scheme)
  coefs = res$coefs
  var_mat = res$var_mat
  pred_mat_part = Reduce(`+`, Map(`*`, coefs[setdiff(nameX, nameX.target)], data_list[setdiff(nameX, nameX.target)])) + coefs[1]
  pred_mat_full = Reduce(`+`, Map(`*`, coefs[nameX], data_list[nameX])) + coefs[1]
  res_mat = data_list[[nameY]] - pred_mat_full

  
  stat_NS = sqrt(n) * coefs[nameX.target]
  stat_S = coefs[nameX.target]/sqrt(diag(var_mat)[nameX.target])
  stat_wald = t(coefs[nameX.target]) %*% solve(var_mat[nameX.target, nameX.target]) %*% coefs[nameX.target]
  
  ## next steps
  record_perm_NS = matrix(0, nrow = num_perms, ncol = length(nameX.target))
  record_perm_S = matrix(0, nrow = num_perms, ncol = length(nameX.target))
  record_wald = rep(0, num_perms)
  
  ## set up a progress bar
  pb = txtProgressBar(min = 0, max = num_perms, style = 3) # set progress bar
  
  ## permutation test
  set.seed(user.seed)
  for (iter in 1:num_perms){
    permInd = sample(1:n)
    data_list_perm = data_list
    
    ### permute the data list
    if (mode == "permute-Y"){
      data_list_perm[[nameY]] = data_list_perm[[nameY]][permInd, permInd]
    }
    else if (mode == "permute-X"){
      for (idn in nameX.target){
        data_list_perm[[idn]] = data_list_perm[[idn]][permInd, permInd]
      }
    }
    else if (mode == "permute-e"){
      res_mat_perm = res_mat[permInd, permInd]
      Y_perm = pred_mat_part + res_mat_perm
      data_list_perm[[nameY]] = Y_perm
    }
    else if (mode == "permute-eX"){
      
    }
    else{
      stop("No such permutation method defined!")
    }
    
    ### regression results
    res = dyadicLM(form = form, data_list = data_list_perm, var_scheme = var_scheme)
    coefs = res$coefs
    var_mat = res$var_mat
    
    record_perm_NS[iter, ] = sqrt(n) * coefs[nameX.target]
    record_perm_S[iter, ] = coefs[nameX.target]/sqrt(diag(var_mat)[nameX.target])
    record_wald[iter] = t(coefs[nameX.target]) %*% solve(var_mat[nameX.target, nameX.target]) %*% coefs[nameX.target]
    
    
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  
  
  ## report p values
  pval_NS = sapply(1:length(nameX.target), function(x){(sum(record_perm_NS[,x] > abs(stat_NS[x])) + sum(record_perm_NS[,x] < -abs(stat_NS[x])))/num_perms})
  pval_S = sapply(1:length(nameX.target), function(x){(sum(record_perm_S[,x] > abs(stat_S[x])) + sum(record_perm_S[,x] < -abs(stat_S[x])))/num_perms})
  pval_W = sum(record_wald > stat_wald[1,1])/num_perms
  
  report = data.frame(
#    X.name = nameX.target,
    coefs = stat_NS/sqrt(n),
    t.stat = stat_S,
    p.val = format(2 * pnorm(abs(stat_S), lower.tail = F), scientific = FALSE, nsmall = 4),
    p.val.coefs.perm = pval_NS,
    p.val.t.perm = pval_S
  )
  
  ## return the p values
  list(
    record_perm_NS = record_perm_NS,
    record_perm_S = record_perm_S,
    report = report, 
    wald = list(wald.stat = stat_wald[1,1], p.val.wald = pval_W)
  )

}

















