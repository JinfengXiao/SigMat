suppressMessages(library(kernlab))

# Calculate the mean prediction accuracy on all classes
mean_class_accuracy = function(truth, pred){
  class_v = sort(unique(truth))
  class_ac = rep(0, length(class_v))
  for(k in 1:length(class_ac)){
    class_ac[k] = sum((pred == truth) & (truth == class_v[k])) / sum(truth == class_v[k])
  }
  return(mean(class_ac))
}

# Train a kernel SVM ensemble
kSVM_train = function(class, kernel){
  cat("Training a KSVM ensemble...\n")
  kernel = as.kernelMatrix(kernel)
  model = ksvm(x = kernel, y = as.factor(class), kernel = "matrix")
  return(model)
}

# Given a kernel matrix, return an n by K(K-1)/2 matrix, storing the predicted f functions returned by the binary classifiers in kSVM. n: Number of observations. K: Number of classes.
kSVM_test_get_f = function(model, kernel){
  n_obs = nrow(kernel)
  K = model@nclass
  f = matrix(0, nrow = n_obs, ncol = K * (K - 1) / 2)
  for(obs_i in 1:n_obs){
    for(bc_i in 1:ncol(f)){  # bc_i: Iterator over binary classifiers
      f[obs_i, bc_i] = sum(model@coef[[bc_i]] * kernel[obs_i, model@alphaindex[[bc_i]]]) - model@b[bc_i]
    }
  }
  return(f)
}

# Return the vote matrix using trained kernel SVM. Rows are observations. Columns are classes.
kSVM_vote_m = function(model, kernel){
  K = model@nclass
  
  first_class = rep(0, K * (K - 1) / 2)
  second_class = rep(0, K * (K - 1) / 2)
  tracker = 1
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      first_class[tracker] = i
      second_class[tracker] = j
      tracker = tracker + 1
    }
  }
  
  f = kSVM_test_get_f(model, kernel)
  
  n_obs = nrow(kernel)
  predicted = rep(0, n_obs)
  v0 = rep(0, K)
  
  score_m = sapply(1:n_obs, function(x) {
    v = rep(0, K)
    vb = rep(0, K * (K - 1) / 2)  # Votes from binary classifiers
    index = (f[x,] <= 0)
    vb[index] = first_class[index]
    index = (f[x,] > 0)
    vb[index] = second_class[index]
    
    for(i in 1:K){
      v[i] = sum(vb == i)
    }
    return(v)
  })
  score_m = t(score_m)
  
  return(score_m)
}

# Rescale a numeric data array to [a, b]
num_rescale = function(dat, a, b){
  dat_res = (b - a) * ((dat - min(dat)) / (max(dat) - min(dat))) + a
  return(dat_res)
}

# Build a cor class score matrix
cor_score_m = function(kernel, class_train){
  score_m = matrix(0, nrow = nrow(kernel), ncol = max(class_train))
  for(k in 1:ncol(score_m)){
    score_m[, k] = apply(kernel[, class_train == k], 1, max)
  }
  return(score_m)
}

# Predict with a trained kernel SVM ensemble
kSVM_test = function(model, kernel){
  v = kSVM_vote_m(model, kernel)
  Y_pred = apply(v, 1, which.max)
  return(Y_pred)
}

# Find the best kernel rescaling parameter alpha by grid search
kSVM_scale = function(model, class, kernel){
  cat("Searching for the best alpha value between 0.2 and 7. This can take long...\n")
  alpha_grid = c(seq(0.2, 1, 0.2), seq(1.5, 7, 0.5))
  max_alpha = 0
  max_ac = 0
  for(alpha in alpha_grid){
    cat(paste0("  Trying alpha = ", alpha, " ...\n"))
    model_alpha = model
    model_alpha@coef = lapply(model@coef, function(x) x * alpha)
    pred = kSVM_test(model_alpha, kernel)
    ac_alpha = mean_class_accuracy(class, pred)
    
    if(ac_alpha > max_ac){
      max_alpha = alpha
      max_ac = ac_alpha
    }
  }
  return(max_alpha)
}

# Given a predicted observation-by-class score matrix and the true label, calculate Success@K.
success_at_K = function(score_m, class, K){
  if(K > 1){
    pred = t(apply(score_m, 1, function(x){
      return(order(x, decreasing = T)[1:K])
    }))
    success = apply(pred, 2, function(x){
      return(x == class)
    })
    success = apply(success, 1, sum)
  }else{
    success = apply(score_m, 1, which.max)
    success = (success == class)
  }
  success = sum(success) / length(success)
  return(success)
}
