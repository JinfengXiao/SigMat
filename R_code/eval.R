# This script evaluates the Success@K on the testing data.
# Run under SigMat/R_code

rm(list = ls())
test_drug = unlist(read.table("../data/test/class_test.csv", colClasses = "character"))
pred_drug = as.matrix(read.table("../data/test/class_pred.csv", sep = ",", colClasses = "character"))
if(nrow(pred_drug) == 1 ){
  if(length(test_drug) == 1){ # Only 1 test sample
    success = (test_drug[1] %in% pred_drug[1, ])
  }
  else{ # topK = 1
    pred_drug = pred_drug[1, ]
    success = sum(test_drug == pred_drug) / length(test_drug)
  }
}else{ # topK > 1
  success = sum(apply(cbind(test_drug, pred_drug), 1, function(x){
    x[1] %in% x[2:length(x)]
  })) / length(test_drug)
}

print(success)
