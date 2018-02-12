# Run under SigMat/R_code

# Version info & Welcome message
cat("\n========================================\n")
cat("Welcome to SigMat v1.0.1 (Updated Feb. 12, 2018).\n")
cat("SigMat matches gene expression signatures to experimental conditions, i.e. classes.\n")
cat("Documentation is available at https://github.com/JinfengXiao/SigMat.\n")
cat("Command for running SigMat: Rscript sigmat_main.R [-k int]\n")
cat("k is a small positive integer indicating how many classes you want a signature to be matched to.\n")
cat("Default: k=10.\n")
cat("========================================\n\n")

# Read in topK from command line arguments
library("optparse")
option_list = list(
  make_option(c("-k", "--topk"), type="numeric", default=10, 
              help="Number of top predicted classes to return", metavar="topK")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
topK = opt$topk
cat(paste0("You have instructed SigMat to match signatures to the top ", topK, " class(es).\n"))

source("functions.R")

# Read in signatures
cat("Reading data...\n")
X_train = as.matrix(read.table("../data/train/sig_train.csv", header = F, sep = ",", colClasses = "numeric"))
X_tune = as.matrix(read.table("../data/tune/sig_tune.csv", header = F, sep = ",", colClasses = "numeric"))
X_test = as.matrix(read.table("../data/test/sig_test.csv", header = F, sep = ",", colClasses = "numeric"))

# Read in drug IDs of signatures
drug_train = unlist(read.table("../data/train/class_train.csv", header = F, colClasses = "character"))
drug_tune = unlist(read.table("../data/tune/class_tune.csv", header = F, colClasses = "character"))

# Convert drug IDs to classes
drug_uniq = unique(drug_train)
class_conv = seq(1, length(drug_uniq), 1)  # A converter between drug IDs and class labels
names(class_conv) = drug_uniq
class_train = class_conv[drug_train]
class_tune = class_conv[drug_tune]

# Compute kernels
kernel_train = exp(cor(t(X_train), method = "spearman"))
corr_tune = cor(t(X_tune), t(X_train), method = "spearman")
kernel_tune = exp(corr_tune)
corr_test = cor(t(X_test), t(X_train), method = "spearman")
kernel_test = exp(corr_test)

# Train SigMat and let it vote
sigmat = kSVM_train(class_train, kernel_train)
alpha = kSVM_scale(sigmat, class_tune, kernel_tune)
cat(paste0("SigMat found best alpha = ", alpha, ".\n"))
sigmat@coef = lapply(sigmat@coef, function(x) x * alpha)
cat("Collecting votes...\n")
vote_tune = kSVM_vote_m(sigmat, kernel_tune)
vote_test = kSVM_vote_m(sigmat, kernel_test)

# Tune beta
cat("Searching for the best beta value...\n")
ac_tune = rep(0, 11)
for(beta_i in 1:11){
  beta = 0.1 * (beta_i - 1)
  score_tune = beta * num_rescale(vote_tune, -1, 1) + (1 - beta) * cor_score_m(corr_tune, class_train)
  ac_tune[beta_i] = success_at_K(score_tune, class_tune, K=topK)
}
beta_best = 0.1 * (which.max(ac_tune) - 1)
cat(paste0("SigMat found best beta = ", beta_best, ".\n"))

# Rank classes on test data and output topK classes
cat("Matching test signatures to classes...\n")
score_test = beta_best * num_rescale(vote_test, -1, 1) + (1 - beta_best) * cor_score_m(corr_test, class_train)
pred_drug = t(apply(score_test, 1, function(x){
  return(names(class_conv)[order(x, decreasing = T)[1:topK]])
}))
write.table(pred_drug, file = "../data/test/class_pred.csv", sep = ",", quote = F, row.names = F, col.names = F)
cat("All done!\n")
