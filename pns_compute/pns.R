# read the momenta
library(shapes)
library(rsample)
library(sdpt3r)
library(geometry)
library(ROCR)

rocplot = function(pred, truth, ...){
  predob = prediction(pred, truth)
  perf = performance(predob, "tpr", "fpr")
  plot(perf,avg="threshold")}

get_header_info = function(filename){
  header <- readLines(filename, n=1)
  header_split <- strsplit(header, ' ')
  header_split <- as.numeric(header_split[[1]])
  
  num_subjs = as.integer(header_split[1])
  ctrl_pts = as.integer(header_split[2])
  return(c(num_subjs, ctrl_pts))
}

get_features = function(num_subjs, ctrl_pts, all_cases) {
  # for each subject
  f_mat <- matrix(, nrow=num_subjs, ncol=3*ctrl_pts)
  for (i in 1:num_subjs) {
    case_1 <- all_cases[((i-1)*ctrl_pts+2+i):(i*ctrl_pts+2+i-1)]
    case_1 <- do.call(rbind, case_1)
    case_1 <- matrix(as.numeric(case_1), ncol=ncol(case_1))
    
    # normalize to unit vector
    norms <- sqrt(rowSums(case_1^2))
    unit_vctrs <- case_1 / norms
    
    # Euclideanize the length of vectors by taking log
    geo_mean = mean(log(norms))
    f_1 = log(norms) - geo_mean
    
    # apply pns, to get new features f_2 and f_3, which are the residuals
    pns_res = pns(t(unit_vctrs), output=FALSE)
    
    f_2 = pns_res$resmat[1,]
    f_3 = pns_res$resmat[2,]
  
    f_mat[i,] <- c(f_1, f_2, f_3)
  }
  return(f_mat)
}


neg_momenta_fn = './momenta/newmean_momenta/neg_Momenta.txt'
pos_momenta_fn = './momenta/newmean_momenta/pos_Momenta.txt'

neg_header = get_header_info(neg_momenta_fn)
pos_header = get_header_info(pos_momenta_fn)

my_file_neg <- readLines(neg_momenta_fn)
my_file_pos <- readLines(pos_momenta_fn)
neg_cases <- strsplit(my_file_neg, ' ')
pos_cases <- strsplit(my_file_pos, ' ')

f_pos = get_features(pos_header[1], pos_header[2], pos_cases)
f_neg = get_features(neg_header[1], neg_header[2], neg_cases)

# pca on features obtained
pos_ctrl_pts = pos_header[2]
for (i in 1:3) {
  my_pca = prcomp(f_pos[,((i-1)*pos_ctrl_pts):(i*pos_ctrl_pts)], center = FALSE)
  if (i == 1){
    pos_f_mat = my_pca$x
  } else {
    pos_f_mat = cbind(pos_f_mat, my_pca$x)
  }
}

for (i in 1:3) {
  my_pca = prcomp(f_neg[,((i-1)*pos_ctrl_pts):(i*pos_ctrl_pts)], center = FALSE)
  if (i == 1){
    neg_f_mat = my_pca$x[,1:34]
  } else {
    neg_f_mat = cbind(neg_f_mat, my_pca$x[,1:34])
  }
}

f_pos = pos_f_mat
f_neg = neg_f_mat

# np_df = as.data.frame(rbind(f_pos, f_neg))
# pos_labels = rep(as.integer(1), pos_header[1])
# neg_labels = rep(as.integer(-1), neg_header[1])
# labels = c(pos_labels, neg_labels)
# np_df = cbind(np_df, labels)

#---------select------
smp_size <- pos_header[1]
set.seed(20)
select_ind <- sample(seq_len(neg_header[1]), size=smp_size)
f_neg_select = neg_f_mat[select_ind,]

np_df = as.data.frame(rbind(pos_f_mat, f_neg_select))
pos_labels = rep(as.integer(1), pos_header[1])
# so that have equal num of rows
neg_labels = rep(as.integer(-1), pos_header[1])
labels = c(pos_labels, neg_labels)
np_df = cbind(np_df, labels)


# np_df$labels = as.factor(np_df$labels)

# split training dataset 
# smp_size <- floor(0.8 * num_subjs)
# set.seed(123)
# train_ind <- sample(seq_len(num_subjs), size=smp_size)
set.seed(10)
k_folds_pos = vfold_cv(np_df[np_df$labels==1,], v=5)
set.seed(5)
k_folds_neg = vfold_cv(np_df[np_df$labels==-1,], v=5)


k_scores = list()
k_labels = list()
for (i in 1:1:nrow(k_folds_pos)) {
  train_set_pos = analysis(k_folds_pos$splits[[i]])
  train_set_neg = analysis(k_folds_neg$splits[[i]])
  
  val_set_pos = assessment(k_folds_pos$splits[[i]])
  val_set_neg = assessment(k_folds_neg$splits[[i]])
  
  val_set = rbind(val_set_pos, val_set_neg)
  
  # train_set_pos = as.matrix(train_set[(train_set$labels == 1),])
  # train_set_pos = train_set_pos[1:25,]
  # train_set_neg = as.matrix(train_set[(train_set$labels == -1),])
  
  train_set_pos = subset(train_set_pos, select=-labels)
  train_set_neg = subset(train_set_neg, select=-labels)
  
  new_val_set = subset(val_set, select=-labels)
  
  out = dwd(as.matrix(train_set_pos), as.matrix(train_set_neg), 1)
  d = ncol(train_set_neg)
  omega = out$X[[1]][2:(d+1)]
  beta =  out$X[[1]][d+3]
  
  norm_omega = sqrt(sum(omega^2))
  scores = vector("numeric", nrow(new_val_set))
  for (i in 1:(nrow(new_val_set))){
    scores[i] = (dot(omega, as.numeric(new_val_set[i,])) + beta) / norm_omega
  }
  k_scores = append(k_scores, list(scores))
  k_labels = append(k_labels, list(val_set$labels))
  print(scores)
  print(val_set$labels)
  # set.seed(1)
  # tune_out = tune(svm,
  #                 labels~.,
  #                 data = train_set,
  #                 kernel = "linear",
  #                 ranges = list(cost = c(0.001, 0.01, 0.1, 1,5,10,100)))
  # print(summary(tune_out))
  # bestmod = tune_out$best.model
  # print(summary(bestmod))
  # class_pred = predict(bestmod, val_set)
  # 
  # svmfit = svm(labels~., data=train_set, kernel='linear',cost=100, scale = FALSE)
  # class_pred = predict(svmfit, val_set)
  # 
  # aa = table(predicted = class_pred, true = val_set$labels)
  # 
  # fitted_opt_train = attributes(predict(svmfit, val_set,decision.values=TRUE))$decision.values
}
rocplot(k_scores, k_labels, main = "Training Data")
