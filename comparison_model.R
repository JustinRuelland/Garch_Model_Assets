test_mariano <- function(pred1,pred2,val,hor){
  e1 = pred1-val
  e2 = pred2-val
  return(dm.test(pred1-val, pred2-val, alternative = c("two.sided"), h = hor, power = 2))}


pred_h1_garch <- function(eps2,cut){
  
  n = length(eps2)
  n_cut = floor(n*cut)
  
  theta =  QML(eps2[1:n_cut])
  pred = double(n-n_cut+1)
  init = theta[1]/(1-theta[2]-theta[3])
  pred[1] = func_sigma2(n_cut,init,eps2[1:n_cut],theta)
  
  l = n-n_cut+1
  for(i in 2:l){
    pred[i] = theta[1]+theta[2]*eps2[n_cut+i-2]+theta[3]*pred[i-1]}
  return(pred[2:l+1])}


rolling_av <- function(eps2,cut,wind,liss_exp){
  
  n = length(eps2)
  n_cut = floor(n*cut)
  if(liss_exp==TRUE){res = movavg(eps2,n = wind, type='e')}
  else{res = rollapply(eps2, width = wind, FUN = mean, align = "right", partial = TRUE)}
  e = length(res)-1
  return(res[n_cut:e])}




HMM_Pred <- function(data_Train) {
  a <- RHmm::HMMFit(data_Train, nStates=5, nMixt=4, dis="MIXTURE", control=list(iter=2000))
  v <- RHmm::viterbi(a,data_Train)
  return(sum(a$HMM$transMat[last(v$states),] * .colSums((matrix(unlist(a$HMM$distribution$mean), nrow=4,ncol=5)) * (matrix(unlist(a$HMM$distribution$proportion), nrow=4,ncol=5)), m=4,n=5)))
}


HMM_Pred_finale <- function(eps2, cut) {
  n = length(eps2)
  n_cut = floor(n*cut)
  df <- data.frame(Col1 = double())
  for(i in n_cut:n){
    df[i+1,1] <- HMM_Pred(eps2[1:i])
  }
  return(df)
}
