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


rolling_av <- function(eps2,cut){
  
  n = length(eps2)
  n_cut = floor(n*cut)
  res = rollapply(eps2, width = 10, FUN = mean, align = "right", partial = TRUE)
  e = length(res)-1
  return(res[n_cut:e])}
