calcAupr <- function(pred, gs) {
   ord.idx <- order(abs(pred), decreasing = T)
   
  prec <- cumsum(gs[ord.idx]) / cumsum(rep(1, length(ord.idx))) #also known as positive predictive value
  rec  <- cumsum(gs[ord.idx]) / sum(gs)                     #also know as true positive rate
  fpr  <- cumsum(gs[ord.idx] == 0) / (length(gs) - sum(gs)) #false positive rate

  prec <- c(prec[1], prec)
  rec <- c(0, rec)
  fpr <- c(0, fpr)
  
  aupr <- areaUnderCurve(rec, prec)
  auroc <- areaUnderCurve(fpr, rec)
  
  return(list(prec=prec, rec=rec, fpr = fpr, AUPR = aupr, AUROC = auroc))
}

areaUnderCurve <- function(x, y) {
 dx <- diff(x)
 my <- y[1:(length(y) - 1)] + diff(y) / 2
 return(sum(dx * my))
}