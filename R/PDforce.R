## function for converting a matrix to pd if it's not ##
PDFORCE = function(mat){
  # browser()
  N = nrow(mat)
  HC = mat
  D = eigen(mat)
  E = D$values
  U = D$vectors
  v = as.numeric(E < 0)
  idx = which(v==1)
  m = sum(v) # number of negative values
  if(m > 0){
    S = sum(v*E)*2
    W = (S*S*100)+1
    P = min(abs(E[idx])) # smallest positive value
    for(i in idx){
      C = E[i]
      E[i] = P * (S-C)*(S-C)/W
    }
    #     HC = U %*% diag(E) %*% t(U)
  }
  return(E) }
