library(MASS)
library(pscl)
library(grpreg)
library(mpath)
library(zic)
library(SGL)


gooogle<-function(data,xvars,zvars,yvar,group=1:ncol(data), 
                  samegrp.overlap=T,
                  penalty=c("grLasso", "grMCP", "grSCAD","gel", "cMCP",
                             "gBridge", "gLasso", "gMCP"),
                  dist=c("poisson","negbin"),
                  nlambda=ifelse(penalty=="SGL",20,100), lambda, 
                  family="gaussian", 
                  lambda.min=ifelse((nrow(data[,unique(c(xvars,zvars))])>ncol(data[,unique(c(xvars,zvars))])),1e-4,.05), lambda.max,
                  crit="BIC",  alpha=ifelse(penalty=="SGL",0.95,1), eps=.001,
                  max.iter=1000, gmax=length(unique(group)), 
                  gamma=ifelse(penalty=="SGL",0.8,ifelse(penalty=="gBridge",0.5,ifelse(penalty == "grSCAD", 4, 3))),tau = 1/3, warn=TRUE, delta=1e-7, 
                  nfolds=10, returnY=FALSE,trace=FALSE, thresh = 0.001, 
                  min.frac = 0.1, standardize = T, verbose = FALSE,step = 1,
                  reset = 10, ...)
{
  y<-data.frame(data[,yvar])
  X<-data.frame(data[,xvars])
  pred.names<-setdiff(names(data),yvar)
  names(X)<-paste(names(X),"_count",sep="")
  names(y)<-yvar
  group.x<-group[which(pred.names %in% xvars)]
  
  if (is.null(zvars))
  {
    group.z<-NULL
    data<-cbind.data.frame(y,X)
    xvars<-names(X)
    fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|1"),sep=""))
  } else {
    Z<-data.frame(data[,zvars])
    names(Z)<-paste(names(Z),"_zero",sep="")
    
    if (samegrp.overlap)
    {
      group.z<-group[which(pred.names %in% zvars)]
    } else {
      group.z<-max(group.x)+group[which(pred.names %in% zvars)]
    }
    data<-cbind.data.frame(y,X,Z)
    xvars<-names(X)
    zvars<-names(Z)
    fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|",paste(zvars,collapse="+")),sep=""))
  }

  fit.zero <- zeroinfl(fit.formula, dist=dist,data = data)

#   yvar<-names(data)[1]
#   xvars<-names(data)[2:(ncol(X)+1)]
#   zvars<-names(data)[-(1:(ncol(X)+1))]
#   p<-dim(X)[2]
#   q<-dim(Z)[2]
#   browser() 
  
  if (dist=="poisson")
  {
    coeff<-c(fit.zero$coefficients$count,fit.zero$coefficients$zero)
    vcov<-fit.zero$vcov
    e<-eigen(vcov)

    if(det(vcov)>0)
    {
      cov.star=e$vectors%*%diag(1/sqrt(abs(e$values)))%*%t(e$vectors)
    } else {
      cov.star=e$vectors%*%diag(1/sqrt(PDFORCE(vcov)))%*%t(e$vectors)
    }
    
    y.star<-cov.star%*%coeff
    cov.star<-data.frame(scale(cov.star))
    y.star<-y.star-mean(y.star)
    names(cov.star)<-c("int_count",xvars,"int_zero",zvars)
    
    group<-c(0,group.x,0,group.z)
    u<-unique(group)
    cov.star_reordered<-rep(0,dim(cov.star)[1])
    for(i in 1:length(u))
    {
      index<-which(group==u[i])
      cov.star_reordered<-cbind(cov.star_reordered,cov.star[,index])
    }
    cov.star_reordered<-cov.star_reordered[,-1]
    group<-sort(group)
    
    # browser()
    
    if (penalty=="gBridge")
    {
      fit<-gBridge(X=cov.star_reordered, y=y.star, group=group, nlambda=nlambda, lambda, lambda.min = lambda.min, lambda.max, delta=delta, max.iter=max.iter, gamma=gamma, family=family, alpha=alpha, eps=eps,  warn=warn)

      fit.final=grpreg::select(fit,crit=crit)
#       fit.final$beta<-fit.final$beta[names(cov.star)[-length(names(cov.star))]]
    } else if (penalty=="SGL") {
      data.sgl<-list(x=cov.star_reordered,y=y.star)
      if(missing(lambda)) lambdas<-NULL else lambdas<-lambda
      fit.final<-SGL(data.sgl, index=group, type = "linear", maxit = max.iter, thresh = thresh, min.frac = min.frac, nlam = nlambda, gamma = gamma, standardize = standardize, verbose = verbose, step = step, reset = reset, alpha = alpha, lambdas)
    }
    else {
      fit=grpreg(cov.star_reordered,y.star,group=group,penalty=penalty,family="gaussian", nlambda=nlambda, lambda, lambda.min=lambda.min,  alpha=alpha, eps=eps, max.iter=max.iter,gmax=gmax, gamma=gamma, tau = tau,  warn=warn)
      fit.final=grpreg::select(fit,crit=crit)
#       fit.final$beta<-fit.final$beta[names(cov.star)[-length(names(cov.star))]]
    }
  } else if(dist=="negbin")
  {
    a<-1/fit.zero$theta
    coeff<-c(a,fit.zero$coefficients$count,fit.zero$coefficients$zero)
    var.a<-(fit.zero$SE.logtheta)^2*exp(2*log(a))
    vcov<-fit.zero$vcov
    vcov<-cbind(c(var.a,rep(0,length(coeff)-1)),rbind(rep(0,length(coeff)-1),vcov))

    e<-eigen(vcov)
    
    if(round(det(vcov),3)>0)
    {
      cov.star=e$vectors%*%diag(1/sqrt(abs(e$values)))%*%t(e$vectors)
    } else {
      cov.star=e$vectors%*%diag(1/sqrt(PDFORCE(vcov)))%*%t(e$vectors)
    }
    # browser()
    y.star<-cov.star%*%coeff
    cov.star<-data.frame(scale(cov.star))
    y.star<-y.star-mean(y.star)
    names(cov.star)<-c("dispersion","int_count",xvars,"int_zero",zvars)
    
    group<-c(0,0,group.x,0,group.z)
    u<-unique(group)
    cov.star_reordered<-rep(0,dim(cov.star)[1])
    for(i in 1:length(u))
        {
          index<-which(group==u[i])
          cov.star_reordered<-cbind(cov.star_reordered,cov.star[,index])
    }
    cov.star_reordered<-cov.star_reordered[,-1]
    group<-sort(group)
    

    if (penalty=="gBridge")
    {
      fit<-gBridge(X=cov.star_reordered, y=y.star, group=group, nlambda=nlambda, lambda, lambda.min = lambda.min, lambda.max, delta=delta, max.iter=max.iter, gamma=gamma, family=family, alpha=alpha, eps=eps,  warn=warn)
      
      fit.final=grpreg::select(fit,crit=crit)
#       fit.final$beta<-fit.final$beta[names(cov.star)[-length(names(cov.star))]]
    } else if (penalty=="SGL") {
      data.sgl<-list(x=cov.star_reordered,y=y.star)
      if(missing(lambda)) lambdas<-NULL else lambdas<-lambda
      fit.final<-SGL(data.sgl, index=group, type = "linear", maxit = max.iter, thresh = thresh, min.frac = min.frac, nlam = nlambda, gamma = gamma, standardize = standardize, verbose = verbose, step = step, reset = reset, alpha = alpha, lambdas)
    }
    else {
      fit=grpreg(cov.star_reordered,y.star,group=group,penalty=penalty,family="gaussian", nlambda=nlambda, lambda, lambda.min=lambda.min,  alpha=alpha, eps=eps, max.iter=max.iter, gmax=gmax, gamma=gamma, tau = tau,  warn=warn)
      fit.final=grpreg::select(fit,crit=crit)
#       fit.final$beta<-fit.final$beta[names(cov.star)[-length(names(cov.star))]]
    }
  } else {
    stop("Gooogle only works with Poisson or Negative Binomial distribution for count data")
  }
  return(fit.final)
}



PDFORCE = function(Q){
  N = nrow(Q)
  HC = Q
  D = eigen(Q)
  E = D$values
  U = D$vectors
  v = as.numeric(E <= 0)
  m = sum(v) # number of negative values
  if(m > 0){
    S = sum(v*E)*2
    W = (S*S*100)+1
    P = E[N - m] # smallest positive value
    k = N - m + 1
    for(i in k:N){
      C = E[i]
      E[i] = P * (S-C)*(S-C)/W
    }
#     HC = U %*% diag(E) %*% t(U)
  }
  return(E) }
