#' Group Regularization for zero Inflated Count Regression Models
#'
#' @param data This is the data frame containing outcome and predictors
#' @param xvars This is the vector of variable names for predictors in count model
#' @param zvars This is the vector of variable names for predictors in zero model
#' @param yvar This is the outcome variable name
#' @param group This is a vector of integers describing the grouping of the coefficients. For greatest efficiency and least ambiguity, it is best if group is a vector of consecutive integers. If any coefficients are to be included in the model without being penalized, their grouping index should be zero
#' @param samegrp.overlap This is a logical argument. If TRUE same grouping indices will be assigned to shared predictors in count and zero model
#' @param penalty This is the penalty to be applied in the model. For group selection, it is one of grLasso, grMCP, grSCAD or SGL while for bi-level selection it is gBridge or cMCP
#' @param dist This is the distribution for count model (Poisson or negative Binomial)
#' @param nlambda This is the number of lambda values. Default is 20 for SGL and 100 for other penalties.
#' @param lambda This is a user specified sequence of lambda values. By default it is left unspecified and the function automatically computes a grid of lambda values that ranges uniformly on the log scale over a relevant range of the values
#' @param lambda.min This is the smallest value of lambda, as a fraction of lambda.max. Default is .0001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param lambda.max This is the maximum value of lambda (needed for gBridge penalty only). Unlike the penalties in grpreg, it is not possible to solve for lambda.max directly with group bridge models. Thus, it must be specified by the user. If it is not specified, gBridge will attempt a guess at lambda.max which mihgt not be very accurate.
#' @param crit This is the selection criteria for the best model. It can be either "BIC" (default) or "AIC"
#' @param alpha This is the tuning parameter for the balance between the group penalty and the L2 penalty, as in grpreg. Default is 0.95 for SGL or 1 for other penalties
#' @param eps This is the convergence threshhold, as in grpreg.
#' @param gamma This is the tuning parameter of the group penalty (see details). Default is 0.8 for SGl, 0.5 for gBridge, 4 for SCAD and 3 for the rest
#' @param tau This is the tuning parameter for the group exponential lasso (GEL); defaults to 1/3
#' @param warn This is the logical argument which if specified as TRUE gives warning when the function fails to converge as in grpreg. DEFAULT is TRUE
#' @param delta This is required for the group bridge penalty only when it is not differentiable at zero, to bound it away from zero. There is typically no need to change this value
#' @param standardize	This is a logical flag for variable standardization prior to fitting the model
#' @param verbose This is the a	logical flag for whether or not step number will be output
#' @param step	This is the fitting parameter used for inital backtracking step size (between 0 and 1)
#' @param reset	This is the fitting parameter used for taking advantage of local strong convexity in nesterov momentum (number of iterations before momentum term is reset)
#' @return The function returns a list containing following outputs
#' @return coefficients A list containing the estimated coefficients for the count and the zero model
#' @return aic The AIC for the fitted model
#' @return bic The BIC for the fitted model
#' @return loglik The log likelihood for the fitted model
#' @export
#' @examples
#' ## Auto Insurance Claim Data
#' library(HDtweedie)
#' data("auto")
#' y<-auto$y
#' y<-round(y)
#' x<-auto$x
#' data<-cbind.data.frame(y,x)
#' group=c(rep(1,5),rep(2,7),rep(3,4),rep(4:14,each=3),15:21)
#' yvar<-names(data)[1]
#' xvars<-names(data)[-1]
#' zvars<-xvars
#'
#' ## ZIP regression
#' fit.poisson<-gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,samegrp.overlap=T,dist="poisson",penalty="gBridge")
#' fit.poisson$aic
#'
#' ## ZINB regression
#' fit.negbin<-gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,samegrp.overlap=T,dist="negbin",penalty="gBridge")
#' fit.negbin$aic



gooogle<-function(data,xvars,zvars,yvar,group=1:ncol(data),
                  samegrp.overlap=T,
                  penalty=c("grLasso", "grMCP", "grSCAD","gel", "cMCP",
                            "gBridge", "gLasso", "gMCP", "gBridge", "SGL"),
                  dist=c("poisson","negbin"),
                  nlambda=ifelse(penalty=="SGL",20,100), lambda,
                  lambda.min=ifelse(penalty=="SGL",0.1,ifelse(nrow(data[,unique(c(xvars,zvars))])>ncol(data[,unique(c(xvars,zvars))])),1e-4,.05), lambda.max,
                  crit="BIC", like="normal", alpha=ifelse(penalty=="SGL",0.95,1), eps=.001,
                  max.iter=1000, gmax=length(unique(group)),
                  gamma=ifelse(penalty=="SGL",0.8,ifelse(penalty=="gBridge",0.5,ifelse(penalty == "grSCAD", 4, 3))),tau = 1/3, warn=TRUE, delta=1e-7,thresh = 0.001,
                  standardize = T, verbose = FALSE,step = 1,
                  reset = 10, ...)
{
  ll.func<-function(beta.count,beta.zero,y,X,Z,dist)
  {
    y<-as.numeric(y[,1])

    if (dist=="negbin") theta<-a

    if (is.null(Z))
    {
      zgam<-rep(beta.zero,length(y))
    } else {
      zgam<-as.matrix(cbind(1,Z))%*%beta.zero
      zgam<-as.numeric(zgam[,1])
    }
    pzero<-exp(zgam)/(1+exp(zgam))

    xbet<-as.matrix(cbind(1,X))%*%beta.count
    xbet<-as.numeric(xbet[,1])
    mu<-exp(xbet)

    if (dist=="poisson"){
      ll<-try(sum(dZIP(y,mu=mu,sigma=pzero,log=T)),silent = T)
    } else {
      ll<-try(sum(dZINBI(y,mu=mu,sigma=1/theta,nu=pzero,log=T)),silent = T)
    }
    if (class(ll)=="try-error")
    {
      ll<-NA
    }
    return(ll)
  }

  y<-data.frame(data[,yvar])
  X<-data.frame(data[,xvars])
  pred.names<-union(xvars,zvars)
  names(X)<-paste(names(X),".count",sep="")
  names(y)<-yvar
  group.x<-group[which(pred.names %in% xvars)] #group of xvars
  n<-nrow(data)
  if (is.null(zvars)) #if there is no covariate in the zero model
  {
    Z<-NULL
    group.z<-NULL
    data<-cbind.data.frame(y,X)
    xvars<-names(X)
    fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|1"),sep=""))
  } else {
    Z<-data.frame(data[,zvars])
    names(Z)<-paste(names(Z),".zero",sep="")

    if (samegrp.overlap) # if X and Z assign same groups for shared covariates
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
  b2.mle<-c(fit.zero$coefficients$count[-1],fit.zero$coefficients$zero[-1])
  #coeff.intercept<-c(fit.zero$coefficients$count[1],fit.zero$coefficients$zero[1])
  vcov<-fit.zero$vcov
  p<-length(xvars)

  if (dist=="poisson")
  {
    sigma.11<-vcov[c(1,p+2),c(1,p+2)]
    sigma.12<-vcov[c(1,p+2),-c(1,p+2)]
    sigma.22<-vcov[-c(1,p+2),-c(1,p+2)]

    vcov.bar<-sigma.22-t(sigma.12)%*%ginv(sigma.11)%*%sigma.12

    e<-eigen(vcov.bar)
    if(det(vcov.bar)>0) # in case vcov is not pd add small values to the diagonal
    {
      cov.star=e$vectors%*%diag(1/sqrt(abs(e$values)))%*%t(e$vectors)
    } else {
      cov.star=e$vectors%*%diag(1/sqrt(PDFORCE(vcov.bar)))%*%t(e$vectors)
    }

    y.star<-cov.star%*%b2.mle # transformed y
    cov.star<-data.frame(cov.star) # scaled x matrix
    names(cov.star)<-c(xvars,zvars)

    ## reorder the group index and the covariates
    group<-c(group.x,group.z)
    u<-unique(group)
    cov.star.reordered<-rep(0,dim(cov.star)[1])
    for(i in 1:length(u))
    {
      index<-which(group==u[i])
      cov.star.reordered<-cbind(cov.star.reordered,cov.star[,index])
    }
    cov.star.reordered<-cov.star.reordered[,-1]
    group<-sort(group)

    if (penalty=="gBridge")
    {
      fit<-gBridge(X=cov.star.reordered, y=y.star, group=group, gamma=gamma)

    } else if (penalty=="SGL") {
      data.sgl<-list(x=cov.star.reordered,y=y.star)
      if(missing(lambda)) lambdas<-NULL else lambdas<-lambda

      #       browser()

      #       fit.sgl<-SGL(data.sgl, index=group, type = "linear", maxit = max.iter, thresh = thresh, min.frac = min.frac, nlam = nlambda, gamma = gamma, standardize = standardize, verbose = verbose, step = step, reset = reset, alpha = alpha, lambdas)
      fit.sgl<-SGL(data.sgl, index=group, type = "linear", maxit = max.iter, thresh = thresh, min.frac = lambda.min, nlam = nlambda, gamma = gamma, standardize = standardize, verbose = verbose, step = step, reset = reset, alpha = alpha, lambdas)
      fit<-list(fit.sgl=fit.sgl,y=y.star,x=cov.star.reordered)
    } else {
      fit=grpreg(cov.star.reordered,y.star,group=group,penalty=penalty,family="gaussian", nlambda=nlambda, lambda, lambda.min=lambda.min,  alpha=alpha, eps=eps, max.iter=max.iter,gmax=gmax, gamma=gamma, tau = tau,  warn=warn)
    }
  } else if(dist=="negbin") {
    # browser()
    sigma.11<-vcov[c(1,p+2),c(1,p+2)]
    sigma.12<-vcov[c(1,p+2),-c(1,p+2)]
    sigma.22<-vcov[-c(1,p+2),-c(1,p+2)]

    vcov.bar<-try(sigma.22-t(sigma.12)%*%ginv(sigma.11)%*%sigma.12,silent=T)

    a<-1/fit.zero$theta
    var.a<-(fit.zero$SE.logtheta)^2*exp(2*log(a))

    e<-eigen(vcov.bar)

    if(round(det(vcov.bar),4)>0) # transform vcov to pd if it's not
    {
      cov.star=e$vectors%*%diag(1/sqrt(abs(e$values)))%*%t(e$vectors)
    } else {
      cov.star=e$vectors%*%diag(1/sqrt(PDFORCE(vcov.bar)))%*%t(e$vectors)
    }
    y.star<-cov.star%*%b2.mle
    cov.star<-data.frame(cov.star) # scaled x matrix
    names(cov.star)<-c(xvars,zvars)

    group<-c(group.x,group.z) # 0 indicator correspond to disperison and intercept terms
    u<-unique(group)
    ## reorder the group index and the covariates
    cov.star.reordered<-rep(0,dim(cov.star)[1])
    for(i in 1:length(u))
    {
      index<-which(group==u[i])
      cov.star.reordered<-cbind(cov.star.reordered,cov.star[,index])
    }
    cov.star.reordered<-cov.star.reordered[,-1]
    group<-sort(group)

    #     browser()
    if (penalty=="gBridge")
    {
      fit<-gBridge(X=cov.star.reordered, y=y.star, group=group, gamma=gamma)

    } else if (penalty=="SGL") {
      data.sgl<-list(x=cov.star.reordered,y=y.star)
      if(missing(lambda)) lambdas<-NULL else lambdas<-lambda
      fit.sgl<-SGL(data.sgl, index=group, type = "linear", maxit = max.iter, thresh = thresh, min.frac = lambda.min, nlam = nlambda, gamma = gamma, standardize = standardize, verbose = verbose, step = step, reset = reset, alpha = alpha, lambdas)
      fit<-list(fit.sgl=fit.sgl,y=y.star,x=cov.star.reordered)
    } else {
      fit=grpreg(cov.star.reordered,y.star,group=group,penalty=penalty,family="gaussian", nlambda=nlambda, lambda, lambda.min=lambda.min,  alpha=alpha, eps=eps, max.iter=max.iter, gmax=gmax, gamma=gamma, tau = tau,  warn=warn)
    }
  } else {
    stop("Gooogle only works with Poisson or Negative Binomial distribution for count data")
  }

  if (like=="normal"){
    fit.final<-fit.final.func(fit=fit,penalty=penalty,dist=dist,crit=crit)
    b2.final<-c(fit.final$coefficients$count,fit.final$coefficients$zero)
    b1.final<-sigma.12%*%solve(sigma.22)%*%(b2.mle-b2.final)
    fit.final.coeff.count<-c(intercept=b1.final[1],fit.final$coefficients$count)
    fit.final.coeff.zero<-c(intercept=b1.final[2],fit.final$coefficients$zero)
    #   if (dist=="negbin")
    #   {
    #     fit.final$coefficients=c(fit.final$coefficients,dispersion=a)
    #   }

    dfc<-sum(fit.final.coeff.count!=0)
    dfz<-sum(fit.final.coeff.count!=0)
    ll<-ll.func(fit.final.coeff.count,fit.final.coeff.zero,y,X,Z,dist=dist)
    aic<--2*ll+2*(dfc+dfz)
    bic<--2*ll+log(n)*(dfc+dfz)


    if (dist=="poisson"){
      coeff.final=list(count=fit.final.coeff.count,zero=fit.final.coeff.zero)
    } else {
      coeff.final=list(count=fit.final.coeff.count,zero=fit.final.coeff.zero,dispersion=a)
    }

    return(list(coefficients=coeff.final,aic=aic, bic=bic, loglik=ll))
  } else if (like=="zip") {
#     browser()
    if (penalty=="SGL"){
      b2.final.mat<-fit.sgl$beta
      rownames(b2.final.mat)<-names(cov.star.reordered)
      b2.final.mat<-b2.final.mat[names(cov.star),]
    } else {
      b2.final.mat<-fit$beta[names(cov.star),]
    }
    b1.final.mat<-apply(b2.final.mat,2,function(x) return(sigma.12%*%solve(sigma.22)%*%(b2.mle-x)))
    count.mat<-rbind(b1.final.mat[1,],b2.final.mat[1:p,])
    rownames(count.mat)[1]<-"Intercept.count"
    rownames(count.mat)<-substr(rownames(count.mat),1,nchar(rownames(count.mat))-6)

    zero.mat<-rbind(b1.final.mat[2,],b2.final.mat[-(1:p),])
    rownames(zero.mat)[1]<-"Intercept.zero"
    rownames(zero.mat)<-substr(rownames(zero.mat),1,nchar(rownames(zero.mat))-5)


    crit.mat<-NULL

    for (i in 1:nlambda)
    {

      beta.count<-count.mat[,i]
      beta.zero<-zero.mat[,i]
      dfc<-sum(beta.count!=0)
      dfz<-sum(beta.zero!=0)
      ll<-ll.func(beta.count,beta.zero,y,X,Z,dist=dist)
      aic<--2*ll+2*(dfc+dfz)
      bic<--2*ll+log(n)*(dfc+dfz)
      crit.mat<-rbind.data.frame(crit.mat,c(ll,aic,bic))
    }
    names(crit.mat)<-c("LogLik","AIC","BIC")
    fit.final.idx<-which(crit.mat[,crit]==min(crit.mat[,crit],na.rm = T))

    fit.final.aic<-crit.mat[fit.final.idx,2]
    fit.final.bic<-crit.mat[fit.final.idx,3]
    fit.final.loglik<-crit.mat[fit.final.idx,1]

    fit.final.coeff.count<-count.mat[,fit.final.idx]
    fit.final.coeff.zero<-zero.mat[,fit.final.idx]

    if (dist=="poisson"){
      coeff.final=list(count=fit.final.coeff.count,zero=fit.final.coeff.zero)
    } else {
      coeff.final=list(count=fit.final.coeff.count,zero=fit.final.coeff.zero,dispersion=a)
    }
    return(list(coefficients=coeff.final,aic=fit.final.aic, bic=fit.final.bic, loglik=fit.final.loglik))
  }
}


