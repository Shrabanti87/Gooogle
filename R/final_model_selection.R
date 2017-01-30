########## selecting lambda and the final model selection ####################

fit.final.func<-function(fit,penalty,dist,crit)
{
  #   browser()
  if (penalty!="SGL")
  {
    fit.final.idx<-which.min(get(crit)(fit))
    fit.final.aic<-AIC(fit)[fit.final.idx]
    fit.final.bic<-BIC(fit)[fit.final.idx]
    fit.final.loglik<-logLik(fit)[fit.final.idx]
    fit.final.coeff<-fit$beta[,fit.final.idx]
  } else {
    fit.sgl<-fit$fit.sgl
    y<-fit$y
    x<-fit$x

    nlambda<-length(fit.sgl$lambdas)
    n<-nrow(fit.sgl$beta)
    H<-as.matrix(fit.sgl$beta)

    # get min BIC
    crit.mat<-rep(NULL,3)
    for (i in 1:nlambda){
      beta<-c(fit.sgl$intercept,H[,i])
      df<-length(beta[beta!=0])
      mu<-as.matrix(cbind(1,x))%*%beta
      mu<-as.numeric(mu[,1])
      loglik<-sum(dnorm(y,mean=mu,log=T))
      aic<--2*loglik + 2*df
      bic= -2*loglik + log(n)*df;
      crit.mat<-rbind.data.frame(crit.mat,c(bic,aic,loglik));
    }
    names(crit.mat)<-c("BIC","AIC","LogLik")
    fit.final.idx<-which.min(crit.mat[,crit])

    fit.final.aic<-crit.mat[fit.final.idx,2]
    fit.final.bic<-crit.mat[fit.final.idx,1]
    fit.final.loglik<-crit.mat[fit.final.idx,3]

    rownames(H)<-names(fit$x)
    fit.final.coeff<-H[,fit.final.idx]
  }

  fit.final.coeff.count<-fit.final.coeff[str_sub(names(fit.final.coeff),start = -5)=="count"]
  names(fit.final.coeff.count)<-substr(names(fit.final.coeff.count),1,nchar(names(fit.final.coeff.count))-6)
  fit.final.coeff.zero<-fit.final.coeff[str_sub(names(fit.final.coeff),start = -4)=="zero"]
  names(fit.final.coeff.zero)<-substr(names(fit.final.coeff.zero),1,nchar(names(fit.final.coeff.zero))-5)

  return(list(coefficients=list(count=fit.final.coeff.count,zero=fit.final.coeff.zero),aic=fit.final.aic, bic=fit.final.bic, loglik=fit.final.loglik))
}
