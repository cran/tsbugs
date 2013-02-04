rv.bugs <-
function(y, ar.order=0, h=NULL, sim=FALSE, 
                  mean.centre=FALSE, st=ar.order+1,
                  mean.prior=ar.prior, ar.prior="dnorm(0,1)",
                  rv.tol0.prior="dgamma(0.000001,0.000001)", rv.eps.prior="dbeta(1, 100)", rv.ilambda2.prior="dgamma(0.01,0.01)"){
  y<-c(y)
  n<-length(y)
  if(!is.null(h)){
    y<-c(y,rep(NA,h))
  }
  h<-length(y)-max(which(!is.na(y)))
  if(st<ar.order)
    stop("The value of st must be at least 1 greater than the number of lags")
  
  bug<-c("model{","")
  #likelihood
  lik<-c("#likelihood",
         paste0("for(t in ",st,":",n+h,"){"),
         "\ty[t] ~ dnorm(y.mean[t], isigma2[t])",
         "\tisigma2[t] <- exp(-2*lsig[t])",
         "\th[t] <- 2*lsig[t]",
         "}")
  bug<-c(bug, lik)
  
  #ymean
  ymean<-c("#mean",
           paste0("for(t in ",st,":",n+h,"){"),
           y.mean<-c("\ty.mean[t] <- 0",
                     "}")
  )
  if(ar.order==0 & mean.centre==T)  ymean[3]<-"\ty.mean[t] <- phi0"
  if(ar.order!=0 & mean.centre==F)  ymean[3]<-paste0("\ty.mean[t] <- ",paste0("phi",1:ar.order,"*y[t-",1:ar.order,"]",collapse=" + "))
  if(ar.order!=0 & mean.centre==T)  ymean[3]<-paste0("\ty.mean[t] <- phi0 + ",paste0("phi",1:ar.order,"*(y[t-",1:ar.order,"]-phi0)",collapse=" + "))
  bug<-c(bug, ymean)
  
  #rv
  rv<-c("#variance",
        paste0("lsig[",st,"] <- -0.5*log(isig02)"),
        paste0("for(t in ",st+1,":",n+h,"){"),
        "\tlsig[t] <- lsig[t-1]+(delta[t]*beta[t])",
        "\tdelta[t] ~ dbern(epsilon)",
        "\tbeta[t] ~ dnorm(0,ilambda2)",
        "}",
        "")
  bug<-c(bug, rv)
  
  #priors
  ar.priors<-paste0("phi",0:ar.order," ~ ",ar.prior)
  ar.priors<-ar.priors[-1]
  if(mean.centre==T)  ar.priors<-c(paste0("phi0 ~ ",mean.prior),ar.priors)  
  rv.priors<-c(paste0("isig02 ~ ",rv.tol0.prior),
               "sig0 <- pow(isig02,-0.5)",
               paste0("epsilon ~ ",rv.eps.prior),
               paste0("ilambda2 ~ ",rv.ilambda2.prior),
               "lambda <- pow(ilambda2,-0.5)")
  bug<-c(bug,"#priors",ar.priors,rv.priors,"")
  
  #forecast
  forc<-NULL
  if(h!=0){
    forc<-c("#forecasts",
            paste("for(t in ",n+1,":",n+h,"){",sep=""),
            "\ty.new[t] <- y[t]",
            "}",
            "")
    bug<-c(bug,forc)
  }
  
  #simulations
  if(sim==TRUE){
    ysim<-c("#simulations",
            paste("for(t in ",st,":",n,"){",sep=""),
            "\ty.mean.c[t] <- cut(y.mean[t])",
            "\tisigma2.c[t] <- cut(isigma2[t])",
            "\ty.sim[t] ~ dnorm(y.mean.c[t],isigma2.c[t])",
            "}",
            "")
    bug<-c(bug,ysim)
  }
  bug<-c(bug,"}","")
  
  p1<-grep("likelihood",bug)
  p2<-grep("prior",bug)
  if(h!=0 & sim==TRUE){
    p3<-grep("forecast",bug); p4<-grep("simulation",bug)
  }
  if(h!=0 & sim==FALSE){
    p3<-grep("forecast",bug); p4<-length(bug)
  }
  if(h==0 & sim==TRUE){
    p3<-grep("simulation",bug); p4<-p3
  } 
  if(h==0 & sim==FALSE){
    p3<-length(bug); p4<-p3
  } 
  p5<-length(bug)
  
  bug<-list(bug=bug,
            data=list(y=y),
            info=list(n=n,h=h,nh=n+h,st=st,
                      args=mget(names(formals()),sys.frame(sys.nframe()))[-1],
                      variance="RV",
                      likelihood=p1:(p2-1),
                      priors=p2:(p3-1),
                      forecasts=NULL,
                      simulations=NULL))
  if(p3!=p4)  bug$info$forecasts<-p3:(p4-1)
  if(p4!=p5)  bug$info$simulations<-p4:(p5-1)
  class(bug)<-"tsbugs"
  return(bug)
  class(bug)<-"tsbugs"
  return(bug)
}
