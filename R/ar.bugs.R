ar.bugs <-
function(y, ar.order=1, h=NULL, sim=FALSE, 
                  mean.centre=FALSE, st=ar.order+1,
                  mean.prior=ar.prior, ar.prior="dnorm(0,1)", tol.prior="dgamma(0.000001,0.000001)", var.prior=NULL, sd.prior=NULL){
  y<-c(y)
  n<-length(y)
  if(!is.null(h)){
    y<-c(y,rep(NA,h))
  }
  h<-length(y)-max(which(!is.na(y)))
  if(st<ar.order)
    stop("The value of st must be at least 1 greater than the number of lags")
  if(!is.null(var.prior) | !is.null(sd.prior)){
    tol.prior<-NULL
  }
  if(length(c(tol.prior,var.prior,sd.prior))>1)
    stop("Only one of tol.prior, var.prior or sd.prior should be given. Set others to null")
  
  bug<-c("model{","")
  #likelihood
  lik<-c("#likelihood",
         paste0("for(t in ",st,":",n+h,"){"),
         "\ty[t] ~ dnorm(y.mean[t], isigma2)",
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
  bug<-c(bug, ymean, "")
  
  #prior
  if(!is.null(tol.prior)){
    ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
                 paste0("isigma2 ~ ",tol.prior),
                 "sigma <- pow(isigma2,-0.5)")
  }
  if(!is.null(var.prior)){
    ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
                 paste0("sigma2 ~ ",var.prior),
                 "isigma2 <- pow(sigma2,-1)")
  }
  if(!is.null(sd.prior)){
    ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
                 paste0("sigma ~ ",sd.prior),
                 "isigma2 <- pow(sigma,-2)")
  }
  ar.priors<-ar.priors[-1]
  if(mean.centre==T)  ar.priors<-c(paste0("phi0 ~ ",mean.prior),ar.priors)
  bug<-c(bug,"#priors",ar.priors,"")
    
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
  
  #simulation
  ysim<-NULL
  if(sim==TRUE){
    ysim<-c("#simulations",
            "isigma2.c <- cut(isigma2)",
            paste("for(t in ",st,":",n,"){",sep=""),
            "\ty.mean.c[t] <- cut(y.mean[t])",
            "\ty.sim[t] ~ dnorm(y.mean.c[t],isigma2.c)",
            "}",
            "")
    bug<-c(bug,ysim)
  }
  bug<-c(bug,"}","")
  #print.tsbugs(list(bug=bug))
  
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
            info=list(n=n,h=h,nh=n+h,
                      variance="CV",
                      likelihood=p1:(p2-1),
                      priors=p2:(p3-1),
                      forecasts=NULL,
                      simulations=NULL))
  if(p3!=p4)  bug$info$forecasts<-p3:(p4-1)
  if(p4!=p5)  bug$info$simulations<-p4:(p5-1)
  class(bug)<-"tsbugs"
  return(bug)
}
