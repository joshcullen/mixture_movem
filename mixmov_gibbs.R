mixture_movement=function(dat,alpha,ngibbs,nmaxclust,nburn){
  nobs=nrow(dat)
  ndata.types=ncol(dat)
  
  #initial values
  z=sample(1:nmaxclust,size=nobs,replace=T)
  ncat.dat=apply(dat,2,max,na.rm=T)
  phi=list()
  for (i in 1:ndata.types){
    phi[[i]]=matrix(1/ncat.dat[i],nmaxclust,ncat.dat[i])
  }
  theta=rep(1/nmaxclust,nmaxclust)
  gamma1=0.1
  
  #get nmat
  nmat=list()
  for (i in 1:ndata.types){
    nmat[[i]]=SummarizeDat(z=z-1, dat=dat[,i]-1, ncateg=ncat.dat[i],nbehav=nmaxclust, nobs=nobs)
  }
  
  #prepare for gibbs
  store.phi=list()
  for (i in 1:ndata.types){
    store.phi[[i]]=matrix(NA,ngibbs,nmaxclust*ncat.dat[i])
  }
  store.theta=matrix(NA,ngibbs,nmaxclust)
  store.loglikel=rep(NA,ngibbs)
  store.gamma1=rep(NA,ngibbs)
  store.z<- list()
  
  #run gibbs sampler
  max.llk=-Inf
  gamma.possib=seq(from=0.1,to=1,by=0.05) #possible values for gamma
  
  #progress bar
  pb <- progress::progress_bar$new(
    format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
    total = ngibbs, clear = FALSE, width = 100)
  
  for (i in 1:ngibbs){
    pb$tick()  #create progress bar
    
    #sample from FCD's 
    lphi=list()
    for (j in 1:ndata.types) lphi[[j]]=log(phi[[j]])
    ltheta=log(theta)
    z=sample.z(nobs=nobs,nmaxclust=nmaxclust,dat=dat,ltheta=ltheta,lphi=lphi,ndata.types=ndata.types)
    for (j in 1:ndata.types){
      nmat[[j]]=SummarizeDat(z=z-1, dat=dat[,j]-1, ncateg=ncat.dat[j],nbehav=nmaxclust, nobs=nobs)
    }

    tmp=sample.v(z=z,gamma1=gamma1,nmaxclust=nmaxclust)
    theta=tmp$theta
    v=tmp$v
    
    phi=sample.phi(alpha=alpha,nmaxclust=nmaxclust,
                   ncat.dat=ncat.dat,ndata.types=ndata.types,nmat=nmat)
    
    gamma1=sample.gamma(v=v,ngroup=nmaxclust,gamma.possib=gamma.possib)    
    
    #calculate log-likelihood
    llk=get.llk(phi=phi,theta=theta,ndata.types=ndata.types,dat=dat,
                nobs=nobs,nmaxclust=nmaxclust)    
    
    #store results
    for (j in 1:ndata.types){
      store.phi[[j]][i,]=phi[[j]]
    }
    store.theta[i,]=theta
    store.loglikel[i]=llk
    store.gamma1[i]=gamma1
    store.z[[i]]=z
    
    #re-order clusters
    if (i < nburn & i%%50==0){
      ordem=order(theta,decreasing=T)
      theta=theta[ordem]
      
      for (j in 1:ndata.types){
        phi[[j]]=phi[[j]][ordem,]
        nmat[[j]]=nmat[[j]][ordem,]
      }
      
      znew=z
      for (j in 1:nmaxclust){
        cond=z==ordem[j]
        znew[cond]=j
      }
      z=znew
    }
  }

  # if (i > nburn & llk>max.llk){
  #   max.llk=llk
  #   
  #   theta.max.llk=theta
  #   names(theta.max.llk)<- 1:length(theta.max.llk)
  #   theta.max.llk<- theta.max.llk %>% sort(decreasing = T)
  #   ord<- names(theta.max.llk)
  #   
  #   phi.max.llk=phi
  #   phi.max.llk<- lapply(phi.max.llk, t)
  #   phi.max.llk<- lapply(phi.max.llk, function(x) x[,as.numeric(ord)])
  #   
  #   z.max.llk=z
  #   z.max.llk<- factor(z.max.llk)
  #   levels(z.max.llk)<- ord
  #   z.max.llk<- as.numeric(as.character(z.max.llk))
  # }
  
  list(phi=store.phi,theta=store.theta,
       loglikel=store.loglikel,z=store.z,
       gamma1=store.gamma1)
  # list(phi=phi.max.llk,theta=theta.max.llk,
  #      loglikel=store.loglikel,z=z.max.llk,
  #      gamma1=store.gamma1)  
}
