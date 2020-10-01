mixture_movement=function(dat,gamma1,alpha,ngibbs,nmaxclust,nburn){
  nobs=nrow(dat)
  ndata.types=ncol(dat)
  
  #initial values
  z=sample(1:nmaxclust,size=nobs,replace=T)
  ncat.dat=apply(dat,2,max)
  phi=list()
  for (i in 1:ndata.types){
    phi[[i]]=matrix(1/ncat.dat[i],nmaxclust,ncat.dat[i])
  }
  theta=rep(1/nmaxclust,nmaxclust)
  
  #prepare for gibbs
  store.phi=list()
  for (i in 1:ndata.types){
    store.phi[[i]]=matrix(NA,ngibbs,nmaxclust*ncat.dat[i])
  }
  store.theta=matrix(NA,ngibbs,nmaxclust)
  store.loglikel=rep(NA,1)
  
  #run gibbs sampler
  for (i in 1:ngibbs){
    print(i)
    print(table(z))
    
    #re-order clusters
    if (i < nburn & i%%50==0){
      ordem=order(theta,decreasing=T)
      theta=theta[ordem]
      
      for (j in 1:ndata.types){
        phi[[j]]=phi[[j]][ordem,]
      }
      
      znew=z
      for (j in 1:nmaxclust){
        cond=z==ordem[j]
        znew[cond]=j
      }
      z=znew
    }
    
    #sample from FCD's 
    lphi=list()
    for (j in 1:ndata.types) lphi[[j]]=log(phi[[j]])
    z=sample.z(ncat.dat=ncat.dat,dat=dat,nmaxclust=nmaxclust,
               lphi=lphi,ltheta=log(theta),ndata.types=ndata.types,nobs=nobs,
               z=z,alpha=alpha,phi=phi)

    theta=sample.v(z=z,gamma1=gamma1,nmaxclust=nmaxclust)
    # theta=theta.true

    phi=sample.phi(alpha=alpha,nmaxclust=nmaxclust,
                   ncat.dat=ncat.dat,ndata.types=ndata.types,dat=dat,z=z)
    
    #calculate log-likelihood
    llk=get.llk(phi=phi,theta=theta,ndata.types=ndata.types,dat=dat,
                nobs=nobs,nmaxclust=nmaxclust)    
    
    #store results
    for (j in 1:ndata.types){
      store.phi[[j]][i,]=phi[[j]]
    }
    store.theta[i,]=theta
    store.loglikel[i]=llk
  }

  list(phi=store.phi,theta=store.theta,
       loglikel=store.loglikel,z=z)  
}
