sample.z=function(nobs,nmaxclust, dat, ltheta, lphi, ndata.types){

  lprob=matrix(ltheta,nobs,nmaxclust,byrow=T)
  for (i in 1:nmaxclust){
    for (j in 1:ndata.types){
      #account for NA
      tmp=lphi[[j]][i,dat[,j]]
      cond=is.na(dat[,j])
      tmp[cond]=0
      
      #finish calculation
      lprob[,i]=lprob[,i]+tmp
    }
  }
  max1=apply(lprob,1,max)
  lprob=lprob-max1
  tmp=exp(lprob)
  prob=tmp/rowSums(tmp)
  
  z=rmultinom1(prob=prob, randu=runif(nobs)) 
  z+1
}
#-----------------------------------
sample.v=function(z,gamma1,nmaxclust){
  tmp=table(z)
  tmp1=rep(0,nmaxclust)
  tmp1[as.numeric(names(tmp))]=tmp
  theta=v=rep(NA,nmaxclust)
  aux=1
  for (i in 1:(nmaxclust-1)){
    seq1=(i+1):nmaxclust
    v[i]=rbeta(1,tmp1[i]+1,sum(tmp1[seq1])+gamma1)
    theta[i]=v[i]*aux
    aux=aux*(1-v[i])
  }
  v[nmaxclust]=1
  theta[nmaxclust]=aux
  theta
}
#-----------------------------------
sample.phi=function(alpha,nmaxclust,ncat.dat,ndata.types,nmat){
  phi=list()
  for (j in 1:ndata.types){
    tmp2=matrix(NA,nmaxclust,ncat.dat[j])
    for (i in 1:nmaxclust){
      tmp2[i,]=rdirichlet(1,nmat[[j]][i,]+alpha)
    }
    phi[[j]]=tmp2
  }
  phi
}
#------------------------------------
get.llk=function(phi,theta,ndata.types,dat,nobs,nmaxclust){
  prob=matrix(theta,nobs,nmaxclust,byrow=T)
  for (i in 1:nmaxclust){
    for (j in 1:ndata.types){
      #account for NA      
      tmp=phi[[j]][i,dat[,j]]
      cond=is.na(dat[,j])
      tmp[cond]=1
      
      #finish calculation
      prob[,i]=prob[,i]*tmp
    }
  }
  sum(log(rowSums(prob)))
}
