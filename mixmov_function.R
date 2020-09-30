sample.z=function(ncat.dat,dat, nmaxclust,
                  lphi,ltheta,ndata.types,nobs){
  
  lprob=matrix(ltheta,nobs,nmaxclust,byrow=T)
  for (i in 1:nmaxclust){
    for (j in 1:ndata.types){
      lprob[,i]=lprob[,i]+lphi[[j]][i,dat[,j]]
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
sample.phi=function(alpha,nmaxclust,ncat.dat,ndata.types,dat,z){
  phi=list()
  for (j in 1:ndata.types){
    tab1=matrix(0,ncat.dat[j],nmaxclust)
    tmp=data.frame(y=dat[,j],z=z)
    tmp1=table(tmp)
    tab1[,as.numeric(colnames(tmp1))]=tmp1
    
    tmp2=matrix(NA,nmaxclust,ncat.dat[j])
    for (i in 1:nmaxclust){
      tmp2[i,]=rdirichlet(1,tab1[,i]+alpha)
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
      prob[,i]=prob[,i]*phi[[j]][i,dat[,j]]
    }
  }
  sum(log(rowSums(prob)))
}
