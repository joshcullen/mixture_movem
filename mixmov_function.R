sample.z=function(ncat.dat,dat, nmaxclust,
                  lphi,ltheta,ndata.types,nobs,z,alpha,phi){
  
  lprob=matrix(ltheta,nobs,nmaxclust,byrow=T)
  for (i in 1:nmaxclust){
    for (j in 1:ndata.types){
      lprob[,i]=lprob[,i]+lphi[[j]][i,dat[,j]]
    }
  }
  
  #get dirichlet densities
  dirichlet.den=matrix(NA,ndata.types,nmaxclust)
  for (j in 1:ndata.types){
    dirichlet.den[j,]=log(ddirichlet(phi[[j]],rep(alpha,ncat.dat[j])))
  }
  dirichlet.den1=colSums(dirichlet.den)
  
  p1=sum(-log(ncat.dat))
  for (i in 1:nobs){
    maxz=max(z)
    if (maxz==nmaxclust) lprob1=lprob[i,]
    if (maxz<nmaxclust){
      lprob1=lprob[i,1:maxz]+dirichlet.den1[maxz+1] #for existing groups
      tmp=p1+ltheta[maxz+1] #for new group
      lprob1=c(lprob1,tmp)
    }
    
    max1=max(lprob1)
    lprob1=lprob1-max1
    tmp=exp(lprob1)
    prob=tmp/sum(tmp)
    
    tmp=rmultinom(1, size=1,prob=prob) 
    z[i]=which(tmp==1)
  }
  z
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
