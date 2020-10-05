sample.z=function(nobs,z,nmat,dat,alpha,ncat.dat,ltheta){

  #this is the same for all data types because it  is the number of observations assigned to each behavior
  ntot=rowSums(nmat[[1]]) 
  randu=runif(nobs)
  
  res.cpp=sampleZ(nmat1=nmat[[1]], nmat2=nmat[[2]], z=z-1, dat=data.matrix(dat)-1,
              ntot=ntot, ltheta=ltheta, randu=randu, NcatDat=ncat.dat,
              nobs=nobs, nmaxclust=nmaxclust, alpha=alpha)
  nmat=list()
  nmat[[1]]=res.cpp$nmat1
  nmat[[2]]=res.cpp$nmat2

  # teste=SummarizeDat(z=res.cpp$z, dat=dat[,1]-1, ncateg=ncat.dat[1],nbehav=nmaxclust, nobs=nobs)
  # unique(teste-res.cpp$nmat1)
  
  list(nmat=nmat,z=res.cpp$z+1)
}
#-----------------------------------
sample.v=function(nmat,gamma1,nmaxclust){
  tmp1=rowSums(nmat[[1]])
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
      prob[,i]=prob[,i]*phi[[j]][i,dat[,j]]
    }
  }
  sum(log(rowSums(prob)))
}
