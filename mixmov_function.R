sample.z=function(nobs,z,nmat,dat,alpha,ncat.dat,ltheta){

  #this is the same for all data types because it  is the number of observations assigned to each behavior
  ntot=rowSums(nmat[[1]]) 
  randu=runif(nobs)
  
  res.cpp=sampleZ(nmat1=nmat[[1]], nmat2=nmat[[2]], z=z-1, dat=data.matrix(dat)-1,
              ntot=ntot, ltheta=ltheta, randu=randu, NcatDat=ncat.dat,
              nobs=nobs, nmaxclust=nmaxclust, alpha=alpha)
  # 
  # for (i in 1:nobs){
  #   #subtract that individual
  #   nmat[[1]][z[i],dat[i,1]]=nmat[[1]][z[i],dat[i,1]]-1
  #   nmat[[2]][z[i],dat[i,2]]=nmat[[2]][z[i],dat[i,2]]-1
  #   ntot[z[i]]=ntot[z[i]]-1 
  #   
  #   #calculate lprob
  #   p1=log(nmat[[1]][,dat[i,1]]+alpha)
  #   p2=log(nmat[[2]][,dat[i,2]]+alpha)
  #   p3=log(ntot+ncat.dat[1]*alpha)
  #   p4=log(ntot+ncat.dat[2]*alpha)
  #   lprob=p1+p2-p3-p4+ltheta
  #   
  #   #sample from multinomial
  #   max1=max(lprob)
  #   lprob1=lprob-max1
  #   tmp=exp(lprob1)
  #   prob=tmp/sum(tmp)
  #   z[i]=cat1(randu[i],prob)+1
  #   
  #   #add that individual
  #   nmat[[1]][z[i],dat[i,1]]=nmat[[1]][z[i],dat[i,1]]+1
  #   nmat[[2]][z[i],dat[i,2]]=nmat[[2]][z[i],dat[i,2]]+1
  #   ntot[z[i]]=ntot[z[i]]+1 
  # }
  
  # fim=data.frame(zcpp=res.cpp$z+1,zr=z)
  # table(fim)
  # which(fim$zcpp!=fim$zr)
  nmat=list()
  nmat[[1]]=res.cpp$nmat1
  nmat[[2]]=res.cpp$nmat2
  z=res.cpp$z+1
  list(nmat=nmat,z=z)
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
