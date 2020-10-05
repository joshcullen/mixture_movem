// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
// [[Rcpp::export]]
IntegerMatrix SummarizeDat(IntegerVector z, IntegerVector dat, int ncateg, 
                           int nbehav, int nobs){
  IntegerMatrix nmat(nbehav,ncateg);
  for (int i = 0; i < nobs; i++){
    nmat(z[i],dat[i])=nmat(z[i],dat[i])+1;
  }
  return nmat;
}

// This function helps with multinomial draws
// [[Rcpp::export]]
int cat1(double value, NumericVector prob) {
  int res=prob.length()-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function samples z's from a categorical distribution
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericMatrix prob, NumericVector randu) {
  
  IntegerVector z(prob.nrow());
  
  for(int i=0; i<prob.nrow();i++){
    z[i]=cat1(randu[i],prob(i,_));
  }
  return z;
}

//' This function samples z's from a categorical distribution
// [[Rcpp::export]]
List sampleZ(IntegerMatrix nmat1, IntegerMatrix nmat2,
                      IntegerVector z, IntegerMatrix dat,
                      IntegerVector ntot, NumericVector ltheta,
                      NumericVector randu, IntegerVector NcatDat,
                      int nobs, int nmaxclust, double alpha) {
  
  IntegerMatrix nmat1a=clone(nmat1);
  IntegerMatrix nmat2a=clone(nmat2);
  IntegerVector ntot1=clone(ntot);
  NumericVector p1(nmaxclust);
  NumericVector p2(nmaxclust);
  NumericVector p3(nmaxclust);
  NumericVector p4(nmaxclust);
  NumericVector lprob(nmaxclust);
  double max1;
  
  for(int i=0; i<nobs;i++){
    //substract that individual
    nmat1a(z[i],dat(i,0))=nmat1a(z[i],dat(i,0))-1;
    nmat2a(z[i],dat(i,1))=nmat2a(z[i],dat(i,1))-1;
    ntot1[z[i]]=ntot1[z[i]]-1;
    
    //calculate lprob
    p1=log(nmat1a(_,dat(i,0))+alpha);
    p2=log(nmat2a(_,dat(i,1))+alpha);
    p3=log(ntot1+NcatDat[0]*alpha);
    p4=log(ntot1+NcatDat[1]*alpha);
    lprob=p1+p2-p3-p4+ltheta;
    
    //sample from multinomial
    max1=max(lprob);
    lprob=lprob-max1;
    lprob=exp(lprob);
    lprob=lprob/sum(lprob);
    z[i]=cat1(randu[i],lprob);
    
    //add that individual
    nmat1a(z[i],dat(i,0))=nmat1a(z[i],dat(i,0))+1;
    nmat2a(z[i],dat(i,1))=nmat2a(z[i],dat(i,1))+1;
    ntot1[z[i]]=ntot1[z[i]]+1;
  }
  
  List L = List::create(Named("nmat1") = nmat1a,
                        Named("nmat2") = nmat2a,
                        Named("ntot") = ntot1,
                        Named("z") =  z);
  return L;
}
