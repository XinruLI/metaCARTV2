#include <Rcpp.h>
using namespace Rcpp;

// This is a function to compute sum of weighted means and sum of weights
// for different values of tau2.
// The first column of the results is the sum of weighted mean
// The seconde column of the results is the sum of weights


// [[Rcpp::export]]
NumericMatrix compute_swy_tau2(NumericVector x1,NumericVector x2, 
                               NumericVector x3,NumericVector x4,
                               NumericVector xuni){
  // x1 is the effect size g
  // x2 is the sampling variance vi
  // x3 is the labels of nodes
  // x4 is tau2
  int i; // the number of nodes
  int j; // the number of values of tau2
  NumericMatrix res(x4.length(),3);
  for (j = 0; j < x4.length(); j++) {
    double SumWYleft = 0;
    double SumWleft = 0;
    double SumWY2byW = 0;
    // Compute the sum of the RE within-subgroups Q
    for (i = 0; i < xuni.length(); i++) {
      double tempWY = 0;
      double tempW = 0;
      int ii;
      for (ii = 0; ii < x3.length(); ii++) {
        if (x3[ii] == xuni[i]) {
          double newW = x2[ii] + x4[j];
          tempWY = tempWY + x1[ii]/newW;
          tempW = tempW + 1/newW;
        }
      }
      SumWYleft = SumWYleft + tempWY;
      SumWleft = SumWleft + tempW;
      SumWY2byW = SumWY2byW + pow(tempWY,2)/tempW;
    } 
    res(j,0) = SumWYleft;
    res(j,1) = SumWleft;
    res(j,2) = SumWY2byW;
  }
  return res;
}

