#include <Rcpp.h>
using namespace Rcpp;

// A function to compute the sum of weighted mean (1st column)
// and the sum of weights (2nd colunm)
// for different values of tau2

// [[Rcpp::export]]
NumericMatrix compute_rl_tau2(NumericVector x1,NumericVector x2, 
                              NumericVector x4){
  
  // x1 is the effect size g
  // x2 is the sampling variance vi
  // x4 is tau2
  int i; // the number of values of tau2
  int ii; // the length of g.sort and vi.sort
  int icum;
  NumericMatrix res(x4.length(),4);
  for (i = 0; i < x4.length(); i++) {
    double SumW = 0;
    double SumWY = 0;
    double CumsumW = 0;
    double CumsumWY = 0;
    for (ii = 0; ii < x2.length(); ii++) {
      double newW = x2[ii] + x4[i];
      SumW = SumW + 1/newW;
      SumWY = SumWY + x1[ii]/newW;
    }
    for (icum = 0; icum <= i; icum++) {
      CumsumW = CumsumW + 1/(x2[icum]+x4[i]);
      CumsumWY = CumsumWY +  x1[icum]/(x2[icum]+x4[i]);
    }
    res(i,1) = SumW;
    res(i,0) = SumWY;
    res(i,2) = CumsumWY;
    res(i,3) = CumsumW;
  }
  return res;
}