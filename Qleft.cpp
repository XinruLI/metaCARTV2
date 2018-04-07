#include <Rcpp.h>
using namespace Rcpp;

// This is a function to compute the within-nodes Q

// [[Rcpp::export]]
NumericVector compute_Q_left(NumericVector x1,NumericVector x2, NumericVector x3, NumericVector xuni){
  int i;
  NumericVector res;
  for (i = 0; i < xuni.length(); i++) {
    double tempWY2 = 0;
    double tempWY = 0;
    double tempW = 0;
    int ii;
    for (ii = 0; ii < x3.length(); ii++) {
      if (x3[ii] == xuni[i]) {
        tempWY2 = tempWY2 + pow(x1[ii], 2)/x2[ii];
        tempWY = tempWY + x1[ii]/x2[ii];
        tempW = tempW + 1/x2[ii];
      }
    }
    res.push_back(tempWY2 - pow(tempWY,2)/tempW);
  }
  return res;
}






