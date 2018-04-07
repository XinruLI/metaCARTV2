#include <Rcpp.h>
using namespace Rcpp;
// A function to compute C for each node

// [[Rcpp::export]]
NumericVector compute_C_left(NumericVector x2, NumericVector x3, NumericVector xuni){
  int i;
  NumericVector res;
  // tapply the group means
  for (i = 0; i < xuni.length(); i++) {
    double tempW2 = 0;
    double tempW = 0;
    int ii;
    for (ii = 0; ii < x3.length(); ii++) {
      if (x3[ii] == xuni[i]) {
        tempW2 = tempW2 + pow(x2[ii], -2);
        tempW = tempW + 1/x2[ii];
      }
    }
    res.push_back(tempW2/tempW);
  }
  return res;
}
