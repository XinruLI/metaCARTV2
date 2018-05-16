#include <Rcpp.h>
using namespace Rcpp;

// This is a function to compute re Q for different values of tau2.
// The first column of the results is the sum of weighted mean
// The seconde column of the results is the sum of weights


// [[Rcpp::export]]
NumericVector compute_re_Q_(NumericVector x1,NumericVector x2, 
                               NumericVector x3,NumericVector x4,
                               NumericVector xuni,
                               NumericVector x5,NumericVector x6){
  // x1 is the effect size g in the unsplit leaves
  // x2 is the sampling variance vi in the unsplit leaves
  // x3 is the labels of nodes in the unsplit leaves
  // x4 is tau2
  // xuni is the unique labels in the unsplit leaves
  // x5 is the effect size g in the parent leaf
  // x6 is the sampling variance vi in the parent leaf
  
  int i; // the number of nodes
  int j; // the number of values of tau2
  NumericVector SWYoleaf;
  NumericVector SWoleaf;
  NumericVector SWY2byWoleaf;
  NumericVector SumWpleaf;
  NumericVector SumWYpleaf;
  NumericVector CumsumWYpleaf;
  NumericVector CumsumWpleaf;
  for (j = 0; j < x4.length(); j++) {
    double SumWYleft = 0;
    double SumWleft = 0;
    double SumWY2byW = 0;
    // in the other leaves
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
    SWYoleaf.push_back(SumWYleft);
    SWoleaf.push_back(SumWleft);
    SWY2byWoleaf.push_back(SumWY2byW);
    // in the parent leaf
    double SumW = 0;
    double SumWY = 0;
    double CumsumW = 0;
    double CumsumWY = 0;
    int k;
    int kcum;
    for (k = 0; k < x6.length(); k++) {
      double newW = x6[k] + x4[j];
      SumW = SumW + 1/newW;
      SumWY = SumWY + x5[k]/newW;
    }
    for (kcum = 0; kcum <= j; kcum++) {
      CumsumW = CumsumW + 1/(x6[kcum]+x4[j]);
      CumsumWY = CumsumWY +  x5[kcum]/(x6[kcum]+x4[j]);
    }
    SumWpleaf.push_back(SumW);
    SumWYpleaf.push_back(SumWY);
    CumsumWYpleaf.push_back(CumsumWY);
    CumsumWpleaf.push_back(CumsumW);
  }
  NumericVector res;
  for (j = 0; j < x4.length(); j++) {
    double tempSSWY = SumWYpleaf[j] + SWYoleaf[j];
    double tempSSW = SumWpleaf[j] + SWoleaf[j];
    res.push_back(
      pow(SumWYpleaf[j]-CumsumWYpleaf[j],2)/(SumWpleaf[j]-CumsumWpleaf[j])+SWY2byWoleaf[j]+pow(CumsumWYpleaf[j],2)/CumsumWpleaf[j]-pow(tempSSWY,2)/tempSSW);
    }
  
  return res;
}

