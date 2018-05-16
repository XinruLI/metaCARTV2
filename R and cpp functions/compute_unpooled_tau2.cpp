#include <Rcpp.h>
using namespace Rcpp;

// A FUNCTION TO COMPUTE UNPOOLED TAU2 (SEPERATE ESTIMATES IN EACH NODE)
//

// [[Rcpp::export]]
NumericVector compute_unpooled_tau2(NumericVector y,
                                    NumericVector vi) {
  NumericVector wy;
  NumericVector wy2;
  NumericVector wts;
  NumericVector w2;
  int df = y.length() - 2;
  int j;
  for (j = 0; j < y.length(); j++) {
    wy.push_back(y[j]/vi[j]);
    wy2.push_back(pow(y[j],2)/vi[j]);
    wts.push_back(1/vi[j]);
    w2.push_back(pow(vi[j], -2));
  }
  NumericVector cwy = cumsum(wy);
  NumericVector cwy2 = cumsum(wy2);
  NumericVector cwts = cumsum(wts);
  NumericVector cw2 = cumsum(w2);
  NumericVector tau2;
  double swy = sum(wy);
  double swy2 = sum(wy2);
  double sw = sum(wts);
  double sw2 = sum(w2);
  for (j = 0; j < y.length() - 1; j++){
    double tempQ = swy2-pow(cwy[j],2)/cwts[j]-pow(swy-cwy[j],2)/(sw-cwts[j]);
    double tempC = sw-cw2[j]/cwts[j]-(sw2 - cw2[j])/(sw - cwts[j]);
    double temptau2 = (tempQ-df)/tempC;
    if (temptau2 < 0) {
      temptau2 = 0;
    }
    tau2.push_back(temptau2);
  }
  return tau2;
}


