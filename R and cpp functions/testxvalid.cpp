#include <Rcpp.h>
using namespace Rcpp;

// This is a function to compute the subgroup effect sizes
// for a tree at each depth
// [[Rcpp::export]]
DataFrame ComputeY(DataFrame x1, NumericVector y,
                NumericVector vi, NumericVector tau2) {
  // x1 is the node labels for each study
  // y is the effect size
  // vi is the sampling variance
  // tau2 is the tau2
  int nsplit;
  int i;
  int j;
  DataFrame res;
  for (nsplit = 0; nsplit < x1.ncol(); nsplit++) {
    IntegerVector Nodes = x1[nsplit];
    IntegerVector uniNodes = Rcpp::sort_unique(Nodes); 
    NumericVector swyNodes;
    for (i = 0; i < uniNodes.length(); i++) {
      // compute the weighted sum for each node
      double sumWY = 0;
      double sumW = 0;
      for (j = 0; j < y.length(); j++) {
        if (Nodes[j] == uniNodes[i]) {
          sumWY = sumWY + y[j]/(vi[j]+tau2[nsplit]);
          sumW = sumW + 1/(tau2[nsplit]+vi[j]);
        }
        
      }
     swyNodes.push_back(sumWY/sumW, std::to_string(uniNodes[i]));
      
    }
    
    res.push_back(swyNodes);
  }
  
  return res;

}

// A function to predict effect size for the test set
// [[Rcpp::export]]
NumericMatrix PredY(List x1, IntegerMatrix x2) {
  //x1 is the list of subgroup means
  //x2 is predicted subgroup membership for the test set
  NumericMatrix res(x2.nrow(), x2.ncol());
  int i;
  int j;
  for (i = 0; i < x2.ncol(); i++){
    NumericVector SubMeans = x1[i];
    LogicalVector tempB = Rcpp::is_na(x2(_,i));
    for (j = 0; j < x2.nrow(); j++) {
      if (tempB[j]) {
        res(j,i) = NA_REAL;
      } else {
        res(j,i) = SubMeans[std::to_string(x2(j,i))];
      }
     
      //res(j,i) = x2(j,i);
    }
  }
  return res;
}

// A function to replace missing values by the overall weighted mean
// [[Rcpp::export]]
NumericMatrix ReplaceNA(IntegerMatrix x1, NumericMatrix x2,
                        NumericVector y, NumericVector vi,
                        NumericVector tau2){
  // x1 is the two-column matrix of the indices of missing values
  // x2 is the matrix of predicted y with missing values
  int j;
  int ii;
  int jj;
  IntegerVector jUni = Rcpp::unique(x1(_,1));
  NumericVector OverallWM;
  for (jj = 0; jj < jUni.length(); jj++) { // compute the weighted means
    double temp1 = 0;
    double temp2 = 0;
    for (ii = 0; ii < y.length(); ii++) {
      temp1 = temp1 + y[ii]/(vi[ii] + tau2[jUni[jj]-1]);
      temp2 = temp2 + 1/(vi[ii] + tau2[jUni[jj]-1]);
    }
    OverallWM.push_back(temp1/temp2, std::to_string(jUni[jj]));
  }
  for (j = 0; j < x1.nrow(); j++) {
  x2(x1(j,0)-1, x1(j,1)-1)  = OverallWM[std::to_string(x1(j,1))];
  }
  return x2;
  
}


