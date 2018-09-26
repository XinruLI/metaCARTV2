#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


struct node {
  NumericVector rows;
  node *left;
  node *right;
} Nnode, *pnode;

// [[Rcpp::export]]
NumericVector LeftChildNode(NumericVector x, double s) {
  node leftson, rightson;
  leftson.rows = x[x < s];
  rightson.rows = x [x >= s];
  node root;
  root.rows = x;
  *root.left = leftson;
  *root.right = rightson;
  return root.rows;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
LeftChildNode(1:10, 4)
*/
