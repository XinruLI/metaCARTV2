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

// [[Rcpp::export]]
int test_rc(int n){
  /* This is called the base condition, it is
   * very important to specify the base condition
   * in recursion, otherwise your program will throw
   * stack overflow error.
   */
  if (n <= 1)
    return 1;
  else 
    return n*test_rc(n-1);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
test_rc(4)
*/
