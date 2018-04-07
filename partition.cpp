#include <Rcpp.h>
using namespace Rcpp;


// A function to test if x1 contains the first element of x2 
bool contain(CharacterVector x1, CharacterVector x2){
  return std::find(x1.begin(), x1.end(), x2[0]) != x1.end();
}

// A function to partition the test set based on a trained tree
// [[Rcpp::export]]
IntegerMatrix partition(DataFrame x1, DataFrame x2, 
                    LogicalVector x3, IntegerVector x4,
                    List x5, DataFrame x6) {
  // x1 is the tree component of the REmrt object
  // x2 is the moderators in the test set
  // x3 indicates whether a moderator is numeric or not
  // x4 is the index vector of the spliting moderators
  // x5 is the list of split points
  // x6 is the moderators in the training set
  IntegerVector pleaf = x1["pleaf"];
  CharacterVector mod = x1["mod"];
  CharacterVector NewModNames = x2.names();
  IntegerVector pnode;
  IntegerMatrix res(x2.nrows(), x1.nrows());
  int i; 
  // int j;
  for (i = 0; i < x2.nrows(); i++) {// put all observations in the root node
    pnode.push_back(1);
    }
  res(_, 0) = pnode;
  int j;
  for (j = 1; j < x1.nrows(); j++) {
    if (x3[j] == true) {
      NumericVector sv = x2[x4[j]-1]; // spliting moderators in new data
      NumericVector tempSP = x5[j-1]; // split points
      for (i = 0; i < x2.nrows(); i++) {
        if (pnode[i] == pleaf[j]) {
          pnode[i] = 2*j;
          if (sv[i] > tempSP[0]) {
            pnode[i] = pnode[i]+1;
          }
        }

      }
    } else {
      CharacterVector sv = x2[x4[j]-1];
      CharacterVector tempSP = x5[j-1];
      CharacterVector tempModOld = x6[x4[j]-1]; //in training set
      for (i = 0; i < x2.nrows(); i++) {
        if (pnode[i] == pleaf[j]) {
          pnode[i] = 2*j;
          CharacterVector tempMod;
          tempMod.push_back(sv[i]);
          if (!contain(tempModOld, tempMod)) {
            pnode[i] = NA_REAL; // if a new category is observed 
            // in the test set, assign NA 
          } else {
            if (!contain(tempSP, tempMod)) {
              pnode[i] = pnode[i]+1;
            }
          }

        }

    }

    }
    res(_, j) = pnode;
  }
  return res;
    
}
  



