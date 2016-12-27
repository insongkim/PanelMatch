#include <Rcpp.h>
using namespace Rcpp;




// [[Rcpp::export]]
bool all_sug(LogicalVector x) {
  // Note the use of is_true to return a bool type.
  return is_true(all(x == TRUE));
}


// [[Rcpp::export]]
NumericMatrix rbind_c (NumericMatrix x, NumericMatrix y){
  NumericMatrix out (x.nrow()+y.nrow(), x.ncol());
  for (int r = 0; r < x.nrow(); r ++) {
    out(r, _) = x(r,_);
  }
  for (int s = x.nrow(); s < out.nrow(); s ++){
    out(s, _) = y(s-x.nrow(),_);
  }
  
  return out;
}

// [[Rcpp::export]]
List findDDmatched2(int L, int F, NumericMatrix x1) {
  int nrow1 = x1.nrow();
 
  List out2(nrow1);
  
  for (int i = L; i < (nrow1 - F); i++) {
    if (x1(i,2) == 1) {
      List out(nrow1);
      for (int j = L; j < (nrow1 - F); j++) {
        
        if (x1(j,2) == 0 & x1(j,0) == x1(i,0) & 
            x1(j-L,1) == x1(j+F,1) &
            x1(i-L,1) == x1(i+F,1) &
            x1(j-L,0) == x1(j,0)-L & x1(i-L,0) == x1(i,0)-L &
            x1(j+F,0) == x1(j,0)+F & x1(i+F,0) == x1(i,0)+F
            ) {
          
          NumericMatrix sm = x1( Range(j-L, j-1) , Range(2,2));
          NumericVector c = sm(_,0);
          NumericMatrix sn = x1( Range(i-L, i-1) , Range(2,2));
          NumericVector d = sn(_,0);
          
          if (all_sug(c == d)){
            out[j] = rbind_c(x1(Range(j-L, j+F), _), x1(Range(i-L,i+F),_));
          } else {
            out[j] = R_NilValue;
          }
     
        } else {
          out[j] = R_NilValue;
        }
        
      }
     
      out2[i] = out;
    } else {
      out2[i] = R_NilValue;
    }
    
    
  }
  return out2;
}

