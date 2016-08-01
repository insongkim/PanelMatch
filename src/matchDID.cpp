# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


//  Test Rcpp code

// [[Rcpp::export()]]
int sumCpp(Rcpp::IntegerVector x) {
  int n = x.size();
  int res = 0;
  for (int i = 0; i < n; i++){
    res += x[i];
  }
  return res;
}


// [[Rcpp::export()]]
NumericMatrix FindMatches(IntegerVector unitIdx,
			  IntegerVector timeIdx,
			  IntegerVector treat) {


  IntegerVector units = unique(unitIdx);
  IntegerVector times = unique(timeIdx);
  int lenUnit = units.size();
  int lenTime = times.size();
  int lenData = treat.size();

  arma::mat Wdid(lenTime, lenUnit);
  Wdid.zeros();
  
  for (int t = 1; t < lenTime; t++) {
    for (int i = 0; i < lenUnit; i++) {

      // checking whether treatment status has changed t vs t-1
      int unit;
      int time;
      int treatT;
      int treatT1;

      int istreatT = 0;
      for (int k = 0; k < lenData; k++){
        unit = unitIdx[k];
	time = timeIdx[k];
      	if( unit == (i+1) && time == (t+1)){
      	  treatT = treat[k];
      	  istreatT = 1;
      	  break;
      	}
      }

      int istreatT1 = 0;
      for (int k = 0; k < lenData; k++){
	unit = unitIdx[k];
	time = timeIdx[k];
      	if( unit == (i+1) && time == (t)){
      	  treatT1 = treat[k];
      	  istreatT1 = 1;
      	  break;
      	}
      }
      
      if(istreatT == 1 && istreatT1){ // when both treatment variable for time t and t-1 exists
      	if(treatT == 1 && treatT1 == 0){ // if treated in time t but control in t-1

	  // Now find matched unit indices for which treat stays in control in t and t-1
	  IntegerVector matchedJ;
	  int unitj;
	  int timej;
	  int unitj2;
	  int timej2;
	  int treatTj;
	  int treatT1j;
	  arma::mat W(lenTime, lenUnit);
	  W.zeros();


	  // finding indices for matched units
	  for (int k = 0; k < lenData; k++){
	    unitj = unitIdx[k];
	    timej = timeIdx[k];
	    if( unitj != (i+1) && timej == (t+1)){
	      treatTj = treat[k];

	      for (int l = 0; l < lenData; l++){
		unitj2 = unitIdx[l];
		timej2 = timeIdx[l];
		if( unitj2 == unitj && timej2 == (t)){
		  treatT1j = treat[l];
		  if(treatTj == 0 && treatT1j == 0){
		    matchedJ.insert(0,unitj);
		    break;
		  }
		}
	      }
	    }
	  }

	  double numMatches = matchedJ.size();
	  // Rcout << "Number of Matched Units: " << numMatches << std::endl;
	  if(numMatches > 0){
	    double vit = 1/numMatches;

	    W(t,i) = 1;
	    W(t-1,i) = 1;

	    for(int m = 0; m < numMatches; m++){
	      int jidx = matchedJ[m]-1;
	      W(t,jidx) = vit;
	      W(t-1,jidx) = -vit;
	    }
	  }
	  
	  Wdid = Wdid + W;

      	  // Rcout << "Unit:" << unit << "  Time:" << time+1 << " matches are: " << matchedJ << std::endl;
      	}
      }
      
    }
  }
  return(wrap(Wdid));
}
