#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix pairMatC(NumericMatrix dat) {
	long nrow = dat.nrow(); long ncol = dat.ncol();
	double newrow = (nrow*(nrow-1))/2; // this is choose(nrow, 2)
	NumericMatrix out(newrow, ncol);
	
	for(long i=0; i < ncol; i++){
		for(long j=0; j < (nrow-1); j++){
			for(long k=j+1; k < nrow; k++){
				long idx = nrow*j - j*(j+1)/2 + (k - j - 1);
				out(idx, i) = dat(j,i) < dat(k,i);
			}
		}	
	}

	/* Each row of out is the pairwise comparison of
	1<2, 1<3,...1<ncol, 2<3, 2<4,...2<ncol,...(ncol-1)<ncol,
	in that order */
	return out;
}