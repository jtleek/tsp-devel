#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
NumericVector fstatsC(NumericMatrix pairMat, NumericMatrix mod, NumericMatrix mod0, NumericVector outcome) {
	int nrow = pairMat.nrow(); int ncol = pairMat.ncol();
	double ss, ss0;
	int df0 = mod0.ncol(); int df = df0+1; // alternative model always has +1 column
	
	arma::mat modc(mod.begin(), mod.nrow(), mod.ncol(), false);
	arma::mat mod0c(mod0.begin(), mod0.nrow(), mod0.ncol(), false);
	arma::colvec outcomec(outcome.begin(), outcome.size(), false);
	arma::vec fstats(nrow);
	arma::vec res = arma::zeros<arma::vec>(outcome.size());
	arma::vec res0 = arma::zeros<arma::vec>(outcome.size());
    
	try{ 
		res0 = outcomec - mod0c*arma::solve(mod0c, outcomec); // The residual for the null model remains the same
		ss0 = arma::as_scalar(arma::trans(res0)*res0);
	} catch(std::exception &ex) {
		stop("WTF???");
	}
	for(int i=0; i < nrow; i++){ // loop through all rows in pairMat
		//modc.insert_cols(mod.ncol(), pairMat(i,_)); can try this later, it's not working
		for(int j=0; j < pairMat.ncol(); j++){ // this loop is for copying the ith row of pairMat into the last column of modc
			modc(j,modc.n_cols-1) = pairMat(i,j); // Here we add the current row to the model
		}
		try{
			res = outcomec - modc*arma::solve(modc, outcomec); // Calculate the residual
			ss = arma::as_scalar(arma::trans(res)*res);
			fstats(i) = ((ss0 - ss)/(df-df0))/(ss/(ncol-df));
		} catch(std::exception &ex) {
			fstats(i) = NA_REAL;
		}
	}
	return Rcpp::wrap(fstats);
}

/* R

pairMat <- matrix(rbinom(100*20,size=1,prob=0.5),nrow=100)
mod0 <- cbind(rep(1,20))
mod <- cbind(mod0,rep(0,20))
outcome <- rnorm(20)

ff = fstatsC(pairMat,mod,mod0,outcome)
*/
