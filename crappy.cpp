#####

library(RcppArmadillo)
library(inline)
library(genefilter)

code2 <- '
arma::mat datc = Rcpp::as<arma::mat>(dat);
arma::mat modc = Rcpp::as<arma::mat>(mod);
arma::mat mod0c = Rcpp::as<arma::mat>(mod0);
int n = datc.n_rows; int m = datc.n_cols;
int df = modc.n_cols; int df0 = mod0c.n_cols;

arma::vec fstats = arma::zeros<arma::vec>(m);
arma::vec res = arma::zeros<arma::vec>(n);
arma::vec res0 = arma::zeros<arma::vec>(n);

double ss,ss0; 

for(int i=0; i < m; i++){
res = datc.col(i) - modc*arma::solve(modc,datc.col(i));
res0 = datc.col(i) - mod0c*arma::solve(mod0c,datc.col(i));
ss = arma::as_scalar( arma::trans(res)*res0);
ss0 = arma::as_scalar( arma::trans(res0)*res0);
fstats(i) = ((ss0 - ss)/(df-df0))/(ss/(n-df));
}
return Rcpp::wrap(fstats);

'

rcppFstats <- cxxfunction(signature(dat="numeric",mod="numeric",mod0="numeric"),
                          code2,plugin="RcppArmadillo",verbose=TRUE)