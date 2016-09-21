#include <Rcpp.h>
using namespace Rcpp;



double q(IntegerVector n_drawn, IntegerVector K_marked, IntegerVector d_unmarked, int N_min,
         int b, double lambda) {
  int N = N_min + b;
  return sum(log(N - n_drawn) - log(N) + log(N-K_marked) - log(N-K_marked - d_unmarked)) - lambda;
}


int rcpp_hgrecap_metropolis_step(int b, NumericVector& Q) {


  double a = R::runif(-1,1);

  int w = copysign(1.0, a);

  if(b == 0 && w == -1) {
    return 0;
  }

  a = std::abs(a);

  double log_P = (w == -1) ? -Q(b-1) : Q(b);

  if(log(a) < log_P) {
    b = b + w;
  }

  return b;

}

// [[Rcpp::export]]
List rcpp_hgrecap_metropolis(IntegerVector n_drawn, IntegerVector K_marked, IntegerVector d_unmarked, int N_min,
                             List prior, List control) {

  int r = n_drawn.length();
  int b = 1;

  // MCMC control constants
  int       M = as<int>(control["mcmc.iterations"]);
  int    thin = as<int>(control["mcmc.thin"]);
  int  burnin = as<int>(control["mcmc.burnin"]);
  int verbose = as<int>(control["verbose"]);

  //Precalculate constants
  double lambda = as<double>(prior["lambda"]);

  NumericVector Q(64);
  for(int i = 0; i < Q.length(); i++)
    Q(i) = q(n_drawn, K_marked, d_unmarked, N_min, i+1, lambda);

  // Allocate output
  IntegerVector outv(M/thin);


  //Burnin
  for(int i = thin; i < burnin; i++) {
    b = rcpp_hgrecap_metropolis_step(b, Q);
    if(b >= Q.length()) {
      Q.push_back( q(n_drawn, K_marked, d_unmarked, N_min, b+1, lambda) );
    }
  }


  //Main loop
  for(int i = 0; i < M/thin; i++) {

    // Run thinned iterations inplace
    for(int j = 0; j < thin; j++) {
      b = rcpp_hgrecap_metropolis_step(b, Q);
      if(b >= Q.length()) {
        Q.push_back(q(n_drawn, K_marked, d_unmarked, N_min, b+1, lambda) );
      }
    }

    if(verbose) {
      ::Rprintf("Iteration %4d: %4d pos\n", i, b);
    }

    // save out
    outv[i] = N_min + b;
  }


  return List::create(outv, Q) ;
}

