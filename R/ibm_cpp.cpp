#include <Rcpp.h>
#include <vector>
using namespace Rcpp;


// Top level function to iterate the IBM n_iterations number of times.
// C++ equivalents of the vital rate functions are defined below as well.
// I think you can actually call R functions from C++, but I don't really
// feel like figuring it out right now.
// Obviously, this is not ready to go just yet, but would like to have it up and
// running soon!

//[[Rcpp::export]]
Rcpp::DataFrame iterate_car_ibm(
  std::vector<int> id,
  std::vector<double> z,
  std::vector<int> age,
  Rcpp::List parameters,
  int yr,
  int n_iterations
) {

  // create some vectors to hold data

  int len = id.size();
  std::vector<int> repr(len);
  std::vector<int> seeds(len);
  std::vector<int> surv(len);
  std::vector<double> z1(len);
  std::vector<bool> alive(len);

  Rcpp::DataFrame temp = Rcpp::DataFrame::create(
    Named("id") = id,
    Named("z")  = z,
    Named("repr") = repr,
    Named("seeds") = seeds,
    Named("surv") = surv,
    Named("z1") = z1,
    Named("age") = age,
    Named("alive") = alive,
    Named("yr") = yr
  );

  // Next, move on to iterating the model

  for(int i = 0; i < n_iterations, i++) {

    p_b = p_bz(z, parameters["flow_int"], parameters["flow_z"]);

    repr = Rcpp::rbinom(n    = pop_size,
                        prob = p_b,
                        size = 1);
    int repr_n = sum(repr);

    seeds = Rcpp::rpois()

    if(i == 0) {

      Rcpp::DataFrame out = temp;

    } else {

      Rcpp::DataFrame::rbind() // come up with something here. Maybe switch to list?
    }

  }


}

//[[Rcpp::export]]
std::vector<double> s_z(
  double z,
  double intercept,
  double slope
) {

  double lin_p = intercept + slope * z;
  double p     = 1 / (1 + Rcpp::exp( -lin_p ));
  return p;

}

//[[Rcpp::export]]
std::vector<double> g_z1z(
  double z1,
  double z,
  double intercept,
  double slope,
  double st_dev
) {

  double mu = intercept + slope * z;
  double pr = Rcpp::dnorm(z1, mean = mu, sd = st_dev);
  return pr;
}

//[[Rcpp::export]]
std::vector<double> p_bz(
  double z,
  double intercept,
  double slope
) {

  double lin_p = intercept + slope * z;
  double p     = 1 / (1 + Rcpp::exp(-lin_p));
  return p;
}

//[[Rcpp::export]]
std::vector<double> b_z(
  double z,
  double intercept,
  double slope
) {

  double lin_n = intercept + slope * z;
  double out   = Rcpp::exp(lin_n);
  return out;
}

//[[Rcpp::export]]
std::vector<double> c_0z1(
  double z1,
  double intercept,
  double st_dev
) {

  double out = Rcpp::dnorm(z1, mean = intercept, sd = st_dev);
  return out;
}

