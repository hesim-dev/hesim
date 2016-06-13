// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>
#include <RcppArmadillo.h>
#include "rand.h"
#include "pv.h"
using namespace Rcpp;

// Random times to next states
arma::vec rtjumps(arma::mat lb, arma::vec loc_x, arma::vec par2,
                    arma::vec trans, std::vector<std::string> d,
                    int nstates){
  arma::vec time_jumps(nstates);
  for (int j = 0; j < nstates; ++j){
    if (trans(j) == 0){
      time_jumps(j) = NAN;
    }
    else{
      double loc = dot(loc_x, lb.row(trans(j) - 1));
      time_jumps(j) = rsurv(loc, par2(trans(j) - 1), d[trans(j) - 1]);
    }
  }
  return time_jumps;
}

// Return Polynomial Constants
std::vector<double> poly_const(std::vector<double> &a, double xb,
                                  arma::vec coef, int start, int p){
  a.push_back(xb + coef(start));
  for (int j = start + 1; j <= start + p; ++j){
    a.push_back(coef(j));
  }
  return a;
}

// Calculate Present Value from t1 to t2 With Piecewise Polynomials
// [[Rcpp::export]]
double pv_splines(double t1, double t2, int state, double r,
               arma::vec x, arma::vec beta, arma::vec poly_beta,
               arma::vec poly_deg, arma::vec knots){

  int n = poly_deg.n_elem;
  knots = knots + t1;
  knots(n) = t2;
  knots = sort(knots);
  double xb = dot(beta, x);
  std::vector<double> a;
  double pv = 0;
  int coef_start = 0;
  for (int i = 0; i < n; ++i){
    if (knots(i + 1) > t2){
      break;
    }
    a = poly_const(a, xb, poly_beta, coef_start, poly_deg(i));
    pv += pv_poly(r, poly_deg(i), knots(i), knots(i + 1),
                 knots(i), a);
    a.clear();
    coef_start = coef_start + poly_deg(i) + 1;
  }
  return pv;
}

// [[Rcpp::export]]
bool notinvec(int item, std::vector<int> v){
  return std::find(v.begin(), v.end(), item) == v.end();
}

// [[Rcpp::export]]
double nxttime(double current_time, double time_jump, double maxt){
  double next_time;
  if (current_time + time_jump > maxt){
    next_time = maxt;
  }
  else {
    next_time = current_time + time_jump;
  }
  return next_time;
}

// [[Rcpp::export]]
arma::mat matS(int index, arma::cube cube, int ntrans, int k){
  arma::mat par(ntrans, k);
  for (int i = 0; i < ntrans; ++i){
    par.row(i) = cube.slice(i).row(index);
  }
  return par;
}

// [[Rcpp::export]]
arma::vec updateAge(arma::vec x, double age, int col){
  x(col) = age;
  return x;
}

//' @export
// [[Rcpp::export]]
List sim_msmC(arma::cube loc_beta, arma::mat loc_x,
            std::vector<std::string> dist, arma::mat tmat,
            arma::vec par2, std::vector<int> absorbing,
            int maxt, int agecol = -1) {
  // Initialize
  int N = loc_x.n_rows;
  int ntrans = loc_beta.n_slices;
  int k_lb = loc_beta.slice(0).n_cols;
  int nsims = loc_beta.slice(0).n_rows;
  int nstates = tmat.n_cols;

  // Storage vectors
  std::vector<int> id;
  std::vector<int> sim;
  std::vector<int> state;
  std::vector<double> time;

  // Begin simulation
  int counter = 0;

  // Simulate for parameter draws 1 to nsims
  for (int s = 0; s < nsims; ++s){
    arma::mat loc_beta_s = matS(s, loc_beta, ntrans, k_lb);

    // Simulate for individuals 1 to N
    for (int i = 0; i < N; ++i){
      // Values for time = 0
      id.push_back(i);
      sim.push_back(s);
      state.push_back(0);
      time.push_back(0.0);

      // Simulate for patient i
      arma::vec loc_xi = loc_x.row(i).t();
      while (notinvec(state[counter], absorbing) && time[counter] < maxt){
        // Current iteration
        arma::vec time_jumps = rtjumps(loc_beta_s, loc_xi, par2,
                                         tmat.row(state[counter]).t(),
                                         dist, nstates);
        arma::uword next_state;
        double time_jump = time_jumps.min(next_state);
        time.push_back(nxttime(time[counter], time_jump, maxt));
        state.push_back(next_state);
        id.push_back(i);
        sim.push_back(s);

        // Move to next iteration and update patient
        ++ counter;
        if (agecol >= 0){
          double age = loc_xi(agecol) + time_jump;
          loc_xi = updateAge(loc_xi, age, agecol);
        }
      }
      ++ counter;
    }
  }
  return List::create(id, sim, state, time);
}

//' @export
// [[Rcpp::export]]
arma::vec sim_msm_pvC(arma::vec id, arma::vec sim, arma::vec state, arma::vec time,
                   double r, arma::mat x, int agecol,
                   arma::mat beta, arma::mat poly_beta, arma::mat poly_deg,
                   arma::mat knots) {
  // Initialize/store
  int N = time.n_elem;
  arma::vec pv(N);

  // Loop
  for (int i = 0; i < N; ++i){
    if (time(i) == 0){
      pv(i) = 0;
    }
    else{
      pv(i) = pv_splines(time(i - 1), time(i), state(i - 1), r,
         x.row(id(i)).t(), beta.row(state(i - 1)).t(),
         poly_beta.row(state(i - 1)).t(),
         poly_deg.row(state(i - 1)).t(), knots.row(state(i - 1)).t());
    }
  }
  return pv;
}

// R code for testing
/*** R
x <- matrix(1, 2, 2)
tmat <- matrix(c(0, 1, 2, 3, 0, 4, 0, 0, 0), 3, 3, byrow = TRUE)
par2 <- rep(1, 4)
beta <- log(array(seq(1, 16), c(2, 2, 4)))

# Costs
cost.poly.beta <- c(1, 2, 3)
cost.poly.degree <- c(1, 0)
cost.knots <- c(0, 3, 0)
t1 <- 5
t2 <- 10
xb <- x[, 1] %*% beta[1, , 1]
Fun(x[, 1], tmat[1, ], par2, x[, 1], t1, t2, 0, .03, beta[1, , 1],
    cost.poly.beta, cost.poly.degree, cost.knots)
pvPoly(.03, cost.poly.degree[1], cost.knots[1] + t1,
       cost.knots[2] + t1, cost.knots[1] + t1, c(xb + 1, 2))
pvPoly(.03, cost.poly.degree[2], cost.knots[2] + t1,
       t2, cost.knots[2] + t1, c(xb + 3))

# Simulation
dist <- rep("weibull", 4)
abosrbing <- 2
maxt <- 50
agecol <- 1
qol <- c(1, 1, 1)
r <- 0
cost.poly.beta.mat <- matrix(rep(cost.poly.beta, 3), 3, 3, byrow = T)
cost.poly.degree.mat <- matrix(rep(cost.poly.degree, 3), 3, 2, byrow = T)
cost.knots.mat <- matrix(rep(c(0, 3, 0), 3), 3, 3, byrow = T)
# simMsmC(beta, x, dist, tmat, par2, abosrbing, maxt, agecol, qol, r, x, agecol, beta[,,1],
#        cost.poly.beta.mat, cost.poly.d,
#        egree.mat, cost.knots.mat)
sim <- sim_msmC(beta, x, dist, tmat, par2, abosrbing, maxt, agecol)
sim_msm_pvC(sim[[1]], sim[[2]], sim[[3]], sim[[4]],
         .03, x, 1, rbind(beta[,,1], beta[1, , 1]), cost.poly.beta.mat,
             cost.poly.degree.mat, cost.knots.mat)
*/
