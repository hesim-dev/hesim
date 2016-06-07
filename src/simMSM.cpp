// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>
#include <RcppArmadillo.h>
#include "rand.h"
#include "pv.h"
using namespace Rcpp;

// Patient Class
class Patient {
  private:
    arma::vec location_x_;
    arma::vec trans_;
    double par2_;
    arma::vec cost_x_;

    public:
      Patient(arma::vec lx, arma::vec tr, double par2, arma::vec cx);
      void setTrans(arma::vec tr);
      void setLocationX(arma::vec lx);
      void setCostX(arma::vec cx);
      arma::vec timeJumps(arma::mat lb, std::vector<std::string> d, int nstates);
      double pvCosts(double t1, double t2, int state, double r,
                     arma::vec cost_beta, arma::vec cost_poly_beta,
                     arma::vec cost_poly_deg, arma::vec cost_knots);
};

// Constructor
Patient::Patient(arma::vec lx, arma::vec tr, double par2, arma::vec cx){
  location_x_ = lx;
  trans_ = tr;
  par2_ = par2;
  cost_x_ = cx;
}

// Public setter for private data member trans_
void Patient::setTrans(arma::vec tr){
  trans_ = tr;
}

// Public setter for private data member location_x_
void Patient::setLocationX(arma::vec lx){
  location_x_ = lx;
}

// Public setter for private data member cost_x_
void Patient::setCostX(arma::vec cx){
  cost_x_ = cx;
}

// Times to next states
arma::vec Patient::timeJumps(arma::mat lb, std::vector<std::string> d, int nstates){
  arma::vec time_jumps(nstates);
    for (int j = 0; j < nstates; ++j){
      if (trans_(j) == 0){
        time_jumps(j) = NAN;
      }
    else{
      double location = dot(location_x_, lb.row(trans_(j) - 1));
      time_jumps(j) = rSurv(location, par2_, d[trans_(j)-1]);
    }
  }
  return time_jumps;
}

// Return Polynomial Constants
std::vector<double> retPolyConst(std::vector<double> &a, double xb,
                                  arma::vec coef, int start, int p){
  a.push_back(xb + coef(start));
  for (int j = start + 1; j <= start + p; ++j){
    a.push_back(coef(j));
  }
  return a;
}

// Present Value of Costs
double Patient::pvCosts(double t1, double t2, int state, double r,
                        arma::vec cost_beta, arma::vec cost_polybeta,
                        arma::vec cost_polydeg, arma::vec cost_knots){

  int n = cost_polydeg.n_elem;
  cost_knots = cost_knots + t1;
  cost_knots(n) = t2;
  cost_knots = sort(cost_knots);
  double xb = dot(cost_beta, cost_x_);
  std::vector<double> a;
  double pv = 0;
  int coef_start = 0;
  for (int i = 0; i < n; ++i){
    if (cost_knots(i + 1) > t2){
      break;
    }
    a = retPolyConst(a, xb, cost_polybeta, coef_start, cost_polydeg(i));
    pv += pvPoly(r, cost_polydeg(i), cost_knots(i), cost_knots(i + 1),
                   cost_knots(i), a);
    a.clear();
    coef_start = coef_start + cost_polydeg(i) + 1;
  }
  return pv;
}

// [[Rcpp::export]]
double Fun(arma::vec lx, arma::vec tr, double p2, arma::vec cx,
           double t1, double t2, int state, double r, arma::vec costBeta,
           arma::vec costPolyBeta, arma::vec costPolyDegree,
           std::vector<double> cost_knots){
  Patient p = Patient(lx, tr, p2, cx);
  return p.pvCosts(t1, t2, state, r, costBeta, costPolyBeta,
                   costPolyDegree, cost_knots);
}

// [[Rcpp::export]]
arma::vec myFun(arma::mat lb, arma::vec lx, int nstates, arma::vec cx,
             arma::vec trans, double par2, std::vector<std::string> d){
   Patient p1 = Patient(lx, trans, par2, cx);
  return p1.timeJumps(lb, d, nstates);
}

// [[Rcpp::export]]
bool notinVec(int item, std::vector<int> v){
  return std::find(v.begin(), v.end(), item) == v.end();
}

// [[Rcpp::export]]
double nextTime(double current_time, double time_jump, double maxt){
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
arma::mat locationBetaMat(int index, arma::cube lb, int ntrans, int k){
  arma::mat beta(ntrans, k);
  for (int i = 0; i < ntrans; ++i){
    beta.row(i) = lb.slice(i).row(index);
  }
  return beta;
}

// [[Rcpp::export]]
arma::vec updateAge(arma::vec x, double age, int col){
  x(col) = age;
  return x;
}

// [[Rcpp::export]]
List simMSM(arma::cube location_beta, arma::mat location_x,
            std::vector<std::string> dist, arma::mat tmat,
            double par2, std::vector<int> absorbing, int maxt,
            int agecol, std::vector<double> qol, double r,
            arma::mat cost_x, arma::mat cost_beta,
            arma::mat cost_polybeta, arma::mat cost_polydeg,
            arma::mat cost_knots) {
  // Initialize
  int N = location_x.n_rows;
  int ntrans = location_beta.n_slices;
  int k_lb = location_beta.slice(0).n_cols;
  int nsims = location_beta.slice(0).n_rows;
  int nstates = tmat.n_cols;

  // Storage vectors
  std::vector<int> id;
  std::vector<int> sim;
  std::vector<int> state;
  std::vector<double> time;
  std::vector<double> qalys;
  std::vector<double> costs;

  // Begin simulation
  int counter = 0;

  // Simulate for parameter draws 1 to nsims
  for (int s = 0; s < nsims; ++s){
    arma::mat location_beta_s = locationBetaMat(s, location_beta, ntrans, k_lb);

    // Simulate for individuals 1 to N
    for (int i = 0; i < N; ++i){
      // Values for time = 0
      id.push_back(i);
      sim.push_back(s);
      state.push_back(0);
      time.push_back(0.0);
      qalys.push_back(0.0);
      costs.push_back(0.0);

      // Simulate for patient i
      arma::vec location_xi = location_x.row(i).t();
      arma::vec cost_xi = cost_x.row(i).t();
      Patient p = Patient(location_xi, tmat.row(0).t(), par2, cost_xi);
      while (notinVec(state[counter], absorbing) && time[counter] < maxt){
        // Current iteration
        arma::vec time_jumps = p.timeJumps(location_beta_s, dist, nstates);
        arma::uword next_state;
        double time_jump = time_jumps.min(next_state);
        time.push_back(nextTime(time[counter], time_jump, maxt));
        state.push_back(next_state);
        id.push_back(i);
        sim.push_back(s);
        qalys.push_back(pv1(time[counter], time[counter + 1],
                        qol[state[counter]], r));
        costs.push_back(p.pvCosts(time[counter], time[counter + 1], state[counter],
                                  r, cost_beta.row(state[counter]).t(),
                                  cost_polybeta.row(state[counter]).t(),
                                  cost_polydeg.row(state[counter]).t(),
                                  cost_knots.row(state[counter]).t()
                                  ));

        // Move to next iteration and update patient
        ++ counter;
        double age = location_xi(agecol) + time_jump;
        p.setTrans(tmat.row(state[counter]).t());
        p.setLocationX(updateAge(location_xi, age, agecol));
      }
      ++ counter;
    }
  }
  return List::create(id, sim, state, time, qalys, costs);
}

// R code for testing
/*** R
x <- matrix(1, 2, 2)
tmat <- matrix(c(0, 1, 2, 3, 0, 4, 0, 0, 0), 3, 3, byrow = TRUE)
par2 <- 1
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
simMSM(beta, x, dist, tmat, par2, abosrbing, maxt, agecol, qol, r, x, beta[,,1],
       cost.poly.beta.mat, cost.poly.degree.mat, cost.knots.mat)
*/
