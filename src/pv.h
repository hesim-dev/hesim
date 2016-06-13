# ifndef PV_H
# define PV_H

unsigned int factorialC(unsigned int n);

double pv1(double t1, double t2, double a, double r);

double exponent_int(double r, unsigned int p, double t, double h);

double poly_int(double r, unsigned int p, double t, double h, std::vector<double> a);

double pv_poly(double r, unsigned int p, double t1,
               double t2, double h, std::vector<double> a);

# endif


