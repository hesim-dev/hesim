# ifndef PV_H
# define PV_H

unsigned int factorialC(unsigned int n);

double pv1(double t1, double t2, double a, double r);

double exponentInt(double r, unsigned int p, double t, double h);

double polyInt(double r, unsigned int p, double t, double h, std::vector<double> a);

double pvPoly(double r, unsigned int p, double t1,
               double t2, double h, std::vector<double> a);

# endif


