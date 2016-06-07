# ifndef RAND_H
# define RAND_H

double qGompertz (double p, double shape, double rate);

double rGompertz (double shape, double rate);

double rSurv(double location, double par2, std::string dist);

# endif


