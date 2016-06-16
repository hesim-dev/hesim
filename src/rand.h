# ifndef RAND_H
# define RAND_H

double qgompertzC (double p, double shape, double rate);

double rgompertzC (double shape, double rate);

double rsurv(double location, double anc1, std::string dist);

# endif


