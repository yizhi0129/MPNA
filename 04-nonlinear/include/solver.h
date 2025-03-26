#ifndef SOLVER_H
#define SOLVER_H

#define KAPPA0 0.01
#define SIGMA 0.1
#define BETA 1.0

#define GAMMA 0.1
//#define GAMMA 1.0
//#define GAMMA 10.0

//#define SIGMA 1.0;
//#define BETA 300.0;

double kappa(double u); 
double heaviside(double x); 
void solve_linearized_implicit(int N, int max_steps);
void solve_newton_method(int N, double tol, int max_iter);

#endif
