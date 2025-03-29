#include <math.h>
#include "../include/solver.h"

double kappa(double u) 
{
    if (u < 0.0 || isnan(u) || isinf(u)) return 0.0;
    //return KAPPA0 * sqrt(u);
    return KAPPA0 * u * u;
}

double heaviside(double x) 
{
    return (x <= 0.2) * BETA;
}

