#include "mathstuff.h"

#include <cmath>
#include <stdlib.h>

double Stuff::degrees_to_rad(const double x)
{
    return x * M_PI / 180.0;
}

double Stuff::rad_to_degrees(const double x)
{
    return x * 180.0 / M_PI;
}

double Stuff::random(const double a, const double b)
{
    return a + (b - a)*((double) rand()/RAND_MAX);
}