#include "mathstuff.h"
#include <cmath>

double MathStuff::degreesToRad(const double x)
{
    return x * M_PI / 180.0;
}

double MathStuff::radToDegrees(const double x)
{
    return x * 180.0 / M_PI;
}
