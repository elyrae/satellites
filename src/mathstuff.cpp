#include "mathstuff.h"

double MathStuff::degreesToRad(const double x)
{
    return x * M_PI / 180.0;
}

double MathStuff::radToDegrees(const double x)
{
    return x * 180.0 / M_PI;
}
