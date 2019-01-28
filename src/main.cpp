#include <iostream>
#include <math.h>
#include <fstream>
#include <ostream>
#include <vector>

// #include "grid.h"
#include "satellitesurface.h"
#include "orbits.h"
// #include "earth.h"
// #include "messages.h"
#include "mathstuff.h"
#include "examples.h"

// #include "cuda_surface.cuh"

// inline double horizon(const double H, const double alpha)
// {
//     const double delta = Stuff::degreesToRad(10.0); // требуемое возвышение спутника над горизонтом
//     const double alpha_star = asin(cos(delta) / H);

//     return (alpha < alpha_star) ? cos(asin(H*sin(alpha)) - alpha) : sin(delta + alpha_star);
// }

int main() {
    Examples::surface_computation();
    
    return 0;
}
