#ifndef SWARM_H
#define SWARM_H

// #include "opt.h"
#include <vector>
#include <utility>

namespace SwarmMethod {
    using Function = double(*)(const std::vector<double>&);

    using Point = std::vector<double>;
    using Region = std::vector<std::pair<double, double>>;

    struct Swarm
    {
        std::vector<Point> p;
        std::vector<Point> v;
    };

    // Параметры метода роя частиц
    struct Parameters {
        int swarm_size;
        double omega;
        double phi;

        int max_iterations;
    };

    struct SwarmPosition
    {
        Point pos;
        double value;
    };

    const Parameters defaults = { .swarm_size = 250, .omega = -0.32, .phi = 2.0, .max_iterations = 10};
    SwarmPosition optimize(const Function F, const Region &reg, const Parameters &params = defaults,
                           const bool debug = false, const bool write_swarm = false);
}

#endif // SWARM_H
