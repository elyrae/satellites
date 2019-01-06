#ifndef SWARM_H
#define SWARM_H

// #include "opt.h"
#include <vector>
#include <utility>

namespace SwarmMethod {
    using Function = double(*)(const std::vector<double>&);

    using Point = std::vector<double>;
    using Swarm = std::vector<Point>;
    using Region = std::vector<std::pair<double, double>>;

    // Параметры метода роя частиц
    struct Parameters {
        int swarmSize;
        double omega;
        double phi;

        int maxIterations;
    };

    struct SwarmPosition
    {
        Point pos;
        double value;
    };

    const Parameters defaults = { .swarmSize = 250, .omega = -0.32, .phi = 2.0, .maxIterations = 10};
    SwarmPosition optimize(const Function F, const Region &reg, const Parameters &params = defaults,
                           const bool debug = false, const bool writeSwarm = false);
}

#endif // SWARM_H
