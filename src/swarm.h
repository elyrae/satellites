#ifndef SWARM_H
#define SWARM_H

#include "opt.h"
#include <utility>

namespace SwarmMethod {
    using Swarm = std::vector<Opt::Point>;
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
        Opt::Point pos;
        double value;
    };

    const Parameters defaults = { .swarmSize = 250, .omega = -0.32, .phi = 2.0, .maxIterations = 10};
    SwarmPosition optimize(const Opt::TargetFunction targetF, const Region& bound, const Parameters& parameters = defaults,
                           const bool debugMode = false, const bool writeSwarm = false);
}

#endif // SWARM_H
