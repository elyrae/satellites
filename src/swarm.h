#ifndef SWARM_H
#define SWARM_H

#include <vector>
#include <utility>

#include "satellitesurface.h"

namespace SwarmMethod {
    using TimeFunction = double (*)(Surface::Timegrid &, const std::vector<double> &);

    using Point = std::vector<double>;
    using Region = std::vector<std::pair<double, double>>;

    struct SwarmPosition
    {
        Point pos;
        double value;
    };

    // Параметры метода роя частиц
    struct Parameters {
        int swarm_size;
        double omega;
        double phi;

        size_t max_iterations;

        bool debug;
    };

    const Parameters defaults = { .swarm_size = 250, .omega = -0.32, .phi = 2.0, .max_iterations = 10, .debug = true};

    class Optimizer {
        Parameters params;
        Region reg;
        Surface::Timegrid tg;
    public:
        Optimizer(const Surface::Timegrid &_tg, const Region &_reg, const Parameters &_par = defaults) : params(_par),
            reg(_reg), tg(_tg) {}
        SwarmPosition optimize(const TimeFunction F);
    };

    class Swarm {
        std::vector<Point> pos;
        std::vector<Point> vel;

        SwarmPosition best;
    public:
        Swarm(const size_t swarm_size, const Region &reg);
        void update_best_position(Surface::Timegrid &tg, const TimeFunction F);
        void update(const SwarmMethod::Parameters &params);

        SwarmPosition best_position() const { return best; }
    };
}

#endif // SWARM_H
