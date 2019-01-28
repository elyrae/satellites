#include "swarm.h"
#include "mathstuff.h"

#include <fstream>
#include <random>
#include <iostream>
#include <ctime>
#include <omp.h>

void print_best_position(const int iteration, const SwarmMethod::SwarmPosition &g, const int elapsed)
{
    std::cout.precision(4);
    std::cout << "Iteration " << iteration+1 << ": [";
    for (size_t i = 0; i < g.pos.size(); i++)
        std::cout << Stuff::rad_to_degrees(g.pos[i]) << ((i < g.pos.size() - 1) ? ", " : "");
    std::cout.precision(6);
    std::cout << "]: " << g.value << "s, " << elapsed << " ms.\n";
    //std::cout << g.second << " ";
    
    //std::cout << g.pos.size() << "\n";
}

// SwarmMethod::SwarmPosition comparePos(const SwarmMethod::SwarmPosition &a, const SwarmMethod::SwarmPosition &b) {
//     return (a.value < b.value) ? a : b;
// }

SwarmMethod::SwarmPosition best_swarm_position(const SwarmMethod::Function F, const SwarmMethod::Swarm &swarm)
{
    SwarmMethod::SwarmPosition p;
    p.pos = swarm[0];
    p.value = F(p.pos);

    double current_f = 0.0;

    // #pragma omp declare reduction(bestPosition: SwarmMethod::SwarmPosition: omp_out = comparePos(omp_out, omp_in))
    // #pragma omp parallel for num_threads(2) reduction(bestPosition : p) private(currentF)
    for (size_t i = 1; i < swarm.size(); i++) {
        current_f = F(swarm[i]);
        if (current_f < p.value) {
            p.pos = swarm[i];
            p.value = current_f;
        }
    }
    return p;
}

void initialize_swarm(const SwarmMethod::Region &reg, SwarmMethod::Swarm &swarm)
{
    for (size_t i = 0; i < swarm.p.size(); i++) {
        swarm.p[i].resize(reg.size());
        swarm.v[i].resize(reg.size());
        
        for (size_t j = 0; j < swarm[i].size(); j++) {
            swarm.p[i][j] = rnd(reg[j].first, reg[j].second);
            swarm.v[i][j] = abs(reg[j].first - reg[j].second)*rnd(-1.0, 1.0);
        }
    }
}

void update_swarm(const SwarmMethod::SwarmPosition &g, const SwarmMethod::Parameters &params, SwarmMethod::Swarm &swarm)
{
    srand(time(NULL));
    for (size_t i = 0; i < swarm.p.size();    i++)
    for (size_t j = 0; j < swarm.p[i].size(); j++) {
        swarm.v[i][j]  = params.omega*swarm.v[i][j] + params.phi*rnd(0.0, 1.0)*(g.pos[j] - swarm.p[i][j]);
        swarm.p[i][j] += swarmVelocities[i][j];
    }
}

SwarmMethod::SwarmPosition SwarmMethod::optimize(const SwarmMethod::Function F, const SwarmMethod::Region &reg,
                                                 const SwarmMethod::Parameters &params,
                                                 const bool debug, const bool writeSwarm)
{
    SwarmMethod::Swarm swarm(params.swarm_size);
    initialize_swarm(reg, swarm);

    clock_t start = clock();
    auto g = best_swarm_position(F, swarm);
    clock_t end = clock();

    if (debug)
        print_bestPosition(0, g, (end - start) / (CLOCKS_PER_SEC / 1000));

    SwarmMethod::SwarmPosition candidate;
    for (size_t iteration = 0; iteration < params.max_iterations; iteration++) {
        update_swarm(g, params, swarm);

        start = clock();
        candidate = best_swarm_position(F, swarm);
        end = clock();

        if (candidate.value < g.value) 
            g = candidate;

        if (debug) {
            print_best_position(iteration, g, (end - start) / (CLOCKS_PER_SEC / 1000));
            std::cout.flush();
        }
    }
    return g;
}
