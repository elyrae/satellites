#include "swarm.h"
#include "mathstuff.h"

#include <fstream>
#include <random>
#include <iostream>
#include <ctime>
//#include <omp.h>

using SwarmMethod::Optimizer;
using SwarmMethod::Swarm;

void print_position(const SwarmMethod::SwarmPosition &g, const int iteration, const int elapsed)
{
    std::cout.precision(4);
    std::cout << "Iteration " << iteration+1 << ": [";
    for (size_t i = 0; i < g.pos.size(); i++)
        std::cout << Stuff::rad_to_degrees(g.pos[i]) << ((i < g.pos.size() - 1) ? ", " : "");
    std::cout.precision(6);
    std::cout << "]: " << g.value << "s, " << elapsed << " ms.\n";
    std::cout.flush();
}

int to_ms(const clock_t start, const clock_t end)
{
    return (end - start) / (CLOCKS_PER_SEC / 1000);
}

SwarmMethod::SwarmPosition Optimizer::optimize(const SwarmMethod::TimeFunction F)
{
    SwarmMethod::Swarm swarm(params.swarm_size, reg);

    clock_t start = clock();
    swarm.update_best_position(tg, F, true);
    clock_t end = clock();
    if (params.debug)
        print_position(swarm.best_position(), 0, to_ms(start, end));

    for (size_t iteration = 1; iteration < params.max_iterations; iteration++) {
        swarm.update(params);

        start = clock();
        swarm.update_best_position(tg, F);
        end = clock();
        if (params.debug)
            print_position(swarm.best_position(), iteration, to_ms(start, end));
    }
    return swarm.best_position();
}

// #pragma omp declare reduction(bestPosition: SwarmMethod::SwarmPosition: omp_out = comparePos(omp_out, omp_in))
// #pragma omp parallel for num_threads(2) reduction(bestPosition : p) private(currentF)

// SwarmMethod::SwarmPosition comparePos(const SwarmMethod::SwarmPosition &a, const SwarmMethod::SwarmPosition &b) {
//     return (a.value < b.value) ? a : b;
// }

Swarm::Swarm(const size_t swarm_size, const Region &reg)
{
    srand(time(NULL));
    pos.resize(swarm_size);
    vel.resize(swarm_size);
    for (size_t i = 0; i < pos.size(); i++) {
        pos[i].resize(reg.size());
        vel[i].resize(reg.size());
        
        for (size_t j = 0; j < pos[i].size(); j++) {
            pos[i][j] = Stuff::random(reg[j].first, reg[j].second);
            vel[i][j] = 0.0; // abs(reg[j].first - reg[j].second)*Stuff::random(-1.0, 1.0);
        }
    }
}

void Swarm::update_best_position(Surface::Timegrid &tg, const SwarmMethod::TimeFunction F, const bool rewrite)
{
    SwarmMethod::SwarmPosition candidate;
    candidate.pos = pos[0];
    candidate.value = F(tg, candidate.pos);

    for (size_t i = 1; i < pos.size(); i++) {
        double current_f = F(tg, pos[i]);
        if (current_f < candidate.value) {
            candidate.pos = pos[i];
            candidate.value = current_f;
        }
    }

    if (rewrite)
        best = candidate;
    else if (candidate.value < best.value)
        best = candidate;
}

void Swarm::update(const SwarmMethod::Parameters &params)
{
    for (size_t i = 0; i < pos.size();    i++)
    for (size_t j = 0; j < pos[i].size(); j++) {
        vel[i][j]  = params.omega*vel[i][j] + params.phi*Stuff::random(0.0, 1.0)*(best.pos[j] - pos[i][j]);
        pos[i][j] += vel[i][j];
    }
}
