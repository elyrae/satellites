#include "swarm.h"
#include "mathstuff.h"

#include <fstream>
#include <random>
#include <iostream>
#include <ctime>
#include <omp.h>

void outSwarmToFile(const std::vector<SwarmMethod::Swarm>& iterations)
{
    const std::string swarmFile = "outSwarm.txt";

    std::ofstream out(swarmFile);
    out << "\n{";
    for (size_t k = 0; k < iterations.size(); k++) {
        out << "{";
        for (size_t i = 0; i < iterations[k].size(); i++) {
            out << "{";
            for (size_t j = 0; j < iterations[k][i].size(); j++)
                out << iterations[k][i][j] << ((j < iterations[k][i].size() - 1) ? "," : "");
            out << "}" << ((i < iterations[k].size() - 1) ? "," : "");
        }
        out << "}" << ((k < iterations.size() - 1) ? "," : "");
    }
    out << "}";
}

void printBestPosition(const int iteration, const SwarmMethod::SwarmPosition &g, const int elapsed)
{
    std::cout.precision(4);
    std::cout << "Iteration " << iteration+1 << ": [";
    for (size_t i = 0; i < g.pos.size(); i++)
        std::cout << MathStuff::radToDegrees(g.pos[i]) << ((i < g.pos.size() - 1) ? ", " : "");
    std::cout.precision(6);
    std::cout << "]: " << g.value << "s, " << elapsed << " ms.\n";
    //std::cout << g.second << " ";
    
    //std::cout << g.pos.size() << "\n";
}

SwarmMethod::SwarmPosition comparePos(const SwarmMethod::SwarmPosition &a, const SwarmMethod::SwarmPosition &b) {
    return (a.value < b.value) ? a : b;
}

SwarmMethod::SwarmPosition bestSwarmPosition(const SwarmMethod::Function F, const SwarmMethod::Swarm& swarm)
{
    SwarmMethod::SwarmPosition p;
    p.pos = swarm[0];
    p.value = F(p.pos);

    double currentF = 0.0;

    // #pragma omp declare reduction(bestPosition: SwarmMethod::SwarmPosition: omp_out = comparePos(omp_out, omp_in))
    // #pragma omp parallel for num_threads(2) reduction(bestPosition : p) private(currentF)
    for (size_t i = 1; i < swarm.size(); i++) {
        currentF = F(swarm[i]);
        if (currentF < p.value) {
            p.pos = swarm[i];
            p.value = currentF;
        }
    }

    return p;
}



// std::pair<Opt::Point, double> bestSwarmPosition(const Opt::TargetFunction targetF, const SwarmMethod::Swarm& swarm)
// {
//     Opt::Point g = swarm[0];
//     double opt = targetF(g), value = 0.0;

//     //#pragma omp parallel for num_threads(2) reduction(opt : min) 
//     for (size_t i = 1; i < swarm.size(); i++) {
//         value = targetF(swarm[i]);
//         if (value < opt) {
//             opt = value;
//             g = swarm[i];
//         }
//     }
//     return {g, opt};
// }

// double mt_rand(const std::pair<double, double>& interval)
// {
//     static std::random_device rd;
//     static std::mt19937 gen(rd());
//     std::uniform_real_distribution<> dist(interval.first, interval.second);
//     return dist(gen);
// }

double rnd(const double a, const double b)
{
    return a + (b - a)*((double) rand()/RAND_MAX);
}


SwarmMethod::SwarmPosition SwarmMethod::optimize(const SwarmMethod::Function F, const SwarmMethod::Region& reg,
                                                 const SwarmMethod::Parameters& params,
                                                 const bool debug, const bool writeSwarm)
{
    std::vector<SwarmMethod::Swarm> iterations(0);
    SwarmMethod::Swarm swarm(params.swarmSize), swarmVelocities(params.swarmSize);
    srand(time(NULL));

    for (size_t i = 0; i < swarm.size(); i++) {
        swarm[i].resize(reg.size());
        swarmVelocities[i].resize(reg.size());
        for (size_t j = 0; j < swarm[i].size(); j++) {
            swarm[i][j] = rnd(reg[j].first, reg[j].second);
            swarmVelocities[i][j] = abs(reg[j].first - reg[j].second)*rnd(-1.0, 1.0);
        }
    }

    clock_t start = clock();
    auto g = bestSwarmPosition(F, swarm);
    clock_t end = clock();

    if (debug)
        printBestPosition(0, g, (end - start) / (CLOCKS_PER_SEC / 1000));
    if (writeSwarm)
        iterations.push_back(swarm);

    int iteration = 0;
    SwarmMethod::SwarmPosition bestPositionCandidate;
    while (iteration < params.maxIterations) {
        for (size_t i = 0; i < swarm.size();    i++)
        for (size_t j = 0; j < swarm[i].size(); j++) {
            swarmVelocities[i][j] = params.omega*swarmVelocities[i][j] 
                                  + params.phi*rnd(0.0, 1.0)*(g.pos[j] - swarm[i][j]);
            swarm[i][j] += swarmVelocities[i][j];
        }

        start = clock();
        bestPositionCandidate = bestSwarmPosition(F, swarm);
        end = clock();

        if (bestPositionCandidate.value < g.value)
            g = bestPositionCandidate;

        if (debug) {
            printBestPosition(iteration, g, (end - start) / (CLOCKS_PER_SEC / 1000));
            std::cout.flush();
        }
        if (writeSwarm)
            iterations.push_back(swarm);
        iteration++;
    }

    if (writeSwarm)
        outSwarmToFile(iterations);
    return g;
}
