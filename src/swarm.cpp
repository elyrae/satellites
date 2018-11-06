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
    // std::cout.precision(4);
    // std::cout << "Iteration " << iteration+1 << ": [";
    // for (size_t i = 0; i < g.pos.size(); i++)
    //     std::cout << MathStuff::radToDegrees(g.pos[i]) << ((i < g.pos.size() - 1) ? ", " : "");
    // std::cout.precision(6);
    // std::cout << "]: " << g.value << "s, " << elapsed << " ms.\n";
    // std::cout << g.second << " ";
    std::cout << g.pos.size() << "\n";
}

SwarmMethod::SwarmPosition comparePos(const SwarmMethod::SwarmPosition &a, const SwarmMethod::SwarmPosition &b) {
    return (a.value < b.value) ? a : b;
}

// #pragma omp declare reduction(minPosition: SwarmMethod::SwarmPosition: omp_out = comparePos(omp_out, omp_in))

SwarmMethod::SwarmPosition bestSwarmPosition(const Opt::TargetFunction targetF, const SwarmMethod::Swarm& swarm)
{
    SwarmMethod::SwarmPosition p;
    p.pos = swarm[0];
    p.value = targetF(p.pos);

    double currentF = 0.0;

    // #pragma omp parallel for num_threads(2) reduction(minPosition : p) private(currentF)
    for (size_t i = 1; i < swarm.size(); i++) {
        currentF = targetF(swarm[i]);
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

double mt_rand(const std::pair<double, double>& interval)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(interval.first, interval.second);
    return dist(gen);
}

SwarmMethod::SwarmPosition SwarmMethod::optimize(const Opt::TargetFunction targetF, const SwarmMethod::Region& reg,
                                                 const SwarmMethod::Parameters& parameters,
                                                 const bool debugMode, const bool writeSwarm)
{
    std::vector<SwarmMethod::Swarm> iterations;
    SwarmMethod::Swarm swarm(parameters.swarmSize), swarmVelocities(parameters.swarmSize);
    for (size_t i = 0; i < swarm.size(); i++) {
        swarm[i].resize(reg.size());
        swarmVelocities[i].resize(reg.size());
        for (size_t j = 0; j < swarm[i].size(); j++) {
            swarm[i][j] = mt_rand(reg[j]);
            swarmVelocities[i][j] = mt_rand({-abs(reg[j].first - reg[j].second),
                                              abs(reg[j].first - reg[j].second)});
        }
    }

    clock_t start = clock();
    auto g = bestSwarmPosition(targetF, swarm);
    clock_t end = clock();

    if (debugMode)
        printBestPosition(0, g, (end - start) / (CLOCKS_PER_SEC / 1000));
    if (writeSwarm)
        iterations.push_back(swarm);

    int iteration = 0;
    SwarmMethod::SwarmPosition bestPositionCandidate;
    while (iteration < parameters.maxIterations) {
        for (size_t i = 0; i < swarm.size();    i++)
        for (size_t j = 0; j < swarm[i].size(); j++) {
            swarmVelocities[i][j] = parameters.omega*swarmVelocities[i][j]
                                  + parameters.phi*mt_rand({0.0, 1.0})*(g.pos[j] - swarm[i][j]);
            swarm[i][j] += swarmVelocities[i][j];
        }

        start = clock();
        bestPositionCandidate = bestSwarmPosition(targetF, swarm);
        end = clock();

        if (bestPositionCandidate.value < g.value)
            g = bestPositionCandidate;

        if (debugMode) {
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
