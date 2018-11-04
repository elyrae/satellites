#include "swarm.h"
#include "mathstuff.h"

#include <fstream>
#include <random>
#include <iostream>
#include <ctime>

void outSwarmToFile(const std::vector<ParticleSwarmMethod::Swarm>& iterations)
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

void printBestPosition(const int iteration, const std::pair<Opt::Point, double>& g, const int elapsed)
{
    std::cout.precision(4);
    std::cout << "Iteration " << iteration+1 << ": [";
    for (size_t i = 0; i < g.first.size(); i++)
        std::cout << MathStuff::radToDegrees(g.first[i]) << ((i < g.first.size() - 1) ? ", " : "");
    std::cout.precision(6);
    std::cout << "]: " << g.second << "s, " << elapsed << " ms.\n";
    // std::cout << g.second << " ";
}

//void printVariance(const ParticleSwarmMethod::Swarm& g, const int elapsed)
//{

//    QTextStream out(stdout);
//    Opt::Point mean(g[0].size()), var(g[0].size());
//    out << "[";
//    for (size_t i = 0; i < g.first.size(); i++)
//        out << MathStuff::radToDegrees(g.first[i]) << ((i < g.first.size() - 1) ? ", " : "");
//    out << "]: " << g.second << "s, " << elapsed << "ms.";
//}

std::pair<Opt::Point, double> bestSwarmPosition(const Opt::TargetFunction targetF,
                                                const ParticleSwarmMethod::Swarm& swarm,
                                                const Opt::SearchType searchType)
{
    Opt::Point g = swarm[0];
    double opt = targetF(g), value = 0.0;

    //#pragma omp parallel for num_threads(2)
    for (size_t i = 1; i < swarm.size(); i++) {
        value = targetF(swarm[i]);
        if ((searchType == Opt::SearchType::SearchMaximum) ? (value > opt) : (opt > value)) {
            opt = value;
            g = swarm[i];
        }
    }
    return {g, opt};
}

double mt_rand(const std::pair<double, double>& interval)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(interval.first, interval.second);
    return dist(gen);
}

std::pair<Opt::Point, double> ParticleSwarmMethod::optimize(const Opt::TargetFunction targetF,
                                                            const ParticleSwarmMethod::Region& reg,
                                                            const Opt::SearchType searchType,
                                                            const ParticleSwarmMethod::Parameters& parameters,
                                                            const bool debugMode, const bool writeSwarm)
{
    std::vector<ParticleSwarmMethod::Swarm> iterations;
    ParticleSwarmMethod::Swarm swarm(parameters.swarmSize), swarmVelocities(parameters.swarmSize);
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
    // auto start = std::chrono::system_clock::now();
    auto g = bestSwarmPosition(targetF, swarm, searchType);
    clock_t end = clock();
    // auto end = std::chrono::system_clock::now();
    // auto elapsed = end - start; elapsed.count()

    if (debugMode)
        printBestPosition(0, g, (end - start) / (CLOCKS_PER_SEC / 1000));
    if (writeSwarm)
        iterations.push_back(swarm);

    int iteration = 0;
    std::pair<Opt::Point, double> bestPositionCandidate;
    while (iteration < parameters.maxIterations) {
        // if (debugMode) {
        //     std::cout << "Iteration " << iteration+1 << ": ";
        //     std::cout.flush();
        // }

        for (size_t i = 0; i < swarm.size(); i++)
            for (size_t j = 0; j < swarm[i].size(); j++) {
                swarmVelocities[i][j] = parameters.omega*swarmVelocities[i][j]
                                      + parameters.phi*mt_rand({0.0, 1.0})*(g.first[j] - swarm[i][j]);
                swarm[i][j] += swarmVelocities[i][j];
            }

        start = clock();
        bestPositionCandidate = bestSwarmPosition(targetF, swarm, searchType);
        end = clock();

        if ((searchType == Opt::SearchType::SearchMaximum) ? (bestPositionCandidate.second > g.second)
                                                           : (g.second > bestPositionCandidate.second))
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
