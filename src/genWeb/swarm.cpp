#include "swarm.h"
#include "mathstuff.h"

#include <QVector>
#include <QTextStream>
#include <QFile>
#include <random>
#include <QTime>
#include <algorithm>

void outSwarmToFile(const std::vector<ParticleSwarmMethod::Swarm>& iterations)
{
    QFile file("outSwarm.txt");
    file.open(QIODevice::WriteOnly | QIODevice::Text);

    QTextStream out(&file);
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
    file.close();
}

void printBestPosition(const std::pair<Opt::Point, double>& g, const int elapsed)
{
    QTextStream out(stdout);
    out.setRealNumberPrecision(4);
    out << "[";
    for (size_t i = 0; i < g.first.size(); i++)
        out << MathStuff::radToDegrees(g.first[i]) << ((i < g.first.size() - 1) ? ", " : "");
    out.setRealNumberPrecision(6);
    out << "]: " << g.second << "s, " << elapsed << "ms.\n";
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
                                                            const ParticleSwarmMethod::Region& bound,
                                                            const Opt::SearchType searchType,
                                                            const ParticleSwarmMethod::Parameters& parameters,
                                                            const bool debugMode, const bool writeSwarm)
{
    QTextStream out(stdout);
    QTime timer;

    std::vector<ParticleSwarmMethod::Swarm> iterations;
    ParticleSwarmMethod::Swarm swarm(parameters.S), swarmVelocities(parameters.S);
    for (size_t i = 0; i < swarm.size(); i++) {
        swarm[i].resize(bound.size());
        swarmVelocities[i].resize(bound.size());
        for (size_t j = 0; j < swarm[i].size(); j++) {
            swarm[i][j] = mt_rand(bound[j]);
            swarmVelocities[i][j] = mt_rand({-abs(bound[j].first - bound[j].second),
                                              abs(bound[j].first - bound[j].second)});
        }
        //std::sort(swarm[i].begin(), swarm[i].begin() + 8);
    }

    if (debugMode)
        timer.start();
    auto g = bestSwarmPosition(targetF, swarm, searchType);
    if (debugMode) {
        printBestPosition(g, timer.elapsed());
        timer.restart();
    }
    if (writeSwarm)
        iterations.push_back(swarm);

    int iteration = 0;
    std::pair<Opt::Point, double> bestPositionCandidate;
    while (iteration < parameters.maxIterations) {
        if (debugMode) {
            out << "Iteration " << iteration+1 << ": ";
            out.flush();
        }
        for (size_t i = 0; i < swarm.size(); i++)
            for (size_t j = 0; j < swarm[i].size(); j++) {
                swarmVelocities[i][j] = parameters.omega*swarmVelocities[i][j]
                                      + parameters.phi*mt_rand({0.0, 1.0})*(g.first[j] - swarm[i][j]);
                swarm[i][j] = swarm[i][j] + swarmVelocities[i][j];
            }

        bestPositionCandidate = bestSwarmPosition(targetF, swarm, searchType);
        if ((searchType == Opt::SearchType::SearchMaximum) ? (bestPositionCandidate.second > g.second)
                                                           : (g.second > bestPositionCandidate.second))
            g = bestPositionCandidate;
        if (debugMode) {
            printBestPosition(g, timer.elapsed());
            out.flush();
            timer.restart();
        }
        if (writeSwarm)
            iterations.push_back(swarm);
        iteration++;
    }
    if (writeSwarm)
        outSwarmToFile(iterations);
    return g;
}
