#include "examples.h"
#include "satellitesurface.h"
#include "mathstuff.h"
#include "swarm.h"
#include "grid.h"

#include <chrono>
#include <iostream>
#include <random>
#include <algorithm>

void Examples::exampleSurfaceComputationCircular()
{
    // Чтение сетки
    Grid::Centroids centroids = Grid::readCentroids("gridCentroidsNew.txt");

    // Чтение и печать орбит
    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    Orbits::printCircularOrbits(orbits);

    Settings::Sets settings = Settings::readSettings("settings.ini");
    Settings::printSettings(settings);

    // Вычисление области покрытия.
    //Surface::Surf surface = Surface::compute(centroids, orbits, settings);

    // Вычисление максимального времени ожидания
    clock_t start = clock();
    double maxTime = Surface::computeTime(centroids, orbits, settings);
    clock_t end = clock();

    std::cout << "\nTiming: " << (end - start) / (CLOCKS_PER_SEC / 1000) << "ms.\n";
    //out << "Result: " << Surface::sumArea(sphereGrid, surface) << ")\n";
    std::cout << "Max time: " << maxTime << "\n";
    //Surface::writeToTextFile(surface, "computedSurface.txt");
}

//void Examples::exampleSurfaceComputationElliptical()
//{
//    QTextStream out(stdout);
//    QTime timer;

//    // Чтение сетки
//    bool ok = false;
//    Grid::SphereGrid sphereGrid = Grid::readSphereGrid({"gridCentroids.txt", "gridAreas.txt"}, &ok);

//    // Чтение и печать орбит
//    std::vector<Orbits::EllipticalOrbit> orbits = Orbits::readEllipticalOrbits("ellipticalOrbits.txt");
//    Orbits::printEllipticalOrbits(orbits);

//    // Вычисление области покрытия и сохранение результата в файл для визуализации в Mathematica.
//    SatelliteSurface::Settings settings = SatelliteSurface::readSettings("settings.ini");
//    SatelliteSurface::printAlgorithmSettings(settings);

//    timer.start();
//    SatelliteSurface::Surface surface = SatelliteSurface::compute(sphereGrid.centroids, orbits, settings);
//    //double maxTime = SatelliteSurface::computeTime(sphereGrid.centroids, orbits, settings);

//    out << "\nTime: " << timer.elapsed() << "ms.\n";
//    //out << "MaxT = " << maxTime << "\n";
//    SatelliteSurface::writeToTextFile(surface, "computedSurface.txt");
//}

double penalty(const double a, const double x, const double b) {
    return std::max(0.0, x - b)*std::max(0.0, x - b) + 
           std::max(0.0, a - x)*std::max(0.0, a - x);
}

double ftime(const Opt::Point& orbitsVector)
{
    static const Grid::Centroids centroids = Grid::readCentroids("gridCentroidsNew.txt");
    static const Grid::Areas areas = Grid::readAreas("gridAreasNew.txt");

    static Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    static const Settings::Sets sets = Settings::readSettings("settings.ini");

    // =======================================

    // orbits[0].ascendingNode = orbitsVector[0];
    orbits[1].ascendingNode = orbitsVector[0];
    orbits[2].ascendingNode = orbitsVector[1];
    orbits[3].ascendingNode = orbitsVector[2];
    orbits[4].ascendingNode = orbitsVector[3];
    orbits[5].ascendingNode = orbitsVector[4];

    orbits[0].initialPhase = orbitsVector[5];
    orbits[1].initialPhase = orbitsVector[6];
    orbits[2].initialPhase = orbitsVector[7];
    orbits[3].initialPhase = orbitsVector[8];
    orbits[4].initialPhase = orbitsVector[9];
    orbits[5].initialPhase = orbitsVector[10];

    double p = 0.0;
    // p += penalty(0.0, orbitsVector[0], M_PI);
    p += penalty(0.0, orbitsVector[0], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[1], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[2], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[3], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[4], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[5], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[6], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[7], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[8], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[9], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[10], 2.0*M_PI);

    p += penalty(orbitsVector[0], orbitsVector[1], orbitsVector[2]);
    p += penalty(orbitsVector[1], orbitsVector[2], orbitsVector[3]);
    p += penalty(orbitsVector[2], orbitsVector[3], orbitsVector[4]);
    // p += penalty(orbitsVector[0], orbitsVector[1], orbitsVector[2]);
    // p += penalty(orbitsVector[1], orbitsVector[2], orbitsVector[3]);
    // p += penalty(orbitsVector[2], orbitsVector[3], orbitsVector[4]);
    // p += penalty(orbitsVector[0], orbitsVector[1], 2.0*M_PI);

    // =======================================

    const double areaPenaltyCoef = 10000.0;
    if (p > 0.001)
        return sets.timeDuration + areaPenaltyCoef + p*1.0E7;

    const double reducedTimeStep = 240.0;
    const auto surf = Surface::compute(centroids, orbits, {sets.coneAngle, sets.timeDuration, reducedTimeStep});
    
    const double fullArea = 12.5513;
    const double area = Surface::sumArea(areas, surf) / fullArea;

    double TStar = 0.0;
    if (area < 0.99) { TStar = sets.timeDuration + (1.0 - area)*areaPenaltyCoef; }
    else             { TStar = Surface::computeTime(centroids, orbits, sets); }

    return TStar;
}

void Examples::SwarmOptimisation()
{
    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    std::cout << "Loaded orbits:\n";
    Orbits::printCircularOrbits(orbits);

    // задаём область поиска в методе роя частиц
    SwarmMethod::Region region;

    // восходящие узлы плоскостей. Восходящий узел одной из них положим равным нулю
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});

    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});

    SwarmMethod::Parameters p;
    p.swarmSize = 20000;
    p.omega = -0.32;
    p.phi = 2.0;
    p.maxIterations = 20;

    const bool showIterations = true;
    const bool writeSwarmToFile = false; 
    const int runs = 1;
    for (int run = 0; run < runs; run++) {
        auto sol = SwarmMethod::optimize(ftime, region, p, showIterations, writeSwarmToFile);
        std::cout << "\n";
        std::cout.flush();
    }
    std::cout << "End of optimisation\n";
    // std::cout << "Optimization result: " << sol.second << "\n";
}

void Examples::exampleGridGeneration(const int iterations)
{
    auto grid = Grid::generate(iterations);
    auto centroids = Grid::centroids(grid);
    auto areas = Grid::areas(grid);

    Grid::writeNodes(grid.nodes, "gridNodesNew.txt");
    Grid::writeCells(grid.cells, "gridCellsNew.txt");
    Grid::writeAreas(areas, "gridAreasNew.txt");
    Grid::writeCentroids(centroids, "gridCentroidsNew.txt");
}
