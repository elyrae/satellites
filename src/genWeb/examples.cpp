#include "examples.h"
#include "satellitesurface.h"
#include "mathstuff.h"
#include "swarm.h"
#include "grid.h"

#include <QTime>
#include <iostream>
#include <random>

void Examples::exampleSurfaceComputationCircular()
{
    QTime timer;

    // Чтение сетки
    //bool ok = false;
    //Grid::SphereGrid sphereGrid = Grid::readGrid({"gridCentr-small.txt", "gridAreas-small.txt"}, &ok);
    Grid::Centroids centroids = Grid::readCentroids("gridCentr-small.txt");

    // Чтение и печать орбит
    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    Orbits::printCircularOrbits(orbits);

    Settings::Sets settings = Settings::readSettings("settings.ini");
    Settings::printSettings(settings);

    timer.start();

    // Вычисление области покрытия и сохранение результата в файл для визуализации в Mathematica.
    //Surface::Surf surface = Surface::compute(sphereGrid.centroids, orbits, settings);
    //Surface::Surf surface = Surface::compute(centroids, orbits, settings);
    // Время ожидания:
    //double maxTime = Surface::computeTime(sphereGrid.centroids, orbits, settings);
    double maxTime = Surface::computeTime(centroids, orbits, settings);
    auto elapsed = timer.elapsed();

    std::cout << "\nTime: " << elapsed << "ms.\n";
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

double max(const std::pair<double, double> p)
{
    return (p.first > p.second) ? p.first : p.second;
}

inline double sqr(const double x)
{
    return x*x;
}

double ftime(const Opt::Point& orbitsVector)
{
    const double fullArea = 12.5513;
    const double reducedTimeStep = 240.0;

    static const Grid::Centroids centroids = Grid::readCentroids("gridCentr-small.txt");
    static const Grid::Areas areas = Grid::readAreas("gridAreas-small.txt");

    static Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    static const Settings::Sets p = Settings::readSettings("settings.ini");

    // orbits[0].ascendingNode = orbitsVector[0];
    // orbits[1].ascendingNode = orbitsVector[0];
    // orbits[2].ascendingNode = orbitsVector[0];
    // orbits[3].ascendingNode = orbitsVector[0];
    // orbits[4].ascendingNode = orbitsVector[0];
    // orbits[5].ascendingNode = orbitsVector[0];

    orbits[5].ascendingNode = orbitsVector[0];
    orbits[6].ascendingNode = orbitsVector[0];
    orbits[7].ascendingNode = orbitsVector[0];
    orbits[8].ascendingNode = orbitsVector[0];
    orbits[9].ascendingNode = orbitsVector[0];
    // orbits[11].ascendingNode = orbitsVector[0];

    orbits[10].ascendingNode = orbitsVector[1];
    orbits[11].ascendingNode = orbitsVector[1];
    orbits[12].ascendingNode = orbitsVector[1];
    orbits[13].ascendingNode = orbitsVector[1];
    orbits[14].ascendingNode = orbitsVector[1];
    // orbits[17].ascendingNode = orbitsVector[1];

    orbits[15].ascendingNode = orbitsVector[2];
    orbits[16].ascendingNode = orbitsVector[2];
    orbits[17].ascendingNode = orbitsVector[2];
    orbits[18].ascendingNode = orbitsVector[2];
    orbits[19].ascendingNode = orbitsVector[2];
    // orbits[23].ascendingNode = orbitsVector[2];
    // ============================================

    orbits[0].initialPhase = orbitsVector[3] + 0.0*orbitsVector[7];
    orbits[1].initialPhase = orbitsVector[3] + 1.0*orbitsVector[7];
    orbits[2].initialPhase = orbitsVector[3] + 2.0*orbitsVector[7];
    orbits[3].initialPhase = orbitsVector[3] + 3.0*orbitsVector[7];
    orbits[4].initialPhase = orbitsVector[3] + 4.0*orbitsVector[7];
    // orbits[5].initialPhase = orbitsVector[3] + 5.0*orbitsVector[7];

    orbits[5].initialPhase = orbitsVector[4] + 0.0*orbitsVector[8];
    orbits[6].initialPhase = orbitsVector[4] + 1.0*orbitsVector[8];
    orbits[7].initialPhase = orbitsVector[4] + 2.0*orbitsVector[8];
    orbits[8].initialPhase = orbitsVector[4] + 3.0*orbitsVector[8];
    orbits[9].initialPhase = orbitsVector[4] + 4.0*orbitsVector[8];
    // orbits[11].initialPhase = orbitsVector[4] + 5.0*orbitsVector[8];

    orbits[10].initialPhase = orbitsVector[5] + 0.0*orbitsVector[9];
    orbits[11].initialPhase = orbitsVector[5] + 1.0*orbitsVector[9];
    orbits[12].initialPhase = orbitsVector[5] + 2.0*orbitsVector[9];
    orbits[13].initialPhase = orbitsVector[5] + 3.0*orbitsVector[9];
    orbits[14].initialPhase = orbitsVector[5] + 4.0*orbitsVector[9];
    // orbits[17].initialPhase = orbitsVector[5] + 5.0*orbitsVector[9];

    orbits[15].initialPhase = orbitsVector[6] + 0.0*orbitsVector[10];
    orbits[16].initialPhase = orbitsVector[6] + 1.0*orbitsVector[10];
    orbits[17].initialPhase = orbitsVector[6] + 2.0*orbitsVector[10];
    orbits[18].initialPhase = orbitsVector[6] + 3.0*orbitsVector[10];
    orbits[19].initialPhase = orbitsVector[6] + 4.0*orbitsVector[10];
    // orbits[23].initialPhase = orbitsVector[6] + 5.0*orbitsVector[10];

    double penalty = 0.0;
    const int params_out_region = 11;
    for (int i = 0; i < params_out_region; i++) {
        penalty += sqr( max({0.0, -orbitsVector[i]}) );
        penalty += sqr( max({0.0,  orbitsVector[i] - 2.0*M_PI}) );
    }
    const int opt_psis = 3;
    for (int i = 0; i < opt_psis - 1; i++)
        penalty += sqr( max({0.0, orbitsVector[i] - orbitsVector[i+1]}) );

    if (penalty > 0.01)
        return penalty*1.0E8;

    const Surface::Surf surf = Surface::compute(centroids, orbits,
        {p.coneAngle, p.timeDuration, reducedTimeStep});
    const double area = Surface::sumArea(areas, surf) / fullArea;

//    const double area = Surface::computeArea(sphereGrid, orbits,
//        {p.coneAngle, p.timeDuration, reducedTimeStep}) / fullArea;

    double TStar = 0.0;
    if (area < 0.99) { TStar = p.timeDuration + (1.0 - area)*10000.0; }
    else             { TStar = Surface::computeTime(centroids, orbits, p); }

    return TStar; //+ penalty*1.0E8;
}

//double ftimeF(const Opt::Point& orbitsVector)
//{
//    static bool ok = false;
//    static const Grid::SphereGrid sphereGrid = Grid::readGrid({"bigGridCentroids.txt", "bigGridAreas.txt"}, &ok);
//    static std::vector<Orbits::CircularOrbit> orbits = Orbits::readCircularOrbits("circularOrbits.txt");
//    static const Settings::Sets p = Settings::readSettings("settings.ini");

//    const int opt_params = 4;

//    orbits[2].ascendingNode = orbitsVector[0];
//    orbits[3].ascendingNode = orbitsVector[0];

//    orbits[4].ascendingNode = orbitsVector[1];
//    orbits[5].ascendingNode = orbitsVector[1];

//    orbits[6].ascendingNode = orbitsVector[2];
//    orbits[7].ascendingNode = orbitsVector[2];

//    orbits[8].ascendingNode = orbitsVector[3];
//    orbits[9].ascendingNode = orbitsVector[3];

//    return Surface::computeTime(sphereGrid.centroids, orbits, p);
//}

//void out_histogram(const Orbits::Constellation& orbits)
//{
//    // Чтение сетки
//    bool ok = false;
//    Grid::SphereGrid sphereGrid = Grid::readGrid({"bigGridCentroids.txt", "bigGridAreas.txt"}, &ok);
//    if (!ok)
//        return;

//    // Чтение и печать орбит
////    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
//    Orbits::printCircularOrbits(orbits);

//    // Вычисление области покрытия и сохранение результата в файл для визуализации в Mathematica.
//    Settings::Sets settings = Settings::readSettings("settings.ini");
//    settings.deltaT = 15.0;
//    settings.timeDuration = 7200.0;
//    //SatelliteSurface::Surface surface = SatelliteSurface::compute(sphereGrid.centroids, orbits, parameters);

//    QFile histoFile("hist_time.txt");
//    histoFile.open(QIODevice::WriteOnly | QIODevice::Text);
//    QTextStream outHist(&histoFile);

//    auto time_data = Surface::computeTimeFull(sphereGrid.centroids, orbits, settings);

//    outHist << "{";
//    for (size_t i = 0; i < time_data.size(); i++) {
//        outHist << "{" << sphereGrid.centroids[i][0] << "," << sphereGrid.centroids[i][1] << ","
//                << sphereGrid.centroids[i][2] << ","
//                << time_data[i] << /*( (i < time_data.size()-1) ? "," : "" ) <<*/ "}";
//        outHist << ((i < time_data.size() - 1) ? "," : "");
//    }
//        //outHist << time_data[i] << ( (i < time_data.size()-1) ? "," : "" );
//    outHist << "}";
//    histoFile.close();
//}

void Examples::SwarmOptimisation()
{
    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    Orbits::printCircularOrbits(orbits);

    ParticleSwarmMethod::Region region;
//    const int optPars = 11;
//    for (int i = 0; i < optPars; i++) {
//        region.push_back({0.0, 2.0*M_PI});
//    }
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});

    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI});

    region.push_back({0.0, MathStuff::degreesToRad(100.0)});
    region.push_back({0.0, MathStuff::degreesToRad(100.0)});
    region.push_back({0.0, MathStuff::degreesToRad(100.0)});
    region.push_back({0.0, MathStuff::degreesToRad(100.0)});

    ParticleSwarmMethod::Parameters p;
    p.S = 1000;
    p.omega = -0.32;
    p.phi = 2.0;
    p.maxIterations = 20;
    auto sol = ParticleSwarmMethod::optimize(ftime, region, Opt::SearchType::SearchMinimum, p, true, false);

    std::cout << "Answer: " << sol.second << "\n";

//    QTextStream out(stdout);
//    out << "Answer: " << ftimeF(sol.first) << "\n";

//    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
//    orbits[2].ascendingNode = orbitsVector[0];
//    orbits[3].ascendingNode = orbitsVector[0];

//    orbits[4].ascendingNode = orbitsVector[1];
//    orbits[5].ascendingNode = orbitsVector[1];

//    orbits[6].ascendingNode = orbitsVector[2];
//    orbits[7].ascendingNode = orbitsVector[2];

//    orbits[8].ascendingNode = orbitsVector[3];
//    orbits[9].ascendingNode = orbitsVector[3];
//    out_histogram(orbits);
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
