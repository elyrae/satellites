#include "examples.h"
#include "satellitesurface.h"
#include "mathstuff.h"
#include "swarm.h"
#include "grid.h"

#include <iostream>
#include <algorithm>

void Examples::surface_computation()
{
    // Чтение сетки
    Grid::Centroids centroids = Grid::read_centroids("../data/grids/centroids_5000.txt");

    // Чтение и печать орбит
    Orbits::Constellation orbits = Orbits::read_circular_orbits("circular_orbits.txt");
    Orbits::print_circular_orbits(orbits);

    Settings::Sets settings = Settings::read_settings("settings.ini");
    Settings::print_settings(settings);

    // Вычисление области покрытия.
    //Surface::Surf surface = Surface::compute(centroids, orbits, settings);

    // Вычисление максимального времени ожидания
    clock_t start = clock();
    double max_time = Surface::compute_time(centroids, orbits, settings);
    clock_t end = clock();

    std::cout << "\nTiming: " << (end - start) / (CLOCKS_PER_SEC / 1000) << "ms.\n";
    //out << "Result: " << Surface::sumArea(sphereGrid, surface) << ")\n";
    std::cout << "Max time: " << max_time << "\n";
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

// double penalty(const double a, const double x, const double b) {
//     return std::max(0.0, x - b)*std::max(0.0, x - b) + 
//            std::max(0.0, a - x)*std::max(0.0, a - x);
// }

double ftime(const SwarmMethod::Point& orbitsVector)
{
    static const Grid::Centroids centroids = Grid::read_centroids("../data/grids/centroids_5000.txt");
    static const Grid::Areas areas = Grid::read_areas("../data/grids/areas_5000.txt");

    static Orbits::Constellation orbits = Orbits::read_circular_orbits("circularOrbits.txt");
    static const Settings::Sets sets = Settings::read_settings("settings.ini");

    // =======================================

    // orbits[0].ascendingNode = orbitsVector[0];
    // orbits[1].ascendingNode = orbitsVector[0];
    // orbits[2].ascendingNode = orbitsVector[0];
    // orbits[3].ascendingNode = orbitsVector[0];
    // orbits[4].ascendingNode = orbitsVector[0];
    // orbits[5].ascendingNode = orbitsVector[0];

    orbits[5].ascending_node  = orbits_vector[0];
    orbits[6].ascending_node  = orbits_vector[0];
    orbits[7].ascending_node  = orbits_vector[0];
    orbits[8].ascending_node  = orbits_vector[0];
    orbits[9].ascending_node  = orbits_vector[0];
    // orbits[10].ascendingNode = orbitsVector[0];

    orbits[10].ascending_node = orbits_vector[1];
    orbits[11].ascending_node = orbits_vector[1];
    orbits[12].ascending_node = orbits_vector[1];
    orbits[13].ascending_node = orbits_vector[1];
    orbits[14].ascending_node = orbits_vector[1];
    // orbits[15].ascendingNode = orbitsVector[1];

    orbits[15].ascending_node = orbits_vector[2];
    orbits[16].ascending_node = orbits_vector[2];
    orbits[17].ascending_node = orbits_vector[2];
    orbits[18].ascending_node = orbits_vector[2];
    orbits[19].ascending_node = orbits_vector[2];
    // orbits[20].ascendingNode = orbitsVector[2];

    // ============================================

    orbits[0].initialPhase  = orbitsVector[3] + 0.0*orbitsVector[7];
    orbits[1].initialPhase  = orbitsVector[3] + 1.0*orbitsVector[7];
    orbits[2].initialPhase  = orbitsVector[3] + 2.0*orbitsVector[7];
    orbits[3].initialPhase  = orbitsVector[3] + 3.0*orbitsVector[7];
    orbits[4].initialPhase  = orbitsVector[3] + 4.0*orbitsVector[7];
    // orbits[5].initialPhase  = orbitsVector[3] + 5.0*orbitsVector[7];

    orbits[5].initialPhase  = orbitsVector[4] + 0.0*orbitsVector[8];
    orbits[6].initialPhase  = orbitsVector[4] + 1.0*orbitsVector[8];
    orbits[7].initialPhase  = orbitsVector[4] + 2.0*orbitsVector[8];
    orbits[8].initialPhase  = orbitsVector[4] + 3.0*orbitsVector[8];
    orbits[9].initialPhase = orbitsVector[4] + 4.0*orbitsVector[8];
    // orbits[10].initialPhase = orbitsVector[4] + 5.0*orbitsVector[8];

    orbits[10].initialPhase = orbitsVector[5] + 0.0*orbitsVector[9];
    orbits[11].initialPhase = orbitsVector[5] + 1.0*orbitsVector[9];
    orbits[12].initialPhase = orbitsVector[5] + 2.0*orbitsVector[9];
    orbits[13].initialPhase = orbitsVector[5] + 3.0*orbitsVector[9];
    orbits[14].initialPhase = orbitsVector[5] + 4.0*orbitsVector[9];
    // orbits[15].initialPhase = orbitsVector[5] + 5.0*orbitsVector[9];

    orbits[15].initialPhase = orbitsVector[6] + 0.0*orbitsVector[10];
    orbits[16].initialPhase = orbitsVector[6] + 1.0*orbitsVector[10];
    orbits[17].initialPhase = orbitsVector[6] + 2.0*orbitsVector[10];
    orbits[18].initialPhase = orbitsVector[6] + 3.0*orbitsVector[10];
    orbits[19].initialPhase = orbitsVector[6] + 4.0*orbitsVector[10];
    // orbits[20].initialPhase = orbitsVector[6] + 5.0*orbitsVector[10];

    double p = 0.0;
    p += penalty(0.0, orbitsVector[0], M_PI);
    p += penalty(0.0, orbitsVector[1], M_PI);
    p += penalty(0.0, orbitsVector[2], M_PI);

    p += penalty(0.0, orbitsVector[3], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[4], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[5], 2.0*M_PI);
    p += penalty(0.0, orbitsVector[6], 2.0*M_PI);

    p += penalty(MathStuff::degrees_to_rad(50.0), orbitsVector[7], MathStuff::degrees_to_rad(80.0));
    p += penalty(MathStuff::degrees_to_rad(50.0), orbitsVector[8], MathStuff::degrees_to_rad(80.0));
    p += penalty(MathStuff::degrees_to_rad(50.0), orbitsVector[9], MathStuff::degrees_to_rad(80.0));
    p += penalty(MathStuff::degrees_to_rad(50.0), orbitsVector[10],MathStuff::degrees_to_rad(80.0));

    p += penalty(orbitsVector[0], orbitsVector[1], orbitsVector[2]);

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

// void Examples::swarm_optimisation()
// {
//     Orbits::Constellation orbits = Orbits::read_circular_orbits("circularOrbits.txt");
//     std::cout << "Loaded orbits:\n";
//     Orbits::print_circular_orbits(orbits);

//     // задаём область поиска в методе роя частиц
//     SwarmMethod::Region region;

//     // восходящие узлы плоскостей. Восходящий узел одной из них положим равным нулю
//     region.push_back({0.0, 1.0*M_PI});
//     region.push_back({0.0, 1.0*M_PI});
//     region.push_back({0.0, 1.0*M_PI});

//     // начальные фазы
//     region.push_back({0.0, 2.0*M_PI});
//     region.push_back({0.0, 2.0*M_PI});
//     region.push_back({0.0, 2.0*M_PI});
//     region.push_back({0.0, 2.0*M_PI});

//     // шаг начальной фазы для каждой плоскости
//     region.push_back({MathStuff::degrees_to_rad(50.0), MathStuff::degrees_to_rad(80.0)});
//     region.push_back({MathStuff::degrees_to_rad(50.0), MathStuff::degrees_to_rad(80.0)});
//     region.push_back({MathStuff::degrees_to_rad(50.0), MathStuff::degrees_to_rad(80.0)});
//     region.push_back({MathStuff::degrees_to_rad(50.0), MathStuff::degrees_to_rad(80.0)});

//     SwarmMethod::Parameters p;
//     p.swarm_size = 500;
//     p.omega = -0.32;
//     p.phi = 2.0;
//     p.max_iterations = 20;

//     const bool show_iterations = true;
//     const bool write_swarm = false; 
//     const int runs = 1;
//     for (int run = 0; run < runs; run++) {
//         auto sol = SwarmMethod::optimize(ftime, region, p, show_iterations, write_swarm);
//         std::cout << "\n";
//         std::cout.flush();
//     }
//     std::cout << "End of optimisation\n";
// }

void Examples::grid_generation(const int iterations)
{
    auto grid = Grid::generate(iterations);
    auto centroids = Grid::centroids(grid);
    auto areas = Grid::areas(grid);

    Grid::write_nodes(grid.nodes,    "gridNodesNew.txt");
    Grid::write_cells(grid.cells,    "gridCellsNew.txt");
    Grid::write_areas(areas,         "gridAreasNew.txt");
    Grid::write_centroids(centroids, "gridCentroidsNew.txt");
}
