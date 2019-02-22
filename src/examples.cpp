#include "examples.h"
#include "surface.h"
#include "stuff.h"
#include "swarm.h"
#include "grid.h"

#include <iostream>
#include <fstream>
#include <algorithm>

void Examples::time_computation()
{
    // Чтение сетки
    Grid::Centroids centroids = Grid::read_centroids("../data/grids/centroids_5000.txt");

    // Чтение и печать орбит
    Orbits::Constellation orbits = Orbits::read_circular_orbits("circular_orbits.txt");
    Orbits::print_circular_orbits(orbits);

    Settings::Sets settings = Settings::read_settings("settings.ini");
    Settings::print_settings(settings);

    // Вычисление максимального времени ожидания
    clock_t start = clock();
    double max_time = Surface::compute_time(centroids, orbits, settings);
    clock_t end = clock();

    std::cout << "\nTiming: " << (end - start) / (CLOCKS_PER_SEC / 1000) << "ms.\n";
    std::cout << "Max time: " << max_time << "\n";
}

void Examples::time_histogram()
{
    Grid::Centroids centroids = Grid::read_centroids("../data/grids/centroids_80000.txt");

    Orbits::Constellation orbits = Orbits::read_circular_orbits("circular_orbits.txt");
    Orbits::print_circular_orbits(orbits);

    Settings::Sets settings = Settings::read_settings("settings.ini");
    Settings::print_settings(settings);

    settings.elevation_angle = 0.0;
    Surface::Times elev0_1500  = Surface::compute_time_full(centroids, orbits, settings);
    settings.elevation_angle = 10.0;
    Surface::Times elev10_1500 = Surface::compute_time_full(centroids, orbits, settings);

    for (Orbits::CircularOrbit &orbit: orbits)
        orbit.height = 2000.0 * 1000.0;

    settings.elevation_angle = 0.0;
    Surface::Times elev0_2000  = Surface::compute_time_full(centroids, orbits, settings);
    settings.elevation_angle = 10.0;
    Surface::Times elev10_2000 = Surface::compute_time_full(centroids, orbits, settings);
    
    std::ofstream out("times.txt");
    for (size_t i = 0; i < elev0_1500.size(); i++)
        out << elev0_1500[i] << " " << elev10_1500[i] << " " << elev0_2000[i] << " " << elev10_2000[i] << "\n";    
}

void Examples::surface_computation()
{
    Grid::Centroids centroids = Grid::read_centroids("../data/grids/centroids_5000.txt");

    Orbits::Constellation orbits = Orbits::read_circular_orbits("circular_orbits.txt");
    Orbits::print_circular_orbits(orbits);

    Settings::Sets settings = Settings::read_settings("settings.ini");
    Settings::print_settings(settings);

    // Вычисление области покрытия.
    clock_t start = clock();
    Surface::Surf surface = Surface::compute(centroids, orbits, settings);
    clock_t end = clock();

    std::cout << "\nTiming: " << (end - start) / (CLOCKS_PER_SEC / 1000) << "ms.\n";
    // out << "Result: " << Surface::sumArea(sphereGrid, surface) << ")\n";
    //Surface::writeToTextFile(surface, "computedSurface.txt");
}

void Examples::compare_different_grids()
{
    Grid::TriangularGrid grid;
    Grid::Centroids centroids;
    Surface::Times times;

    Orbits::Constellation orbits = Orbits::read_circular_orbits("circular_orbits.txt");
    Orbits::print_circular_orbits(orbits);

    Settings::Sets settings = Settings::read_settings("settings.ini");
    Settings::print_settings(settings);

    std::ofstream out("different_grids.txt");
    for (int it = 0; it < 7; it++) {
        grid = Grid::generate(it);
        centroids = Grid::centroids(grid);
        std::cout << it << ": " << centroids.X.size() << "\n";

        times = Surface::compute_time_full(centroids, orbits, settings);
        for (size_t i = 0; i < times.size(); i++)
            out << times[i] << " ";
        out << "\n"; 
    }    
}

void Examples::destroy_satellites()
{
    Grid::TriangularGrid grid;
    Grid::Centroids centroids;

    const size_t iters = 3;
    centroids = Grid::centroids(Grid::generate(iters));
    std::cout << centroids.X.size() << " centroids\n";

    Orbits::Constellation orbits = Orbits::read_circular_orbits("circular_orbits.txt");
    Orbits::print_circular_orbits(orbits);

    Settings::Sets settings = Settings::read_settings("settings.ini");
    Settings::print_settings(settings);

    auto full_constellation = orbits;
    std::cout << "Full constellation: " << Surface::compute_time(centroids, orbits, settings) << "\n";
    for (size_t i = 0; i < orbits.size(); i++) {
        orbits.erase(orbits.begin() + i);
        Orbits::print_circular_orbits(orbits);

        std::cout << i << ": " << Surface::compute_time(centroids, orbits, settings) << "\n\n";
        orbits = full_constellation;
    }
}

double penalty(const double a, const double x, const double b) {
    return std::max(0.0, x - b)*std::max(0.0, x - b) + 
           std::max(0.0, a - x)*std::max(0.0, a - x);
}

double time_function(Surface::Timegrid &tg, const SwarmMethod::Point &vect)
{
    Orbits::Constellation orbits(tg.default_constellation());

    orbits[0].ascending_node = 0.0;
    orbits[1].ascending_node = 0.0;
    orbits[2].ascending_node = 0.0;

    orbits[3].ascending_node = vect[0];
    orbits[4].ascending_node = vect[0];
    orbits[5].ascending_node = vect[0];

    orbits[6].ascending_node = vect[1];
    orbits[7].ascending_node = vect[1];
    orbits[8].ascending_node = vect[1];

    orbits[9].ascending_node  = vect[2];
    orbits[10].ascending_node = vect[2];
    orbits[11].ascending_node = vect[2];    

    orbits[0].initial_phase = vect[3] + 0.0*vect[4];
    orbits[1].initial_phase = vect[3] + 1.0*vect[4];
    orbits[2].initial_phase = vect[3] + 2.0*vect[4];

    orbits[3].initial_phase = vect[5] + 0.0*vect[6];
    orbits[4].initial_phase = vect[5] + 1.0*vect[6];
    orbits[5].initial_phase = vect[5] + 2.0*vect[6];

    orbits[6].initial_phase = vect[7] + 0.0*vect[8];
    orbits[7].initial_phase = vect[7] + 1.0*vect[8];
    orbits[8].initial_phase = vect[7] + 2.0*vect[8];

    orbits[9].initial_phase  = vect[9] + 0.0*vect[10];
    orbits[10].initial_phase = vect[9] + 1.0*vect[10];
    orbits[11].initial_phase = vect[9] + 2.0*vect[10];  

    double p = 0.0;
    p += penalty(0.0,     vect[0], vect[1]);
    p += penalty(vect[0], vect[1], vect[2]);
    p += penalty(vect[1], vect[2], M_PI);

    p += penalty(0.0, vect[3], 2.0*M_PI);
    p += penalty(0.0, vect[4], 2.0*M_PI/3.0);

    p += penalty(0.0, vect[5], 2.0*M_PI);
    p += penalty(0.0, vect[6], 2.0*M_PI/3.0);
    
    p += penalty(0.0, vect[7], 2.0*M_PI);
    p += penalty(0.0, vect[8], 2.0*M_PI/3.0);
    
    p += penalty(0.0, vect[9], 2.0*M_PI);
    p += penalty(0.0, vect[10], 2.0*M_PI/3.0);

    const double area_penalty = 1.0E4;
    if (p > 0.001)
        return tg.settings().time_duration + area_penalty + p*1.0E7;

    const double area = tg.compute_covered_area(orbits) / tg.full_area();
    if (area < 0.99)
         return tg.settings().time_duration + (1.0 - area)*area_penalty;
    else return tg.compute_max_time(orbits);
}

void Examples::swarm_optimisation()
{
    // задаём область поиска в методе роя частиц
    SwarmMethod::Region region;
    region.push_back({0.0, M_PI});
    region.push_back({0.0, M_PI});
    region.push_back({0.0, M_PI});
    
    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI/3.0});

    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI/3.0});

    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI/3.0});

    region.push_back({0.0, 2.0*M_PI});
    region.push_back({0.0, 2.0*M_PI/3.0});

    // ##############################
    SwarmMethod::Parameters swarm_params;
    swarm_params.swarm_size = 5000;
    swarm_params.omega = SwarmMethod::defaults.omega;
    swarm_params.phi = SwarmMethod::defaults.phi;

    swarm_params.max_iterations = 15;
    swarm_params.debug = true;

    // ##############################
    Orbits::Constellation orbits = Orbits::read_circular_orbits("circular_orbits.txt");
    Orbits::print_circular_orbits(orbits);

    Settings::Sets satellite_settings = Settings::read_settings("settings.ini");
    Settings::print_settings(satellite_settings);

    const int grid_iters = 3;
    Surface::Timegrid tg(grid_iters, orbits, satellite_settings);
    // ##############################

    const int runs = 1;
    SwarmMethod::Optimizer opt(tg, region, swarm_params);
    for (int run = 0; run < runs; run++) {
        auto sol = opt.optimize(time_function); // SwarmMethod::optimize(ftime, region, params, show_iterations);
        std::cout << "\n";
        std::cout.flush();
    }
    std::cout << "End of optimisation\n";
}

void Examples::grid_generation(const int iterations)
{
    auto grid = Grid::generate(iterations);
    auto centroids = Grid::centroids(grid);
    auto areas = Grid::areas(grid);

    Grid::write_nodes(grid.nodes,    "nodes_new.txt");
    Grid::write_cells(grid.cells,    "cells_new.txt");
    Grid::write_areas(areas,         "areas_new.txt");
    Grid::write_centroids(centroids, "centroids_new.txt");
}
