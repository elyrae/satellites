#include "cudatest.h"
#include "orbits.h"

#include <stdlib.h>
#include <cmath>
#include <vector>

const int SURFACE_SIZE       = 5120;
const int CONFIGURATION_SIZE = 20;
const int CONFIGURATIONS     = 256;
const int TIMESTEPS          = 240;

double rnd_d(const double a, const double b)
{
    return a + (b - a)*((double) rand()/RAND_MAX);
}

void random_fill(Orbits::ConstellationUnion &configs)
{
    for (size_t iconf = 0; iconf < CONFIGURATIONS; ++iconf) {
        configs[iconf].resize(CONFIGURATION_SIZE);
        for (size_t iorb = 0; iorb < CONFIGURATION_SIZE; ++iorb) {
            configs[iconf][iorb].ascending_node = rnd_d(0.0, 2.0*M_PI);
            configs[iconf][iorb].inclination    = rnd_d(0.0, 2.0*M_PI);
            configs[iconf][iorb].initial_phase  = rnd_d(0.0, 2.0*M_PI);
            configs[iconf][iorb].height         = 1500.0 * 1000.0;
        }
    }
}

// void measure_time()
// {
//     Grid::Centroids_f centroids = Grid::readCentroids_f("../data/grids/centroids_5000.txt");
//     Settings::Sets sets = Settings::read_settings("settings.ini");
//     Settings::printSettings(sets);

//     const size_t sat_positions_size = CONFIGURATION_SIZE*CONFIGURATIONS*TIMESTEPS;
//     const size_t surf_size          =       SURFACE_SIZE*CONFIGURATIONS*TIMESTEPS;
//     std::vector<float> x(sat_positions_size), y(sat_positions_size), z(sat_positions_size), h(sat_positions_size);
//     std::vector<float> surf(surf_size);

//     std::vector<Orbits::Constellation> configs(CONFIGURATIONS);
//     random_fill(configs);
//     fill_orbits(configs, sets, x, y, z, h);

//     printf("\n");
//     printf("%f MB for satellite positions\n", (3.0*sizeof(float)*sat_positions_size) / (1024.0*1024.0));
//     printf("%f MB for surface flags\n",                    (sizeof(float)*surf_size) / (1024.0*1024.0));
//     printf("%f MB for centroids\n",           (3.0*sizeof(float)*SURFACE_SIZE) / (1024.0*1024.0));
//     printf("\n");

//     CUDA_Surface::Points3D gpu_pos,       cpu_pos(sat_positions_size, x.data(), y.data(), z.data()),
//                            gpu_centroids, cpu_centroids(centroids.X.size(), centroids.X.data(), centroids.Y.data(), centroids.Z.data()); 
//     CUDA_Surface::Points1D gpu_surf,      cpu_surf(surf_size, surf.data());
    
//     gpu_pos.allocate(sat_positions_size);
//     gpu_surf.allocate(surf_size);
//     gpu_centroids.allocate(SURFACE_SIZE);
//     printf("Success\n");

//     gpu_pos.load_from(cpu_pos);
//     gpu_centroids.load_from(cpu_centroids);
//     printf("Success load to GPU\n");

//     const size_t max_time_size = SURFACE_SIZE*CONFIGURATIONS;
//     CUDA_Surface::Points1D gpu_max_time;
//     gpu_max_time.allocate(max_time_size);

//     float gpu_time = CUDA_Surface::compute_surface(gpu_surf, gpu_max_time, gpu_centroids, gpu_pos);
//     printf("Success. gpu time = %f\n", gpu_time);

//     gpu_pos.free();
//     gpu_surf.free();
//     gpu_centroids.free();
//     gpu_max_time.free();
//     printf("Memory was successfully deallocated\n");    
// }