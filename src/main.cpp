#include <iostream>
#include <math.h>
#include <fstream>
#include <ostream>
#include <vector>

#include "grid.h"
#include "satellitesurface.h"
#include "orbits.h"
#include "earth.h"
#include "messages.h"
#include "mathstuff.h"
#include "examples.h"

#include "cuda_surface.cuh"

const int SURFACE_SIZE       = 5120;
const int CONFIGURATION_SIZE = 20;
const int CONFIGURATIONS     = 256;
const int TIMESTEPS          = 240;

inline double horizon(const double H, const double alpha)
{
    const double delta = MathStuff::degreesToRad(10.0); // требуемое возвышение спутника над горизонтом
    const double alpha_star = asin(cos(delta) / H);

    return (alpha < alpha_star) ? cos(asin(H*sin(alpha)) - alpha) : sin(delta + alpha_star);
    // return sin(delta + alpha_star);
}

// void fill_orbits(const std::vector<Orbits::Constellation> &configurations, const Settings::Sets &settings, 
//                  std::vector<float4> &p)
// {
//     // const double alpha = MathStuff::degreesToRad(settings.coneAngle) / 2.0;
//     const size_t configs = configurations.size();
//     for (size_t iconf = 0; iconf < configs; ++iconf) {
//         Orbits::Constellation orbits = configurations[iconf];

//         for (size_t iorb = 0; iorb < orbits.size(); ++iorb) {
//             double semi_major_axis = orbits[iorb].semiMajorAxis();
//             double mean_angular_velocity = sqrt(Earth::mu / (semi_major_axis*semi_major_axis*semi_major_axis));        
//             double cos_i = cos(orbits[iorb].inclination);
//             double sin_i = sin(orbits[iorb].inclination);

//             for (size_t timestep = 0; timestep < TIMESTEPS; ++timestep) {
//                 double t = timestep*settings.deltaT;
//                 double cos_node = cos(orbits[iorb].ascendingNode - Earth::angularVelocity*t);
//                 double sin_node = sin(orbits[iorb].ascendingNode - Earth::angularVelocity*t);
//                 double cos_anomaly_plus_phase = cos(mean_angular_velocity*t + orbits[iorb].initialPhase);
//                 double sin_anomaly_plus_phase = sin(mean_angular_velocity*t + orbits[iorb].initialPhase);

//                 p[(timestep*orbits.size() + iorb)*configs + iconf].x = (float) (      cos_node*cos_anomaly_plus_phase - cos_i*sin_anomaly_plus_phase*sin_node); 
//                 p[(timestep*orbits.size() + iorb)*configs + iconf].y = (float) (cos_i*cos_node*sin_anomaly_plus_phase +       cos_anomaly_plus_phase*sin_node); 
//                 p[(timestep*orbits.size() + iorb)*configs + iconf].z = (float) (sin_i*sin_anomaly_plus_phase                                                 );
//                 // h[(timestep*orbits.size() + iorb)*configs + iconf] = (float) (horizon(semi_major_axis / Earth::radius, alpha)                              );        
//             }
//         }
//     }
// }

void fill_orbits(const std::vector<Orbits::Constellation> &configurations, const Settings::Sets &settings, 
                 std::vector<float> &x, std::vector<float> &y, std::vector<float> &z, std::vector<float> &h)
{
    const double alpha = MathStuff::degreesToRad(settings.coneAngle) / 2.0;
    const size_t configs = configurations.size();
    for (size_t iconf = 0; iconf < configs; ++iconf) {
        Orbits::Constellation orbits = configurations[iconf];

        for (size_t iorb = 0; iorb < orbits.size(); ++iorb) {
            double semi_major_axis = orbits[iorb].semiMajorAxis();
            double mean_angular_velocity = sqrt(Earth::mu / (semi_major_axis*semi_major_axis*semi_major_axis));        
            double cos_i = cos(orbits[iorb].inclination);
            double sin_i = sin(orbits[iorb].inclination);

            for (size_t timestep = 0; timestep < TIMESTEPS; ++timestep) {
                double t = timestep*settings.deltaT;
                double cos_node = cos(orbits[iorb].ascendingNode - Earth::angularVelocity*t);
                double sin_node = sin(orbits[iorb].ascendingNode - Earth::angularVelocity*t);
                double cos_anomaly_plus_phase = cos(mean_angular_velocity*t + orbits[iorb].initialPhase);
                double sin_anomaly_plus_phase = sin(mean_angular_velocity*t + orbits[iorb].initialPhase);

                x[(timestep*orbits.size() + iorb)*configs + iconf] = (float) (      cos_node*cos_anomaly_plus_phase - cos_i*sin_anomaly_plus_phase*sin_node); 
                y[(timestep*orbits.size() + iorb)*configs + iconf] = (float) (cos_i*cos_node*sin_anomaly_plus_phase +       cos_anomaly_plus_phase*sin_node); 
                z[(timestep*orbits.size() + iorb)*configs + iconf] = (float) (sin_i*sin_anomaly_plus_phase                                                 );
                h[(timestep*orbits.size() + iorb)*configs + iconf] = (float) (horizon(semi_major_axis / Earth::radius, alpha)                              );        
            }
        }
    }
}

double rnd_d(const double a, const double b)
{
    return a + (b - a)*((double) rand()/RAND_MAX);
}

void random_fill(std::vector<Orbits::Constellation> &configs)
{
    for (size_t iconf = 0; iconf < CONFIGURATIONS; ++iconf) {
        configs[iconf].resize(CONFIGURATION_SIZE);
        for (size_t iorb = 0; iorb < CONFIGURATION_SIZE; ++iorb) {
            configs[iconf][iorb].ascendingNode = rnd_d(0.0, 2.0*M_PI);
            configs[iconf][iorb].inclination   = rnd_d(0.0, 2.0*M_PI);
            configs[iconf][iorb].initialPhase  = rnd_d(0.0, 2.0*M_PI);
            configs[iconf][iorb].height        = 1500.0 * 1000.0;
        }
    }
}

void print_first(const size_t n, const std::vector<float> &x, const std::vector<float> &y, const std::vector<float> &z)
{
    printf("x: ");
    for (size_t i = 0; i < ((n < x.size()) ? n : x.size()); ++i)
        printf("%5.2f ", x[i]);
    printf("\ny: ");
    for (size_t i = 0; i < ((n < x.size()) ? n : x.size()); ++i)
        printf("%5.2f ", y[i]);
    printf("\nz: ");
    for (size_t i = 0; i < ((n < x.size()) ? n : x.size()); ++i)
        printf("%5.2f ", z[i]);
    printf("\n");
}

void centroids_to_float3(const Grid::Centroids_f &c, std::vector<float4> &v)
{
    for (size_t i = 0; i < c.X.size(); ++i) {
        v[i].x = c.X[i];
        v[i].y = c.Y[i];
        v[i].z = c.Z[i];
    }
}

int main() {
    Grid::Centroids_f centroids = Grid::readCentroids_f("../data/grids/centroids_5000.txt");
    Settings::Sets sets = Settings::readSettings("settings.ini");
    Settings::printSettings(sets);

    const size_t sat_positions_size = CONFIGURATION_SIZE*CONFIGURATIONS*TIMESTEPS;
    const size_t surf_size          =       SURFACE_SIZE*CONFIGURATIONS*TIMESTEPS;
    std::vector<float> x(sat_positions_size), y(sat_positions_size), z(sat_positions_size), h(sat_positions_size);
    std::vector<float> surf(surf_size);

    std::vector<Orbits::Constellation> configs(CONFIGURATIONS);
    random_fill(configs);
    fill_orbits(configs, sets, x, y, z, h);

    printf("\n");
    printf("%f MB for satellite positions\n", (3.0*sizeof(float)*sat_positions_size) / (1024.0*1024.0));
    printf("%f MB for surface flags\n",                    (sizeof(float)*surf_size) / (1024.0*1024.0));
    printf("%f MB for centroids\n",           (3.0*sizeof(float)*SURFACE_SIZE) / (1024.0*1024.0));
    printf("\n");

    CUDA_Surface::Points3D gpu_pos,       cpu_pos(sat_positions_size, x.data(), y.data(), z.data()),
                           gpu_centroids, cpu_centroids(centroids.X.size(), centroids.X.data(), centroids.Y.data(), centroids.Z.data()); 
    CUDA_Surface::Points1D gpu_surf,      cpu_surf(surf_size, surf.data());
    
    gpu_pos.allocate(sat_positions_size);
    gpu_surf.allocate(surf_size);
    gpu_centroids.allocate(SURFACE_SIZE);
    printf("Success\n");

    gpu_pos.load_from(cpu_pos);
    gpu_centroids.load_from(cpu_centroids);
    printf("Success load to GPU\n");

    const size_t max_time_size = SURFACE_SIZE*CONFIGURATIONS;
    CUDA_Surface::Points1D gpu_max_time;
    gpu_max_time.allocate(max_time_size);

    float gpu_time = CUDA_Surface::compute_surface(gpu_surf, gpu_max_time, gpu_centroids, gpu_pos);
    printf("Success. gpu time = %f\n", gpu_time);

    gpu_pos.free();
    gpu_surf.free();
    gpu_centroids.free();
    gpu_max_time.free();
    printf("Memory was successfully deallocated\n");
    
    return 0;
}
