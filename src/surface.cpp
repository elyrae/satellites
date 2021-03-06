#include "surface.h"
#include "earth.h"
#include "messages.h"
#include "stuff.h"
#include "grid.h"

#include <fstream>
#include <cmath>
#include <algorithm>
// #include <iostream>

using Surface::Timegrid;

Timegrid::Timegrid(const size_t grid_iterations, const Orbits::Constellation &def_cons, const Settings::Sets &_sets)
{
    def_constellation = def_cons;

    sets = _sets;
    sets_for_area = _sets;
    sets_for_area.delta_t = area_delta_t;

    Grid::TriangularGrid grid = Grid::generate(grid_iterations);
    centroids = Grid::centroids(grid);
    areas = Grid::areas(grid);

    area = 0.0;
    for (const double a: areas)
        area += a;
}

double Timegrid::compute_max_time(const Orbits::Constellation &orb) const
{
    return Surface::compute_time(centroids, orb, sets);
}

double Timegrid::compute_covered_area(const Orbits::Constellation &orb) const
{
    return Surface::compute_area(centroids, areas, orb, sets_for_area);
}

// старая формула для вычисления горизонта без учета возвышения спутника над горизонтом
double horizon(const double H, const double alpha)
{
    return (H*sin(alpha) < 1.0) ? cos(asin(H*sin(alpha)) - alpha) : (1.0 / H);
}

double horizon(const double H, const double alpha, const double delta)
{
    const double alpha_star = asin(cos(delta) / H);
    return (alpha < alpha_star) ? cos(asin(H*sin(alpha)) - alpha) : sin(delta + alpha_star);
}

void Surface::fill_orbits(const Orbits::ConstellationUnion &configurations, const Settings::Sets &settings, 
                          std::vector<float> &x, std::vector<float> &y, std::vector<float> &z, std::vector<float> &h)
{
    const int TIMESTEPS = 240;

    const double alpha = Stuff::degrees_to_rad(settings.cone_angle) / 2.0;
    const size_t configs = configurations.size();
    for (size_t iconf = 0; iconf < configs; ++iconf) {
        Orbits::Constellation orbits = configurations[iconf];

        for (size_t iorb = 0; iorb < orbits.size(); ++iorb) {
            double semi_major_axis = orbits[iorb].semi_major_axis();
            double mean_angular_velocity = sqrt(Earth::mu / (semi_major_axis*semi_major_axis*semi_major_axis));
            double cos_i = cos(orbits[iorb].inclination);
            double sin_i = sin(orbits[iorb].inclination);

            for (size_t timestep = 0; timestep < TIMESTEPS; ++timestep) {
                double t = timestep*settings.delta_t;
                double cos_node = cos(orbits[iorb].ascending_node - Earth::angular_velocity*t);
                double sin_node = sin(orbits[iorb].ascending_node - Earth::angular_velocity*t);
                double cos_anomaly = cos(mean_angular_velocity*t + orbits[iorb].initial_phase);
                double sin_anomaly = sin(mean_angular_velocity*t + orbits[iorb].initial_phase);

                x[(timestep*orbits.size() + iorb)*configs + iconf] = (float) (      cos_node*cos_anomaly - cos_i*sin_anomaly*sin_node); 
                y[(timestep*orbits.size() + iorb)*configs + iconf] = (float) (cos_i*cos_node*sin_anomaly +       cos_anomaly*sin_node); 
                z[(timestep*orbits.size() + iorb)*configs + iconf] = (float) (sin_i*sin_anomaly);
                h[(timestep*orbits.size() + iorb)*configs + iconf] = (float) horizon(semi_major_axis / Earth::radius, alpha);        
            }
        }
    }
}

Surface::Surf Surface::compute(const Grid::Centroids &centroids, const Orbits::Constellation &orbits, const Settings::Sets &settings)
{
    const double alpha = Stuff::degrees_to_rad(settings.cone_angle) / 2.0;
    const double elevation = Stuff::degrees_to_rad(settings.elevation_angle);
    double cos_node        = 0.0, sin_node        = 0.0, x               = 0.0, y           = 0.0, z     = 0.0,
           cos_inclination = 0.0, sin_inclination = 0.0, cos_anomaly     = 0.0, sin_anomaly = 0.0, horiz = 0.0, 
           t               = 0.0, semi_major_axis = 0.0, mean_motion = 0.0;
    Surface::Surf surface(centroids.X.size(), 0);

    for (const Orbits::CircularOrbit &orbit: orbits) {
        semi_major_axis = orbit.semi_major_axis();
        mean_motion = sqrt(Earth::mu / (semi_major_axis*semi_major_axis*semi_major_axis));
        horiz = horizon(semi_major_axis / Earth::radius, alpha, elevation);

        cos_inclination = cos(orbit.inclination);
        sin_inclination = sin(orbit.inclination);
        t = 0.0;
        while (t < settings.time_duration + settings.delta_t/2.0) {
            cos_node = cos(orbit.ascending_node - Earth::angular_velocity*t);
            sin_node = sin(orbit.ascending_node - Earth::angular_velocity*t);
            cos_anomaly = cos(mean_motion*t + orbit.initial_phase);
            sin_anomaly = sin(mean_motion*t + orbit.initial_phase);

            x =                 cos_node*cos_anomaly - cos_inclination*sin_anomaly*sin_node;
            y = cos_inclination*cos_node*sin_anomaly +                 cos_anomaly*sin_node;
            z = sin_inclination*sin_anomaly;
            for (size_t i = 0; i < surface.size(); i++)
                surface[i] |= ((centroids.X[i]*x + centroids.Y[i]*y + centroids.Z[i]*z) > horiz);
            t = t + settings.delta_t;
        }
    }
    return surface;
}

double Surface::compute_area(const Grid::Centroids &centroids, const Grid::Areas &areas, 
                             const Orbits::Constellation &orbits, const Settings::Sets &settings)
{
    const double alpha = Stuff::degrees_to_rad(settings.cone_angle) / 2.0;
    const double elevation = Stuff::degrees_to_rad(settings.elevation_angle);
    double cos_node        = 0.0, sin_node        = 0.0, x           = 0.0, y           = 0.0, z     = 0.0,
           cos_inclination = 0.0, sin_inclination = 0.0, cos_anomaly = 0.0, sin_anomaly = 0.0, horiz = 0.0, 
           t               = 0.0, semi_major_axis = 0.0, mean_motion = 0.0;
    Surface::Surf surface(centroids.X.size(), 0);

    for (const Orbits::CircularOrbit &orbit: orbits) {
        semi_major_axis = orbit.semi_major_axis();
        mean_motion = sqrt(Earth::mu / (semi_major_axis*semi_major_axis*semi_major_axis));
        horiz = horizon(semi_major_axis / Earth::radius, alpha, elevation);

        cos_inclination = cos(orbit.inclination);
        sin_inclination = sin(orbit.inclination);
        t = 0.0;
        while (t < settings.time_duration + settings.delta_t/2.0) {
            cos_node = cos(orbit.ascending_node - Earth::angular_velocity*t);
            sin_node = sin(orbit.ascending_node - Earth::angular_velocity*t);
            cos_anomaly = cos(mean_motion*t + orbit.initial_phase);
            sin_anomaly = sin(mean_motion*t + orbit.initial_phase);

            x =                 cos_node*cos_anomaly - cos_inclination*sin_anomaly*sin_node;
            y = cos_inclination*cos_node*sin_anomaly +                 cos_anomaly*sin_node;
            z = sin_inclination*sin_anomaly;
            for (size_t i = 0; i < surface.size(); i++)
                surface[i] |= ((centroids.X[i]*x + centroids.Y[i]*y + centroids.Z[i]*z) > horiz);
            t = t + settings.delta_t;
        }
    }

    double area = 0.0;
    for (size_t i = 0; i < surface.size(); i++)
        area += surface[i] ? areas[i] : 0.0;
    return area;    
}

// double Surface::compute_time_parallel(const Grid::Centroids &centroids, const Orbits::Constellation &orbits, const Settings::Sets &settings)
// {
//     const double alpha = Stuff::degrees_to_rad(settings.cone_angle) / 2.0;

//     Surface::Surf    surface(centroids.X.size(), 0);
//     std::vector<double> time(centroids.X.size(), 0.0); // Time, Dr. Freeman? Is it really that time again?

//     double max_time = 0.0;
//     size_t timesteps = floor(settings.time_duration / settings.delta_t);
//     std::vector<double>   x(orbits.size()*timesteps, 0.0);
//     std::vector<double>   y(orbits.size()*timesteps, 0.0);
//     std::vector<double>   z(orbits.size()*timesteps, 0.0);
//     std::vector<double> hor(orbits.size()*timesteps, 0.0);

//     const int THREADS = 4;
//     #pragma omp parallel num_threads(THREADS)
//     {
//         #pragma omp for
//         for (size_t i = 0; i < orbits.size(); ++i) {
//             double semi_major_axis = orbits[i].semi_major_axis();
//             double mean_angular_velocity = sqrt(Earth::mu / (semi_major_axis*semi_major_axis*semi_major_axis));        
//             double cos_i = cos(orbits[i].inclination);
//             double sin_i = sin(orbits[i].inclination);        
        
//             for (size_t j = 0; j < timesteps; ++j) {
//                 double t = j*settings.delta_t;

//                 double cos_node = cos(orbits[i].ascending_node - Earth::angular_velocity*t);
//                 double sin_node = sin(orbits[i].ascending_node - Earth::angular_velocity*t);
//                 double cos_anomaly = cos(mean_angular_velocity*t + orbits[i].initial_phase);
//                 double sin_anomaly = sin(mean_angular_velocity*t + orbits[i].initial_phase);

//                   x[j*orbits.size() + i] =       cos_node*cos_anomaly - cos_i*sin_anomaly*sin_node;
//                   y[j*orbits.size() + i] = cos_i*cos_node*sin_anomaly +       cos_anomaly*sin_node;
//                   z[j*orbits.size() + i] = sin_i*sin_anomaly;
//                 hor[j*orbits.size() + i] = Surface::horizon(semi_major_axis / Earth::radius, alpha);            
//             }
//         }

//         for (size_t i = 0; i < timesteps; ++i) { // timestep
//             #pragma omp for
//             for (size_t k = 0; k < surface.size(); ++k) // point              
//             for (size_t j = 0; j <  orbits.size(); ++j) // orbit           
//                 surface[k] |= ((centroids.X[k]*x[i*orbits.size() + j] 
//                               + centroids.Y[k]*y[i*orbits.size() + j] 
//                               + centroids.Z[k]*z[i*orbits.size() + j]) > hor[i*orbits.size() + j]);
               
//             #pragma omp for reduction(max : max_time)                 
//             for (size_t i = 0; i < surface.size(); i++) {
//                 if (surface[i]) {
//                     max_time = std::max(time[i], max_time);
//                     time[i] = 0.0;
//                 }
//                 else time[i] += settings.delta_t;
//                 surface[i] = 0;
//             }          
//         }

//         #pragma omp for reduction(max : max_time)
//         for (size_t i = 0; i < surface.size(); i++)
//             max_time = std::max(max_time, time[i]);        
//     }
//     return max_time;
// }

void satellite_position(const Orbits::CircularOrbit &orbit, const double t, double *x, double *y, double *z)
{
    double cos_node = 0.0, sin_node = 0.0, cos_incl = 0.0, sin_incl = 0.0, cos_anomaly = 0.0, sin_anomaly = 0.0;
    double semi_major_axis = orbit.semi_major_axis();
    double mean_angular_velocity = sqrt(Earth::mu / (semi_major_axis*semi_major_axis*semi_major_axis));
    sincos(orbit.inclination,                                &sin_incl ,   &cos_incl);
    sincos(orbit.ascending_node - Earth::angular_velocity*t, &sin_node,    &cos_node);
    sincos(mean_angular_velocity*t + orbit.initial_phase   , &sin_anomaly, &cos_anomaly);

    *x = cos_node*cos_anomaly - cos_incl*sin_anomaly*sin_node;
    *y = sin_node*cos_anomaly + cos_incl*sin_anomaly*cos_node;
    *z = sin_incl*sin_anomaly;
}

double Surface::compute_time(const Grid::Centroids &centroids, const Orbits::Constellation &orbits, const Settings::Sets &settings)
{
    const double alpha     = Stuff::degrees_to_rad(settings.cone_angle) / 2.0;
    const double elevation = Stuff::degrees_to_rad(settings.elevation_angle);
    double cos_node = 0.0, sin_node        = 0.0, x                     = 0.0, y           = 0.0, z     = 0.0,
              cos_i = 0.0, sin_i           = 0.0, cos_anomaly           = 0.0, sin_anomaly = 0.0, horiz = 0.0, 
                  t = 0.0, semi_major_axis = 0.0, mean_angular_velocity = 0.0;

    Surface::Surf    surface(centroids.X.size(), 0);
    std::vector<double> time(centroids.X.size(), 0.0); // Time, Dr. Freeman? Is it really that time again?
    
    double max_time = 0.0;
    while (t < settings.time_duration + settings.delta_t/2.0) {
        for (const Orbits::CircularOrbit &orbit: orbits) {
            semi_major_axis = orbit.semi_major_axis();
            mean_angular_velocity = sqrt(Earth::mu / (semi_major_axis*semi_major_axis*semi_major_axis));
            horiz = horizon(semi_major_axis / Earth::radius, alpha, elevation);

            cos_i = cos(orbit.inclination);
            sin_i = sin(orbit.inclination);
            cos_node = cos(orbit.ascending_node - Earth::angular_velocity*t);
            sin_node = sin(orbit.ascending_node - Earth::angular_velocity*t);
            cos_anomaly = cos(mean_angular_velocity*t + orbit.initial_phase);
            sin_anomaly = sin(mean_angular_velocity*t + orbit.initial_phase);

            x =       cos_node*cos_anomaly - cos_i*sin_anomaly*sin_node;
            y = cos_i*cos_node*sin_anomaly +       cos_anomaly*sin_node;
            z = sin_i*sin_anomaly;
            for (size_t i = 0; i < surface.size(); i++)
                surface[i] |= ((centroids.X[i]*x + centroids.Y[i]*y + centroids.Z[i]*z) > horiz);
        }
        for (size_t i = 0; i < surface.size(); i++) {
            if (surface[i]) {
                max_time = std::max(time[i], max_time);
                time[i] = 0.0;
            }
            else time[i] += settings.delta_t;
            surface[i] = 0;
        }            
        t += settings.delta_t;
    }
    for (size_t i = 0; i < surface.size(); i++)
        max_time = std::max(max_time, time[i]);
    return max_time;
}

Surface::Times Surface::compute_time_full(const Grid::Centroids &centroids, const Orbits::Constellation &orbits, const Settings::Sets &settings)
{
    const double alpha     = Stuff::degrees_to_rad(settings.cone_angle) / 2.0;
    const double elevation = Stuff::degrees_to_rad(settings.elevation_angle);
    double cos_node = 0.0, sin_node = 0.0, cos_incl = 0.0, cos_anomaly  = 0.0, sin_incl = 0.0, sin_anomaly = 0.0, 
           x        = 0.0, y        = 0.0, z     = 0.0, t            = 0.0, horiz = 0.0, 
           semi_major_axis = 0.0, mean_angular_velocity = 0.0;

    Surface::Surf surface(centroids.X.size(), 0);
    std::vector<double>     time(centroids.X.size(), 0.0); // Time, Dr. Freeman? Is it really that time again?
    std::vector<double> max_time(centroids.X.size(), 0.0);

    while (t < settings.time_duration + settings.delta_t/2.0) {
        for (const Orbits::CircularOrbit &orbit: orbits) {
            semi_major_axis = orbit.semi_major_axis();
            mean_angular_velocity = sqrt(Earth::mu / (semi_major_axis*semi_major_axis*semi_major_axis));
            horiz = horizon(orbit.semi_major_axis() / Earth::radius, alpha, elevation);

            cos_incl = cos(orbit.inclination);
            sin_incl = sin(orbit.inclination);
            cos_node = cos(orbit.ascending_node - Earth::angular_velocity*t);
            sin_node = sin(orbit.ascending_node - Earth::angular_velocity*t);
            cos_anomaly = cos(mean_angular_velocity*t + orbit.initial_phase);
            sin_anomaly = sin(mean_angular_velocity*t + orbit.initial_phase);

            x = cos_node*cos_anomaly - cos_incl*sin_anomaly*sin_node;
            y = sin_node*cos_anomaly + cos_incl*sin_anomaly*cos_node;
            z = sin_incl*sin_anomaly;
            for (size_t i = 0; i < surface.size(); i++)
                surface[i] |= ((centroids.X[i]*x + centroids.Y[i]*y + centroids.Z[i]*z) > horiz);
        }
        for (size_t i = 0; i < surface.size(); i++) {
            if (surface[i]) {
                max_time[i] = std::max(time[i], max_time[i]);
                time[i] = 0.0;
            }
            else time[i] += settings.delta_t;
            surface[i] = 0;
        }
        t += settings.delta_t;
    }
    for (size_t i = 0; i < surface.size(); i++)
        max_time[i] = std::max(time[i], max_time[i]);
    return max_time;
}

void Surface::write(const Surface::Surf &surface, const std::string &filepath)
{
    std::ofstream out(filepath);
    for (const char surf_element: surface)
        out << (char)(surf_element + '0') << "\n";
}

void Surface::write(const Times &time, const std::string &filepath)
{
    std::ofstream out(filepath);
    for (const double t: time)
        out << t << "\n";    
}