#include "orbits.h"
#include "earth.h"
#include "messages.h"
#include "stuff.h"

#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

// Чтение текстового файла параметров круговых орбит
Orbits::Constellation Orbits::read_circular_orbits(const std::string &filepath)
{
    std::ifstream in(filepath);
    Orbits::Constellation orbits;
    Orbits::CircularOrbit orbit;
    while (in.good()) {
        in >> orbit.ascending_node;
        in >> orbit.inclination;
        in >> orbit.initial_phase;
        in >> orbit.height;

        orbit.ascending_node = Stuff::degrees_to_rad(orbit.ascending_node);
        orbit.inclination    = Stuff::degrees_to_rad(orbit.inclination);
        orbit.initial_phase  = Stuff::degrees_to_rad(orbit.initial_phase);
        orbit.height         = orbit.height*1000.0;

        orbits.push_back(orbit);
    }
    return orbits;
}

void Orbits::print_circular_orbits(const Orbits::Constellation &orbits)
{
    printf("%s", Messages::circular_orbit_header.c_str());
    for (const Orbits::CircularOrbit &orbit: orbits)
        printf(Messages::circular_orbit_message.c_str(),
               Stuff::rad_to_degrees(orbit.ascending_node),
               Stuff::rad_to_degrees(orbit.inclination),
               Stuff::rad_to_degrees(orbit.initial_phase),
               orbit.height / 1000.0);
}

double Orbits::solve_keplers_equation(const double t, const Orbits::EllipticalOrbit &orbit, const double eps)
{
    const double semi_major_axis = orbit.semi_major_axis();
    const double M = orbit.mean_anomaly + t*sqrt(Earth::mu / (semi_major_axis*semi_major_axis*semi_major_axis));
    const int iterations = 5;

    double E = M, delta = 0.0;
    int iteration = 0;
    do {
        delta = (E - orbit.eccentricity()*sin(E) - M) / (1.0 - orbit.eccentricity()*cos(E));
        E = E - delta;

        iteration++;
    } while ((fabs(delta) > eps) & (iteration < iterations));
    return E;
}

double Orbits::EllipticalOrbit::eccentricity() const
{
    return (apocenter - pericenter) / (apocenter + pericenter + 2.0*Earth::radius);
}

double Orbits::EllipticalOrbit::semi_major_axis() const
{
    return Earth::radius + (apocenter + pericenter)/2.0;
}

double Orbits::EllipticalOrbit::period() const
{
    return 2.0*M_PI*sqrt(semi_major_axis()*semi_major_axis()*semi_major_axis() / Earth::mu);
}

double Orbits::true_anomaly(const double E, const double e)
{
    return 2.0*atan2(sqrt(1 + e)*sin(E / 2.0), sqrt(1 - e)*cos(E / 2.0));
}

double Orbits::EllipticalOrbit::height(const double E) const
{
    return semi_major_axis()*(1.0 - eccentricity()*cos(E));
}

double Orbits::CircularOrbit::semi_major_axis() const
{
    return height + Earth::radius;
}

double Orbits::CircularOrbit::period() const
{
    return 2.0*M_PI*pow(height*1000.0 + Earth::radius, 1.5) / sqrt(Earth::mu);
}

double Orbits::CircularOrbit::mean_motion() const
{
    return sqrt(Earth::mu / (semi_major_axis()*semi_major_axis()*semi_major_axis()));
}
