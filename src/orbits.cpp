#include "orbits.h"
#include "earth.h"
#include "messages.h"
#include "mathstuff.h"

#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

//bool isCorrectCircularOrbit(const Orbits::CircularOrbit& orbit)
//{
//    return (0.0 <= orbit.psi);
//}

// Чтение из текстового файла параметров круговых орбит.
Orbits::Constellation Orbits::read_circular_orbits(const std::string& filepath)
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

        //if (in.good())
            orbits.push_back(orbit);
        //else break;
    }
    return orbits;
}

void Orbits::print_circular_orbits(const Orbits::Constellation& orbits)
{
    // нехорошо, конечно, смешивать std::cout << и printf
    std::cout << Messages::circular_orbit_header;
    for (const Orbits::CircularOrbit& orbit: orbits)
        printf(Messages::circular_orbit_message.c_str(),
               Stuff::rad_to_degrees(orbit.ascending_node),
               Stuff::rad_to_degrees(orbit.inclination),
               Stuff::rad_to_degrees(orbit.initial_phase),
               orbit.height / 1000.0);
}

//Orbits::EllipticalConstellation Orbits::readEllipticalOrbits(const QString &filepath)
//{
//    Orbits::EllipticalConstellation orbits;
//    QFile orbitsFile(filepath);
//    if (!orbitsFile.open(QIODevice::ReadOnly | QIODevice::Text))
//        return {};

//    QTextStream in(&orbitsFile);
//    Orbits::EllipticalOrbit orbit;
//    while (!in.atEnd()) {
//        in >> orbit.ascendingNode;
//        in >> orbit.inclination;
//        in >> orbit.argumentOfPeriapsis;

//        in >> orbit.pericenter;
//        in >> orbit.apocenter;

//        in >> orbit.meanAnomaly;

//        orbit.ascendingNode = Stuff::degreesToRad(orbit.ascendingNode);
//        orbit.inclination = Stuff::degreesToRad(orbit.inclination);
//        orbit.argumentOfPeriapsis = Stuff::degreesToRad(orbit.argumentOfPeriapsis);
//        orbit.meanAnomaly = Stuff::degreesToRad(orbit.meanAnomaly);

//        orbit.apocenter  = orbit.apocenter  * 1000.0;
//        orbit.pericenter = orbit.pericenter * 1000.0;

//        if(in.status() == QTextStream::Ok)
//            orbits.push_back(orbit);
//    }
//    return orbits;
//}

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
    } while ((fabs(delta) > eps) && (iteration < iterations));
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

//void Orbits::printEllipticalOrbits(const std::vector<Orbits::EllipticalOrbit> &orbits)
//{
//    QTextStream out(stdout);
//    out << "Elliptical orbits:\n";
//    for (const Orbits::EllipticalOrbit& orbit: orbits)
//        out << Messages::ellipticalOrbitMessage.arg(Stuff::radToDegrees(orbit.ascendingNode),   -5, 'g', 4)
//                                               .arg(Stuff::radToDegrees(orbit.inclination), -5, 'g', 4)
//                                               .arg(Stuff::radToDegrees(orbit.argumentOfPeriapsis), -5, 'g', 4)
//                                               .arg(orbit.pericenter / 1000.0, -7, 'g', 5)
//                                               .arg(orbit.apocenter / 1000.0 , -7, 'g', 5)
//                                               .arg(Stuff::radToDegrees(orbit.meanAnomaly), -5, 'g', 2)
//                                               .arg(orbit.period() / 60.0, -5, 'g', 2);
//}
