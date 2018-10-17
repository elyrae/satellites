#include "orbits.h"
#include "earth.h"
#include "messages.h"
#include "mathstuff.h"

#include <cstdio>
#include <cmath>
#include <fstream>

//bool isCorrectCircularOrbit(const Orbits::CircularOrbit& orbit)
//{
//    return (0.0 <= orbit.psi);
//}

// Чтение из текстового файла параметров круговых орбит.
Orbits::Constellation Orbits::readCircularOrbits(const std::string& filepath)
{
    std::ifstream in(filepath);
    Orbits::Constellation orbits;
    Orbits::CircularOrbit orbit;
    while (in.good()) {
        in >> orbit.ascendingNode;
        in >> orbit.inclination;
        in >> orbit.initialPhase;
        in >> orbit.height;

        orbit.ascendingNode = MathStuff::degreesToRad(orbit.ascendingNode);
        orbit.inclination   = MathStuff::degreesToRad(orbit.inclination);
        orbit.initialPhase  = MathStuff::degreesToRad(orbit.initialPhase);
        orbit.height        = orbit.height*1000.0;

        if (in.good())
            orbits.push_back(orbit);
        else break;
    }
    return orbits;
}

void Orbits::printCircularOrbits(const Orbits::Constellation& orbits)
{
    printf(Messages::circularOrbitHeader.c_str());
    for (const Orbits::CircularOrbit& orbit: orbits)
        printf(Messages::circularOrbitMessage.c_str(),
               MathStuff::radToDegrees(orbit.ascendingNode),
               MathStuff::radToDegrees(orbit.inclination),
               MathStuff::radToDegrees(orbit.initialPhase),
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

//        orbit.ascendingNode = MathStuff::degreesToRad(orbit.ascendingNode);
//        orbit.inclination = MathStuff::degreesToRad(orbit.inclination);
//        orbit.argumentOfPeriapsis = MathStuff::degreesToRad(orbit.argumentOfPeriapsis);
//        orbit.meanAnomaly = MathStuff::degreesToRad(orbit.meanAnomaly);

//        orbit.apocenter  = orbit.apocenter  * 1000.0;
//        orbit.pericenter = orbit.pericenter * 1000.0;

//        if(in.status() == QTextStream::Ok)
//            orbits.push_back(orbit);
//    }
//    return orbits;
//}

double Orbits::solveKeplersEquation(const double t, const Orbits::EllipticalOrbit &orbit, const double eps)
{
    const double semiMajorAxis = orbit.semiMajorAxis();
    const double M = orbit.meanAnomaly + t*sqrt(Earth::mu / (semiMajorAxis*semiMajorAxis*semiMajorAxis));
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

double Orbits::EllipticalOrbit::semiMajorAxis() const
{
    return Earth::radius + (apocenter + pericenter)/2.0;
}

double Orbits::EllipticalOrbit::period() const
{
    return 2.0*M_PI*sqrt(semiMajorAxis()*semiMajorAxis()*semiMajorAxis() / Earth::mu);
}

double Orbits::trueAnomaly(const double E, const double e)
{
    return 2.0*atan2( sqrt(1 + e)*sin(E / 2.0), sqrt(1 - e)*cos(E / 2.0) );
}

double Orbits::EllipticalOrbit::height(const double E) const
{
    return semiMajorAxis()*(1.0 - eccentricity()*cos(E));
}

double Orbits::CircularOrbit::semiMajorAxis() const
{
    return height + Earth::radius;
}

double Orbits::CircularOrbit::period() const
{
    return 2.0*M_PI*pow(height*1000.0 + Earth::radius, 1.5) / sqrt(Earth::mu);
}

double Orbits::CircularOrbit::meanMotion() const
{
    return sqrt(Earth::mu / (semiMajorAxis()*semiMajorAxis()*semiMajorAxis()));
}

//void Orbits::printEllipticalOrbits(const std::vector<Orbits::EllipticalOrbit> &orbits)
//{
//    QTextStream out(stdout);
//    out << "Elliptical orbits:\n";
//    for (const Orbits::EllipticalOrbit& orbit: orbits)
//        out << Messages::ellipticalOrbitMessage.arg(MathStuff::radToDegrees(orbit.ascendingNode),   -5, 'g', 4)
//                                               .arg(MathStuff::radToDegrees(orbit.inclination), -5, 'g', 4)
//                                               .arg(MathStuff::radToDegrees(orbit.argumentOfPeriapsis), -5, 'g', 4)
//                                               .arg(orbit.pericenter / 1000.0, -7, 'g', 5)
//                                               .arg(orbit.apocenter / 1000.0 , -7, 'g', 5)
//                                               .arg(MathStuff::radToDegrees(orbit.meanAnomaly), -5, 'g', 2)
//                                               .arg(orbit.period() / 60.0, -5, 'g', 2);
//}
