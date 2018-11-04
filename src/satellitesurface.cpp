#include "satellitesurface.h"
#include "earth.h"
#include "messages.h"
#include "mathstuff.h"

#include <fstream>
#include <cmath>

// старая формула для вычисления горизонта без учета возвышения спутника над горизонтом
// inline double horizon(const double H, const double alpha)
// {
//     return (H*sin(alpha) < 1.0) ? cos(asin(H*sin(alpha)) - alpha) : (1.0 / H);
// }

inline double horizon(const double H, const double alpha)
{
    const double delta = MathStuff::degreesToRad(0.0); // требуемое возвышение спутника над горизонтом
    const double alpha_star = asin(cos(delta) / H);

    // return (alpha < alpha_star) ? cos(asin(H*sin(alpha)) - alpha) : sin(delta + alpha_star);
    return sin(delta + alpha_star);
}

Surface::Surf Surface::compute(const Grid::Centroids& centroids, const Orbits::Constellation& orbits, const Settings::Sets& settings)
{
    const double alpha = MathStuff::degreesToRad(settings.coneAngle) / 2.0;
    double cosNode        = 0.0, sinNode        = 0.0, px         = 0.0, py         = 0.0, pz = 0.0,
           cosInclination = 0.0, sinInclination = 0.0, cosAnomaly = 0.0, sinAnomaly = 0.0;

    Surface::Surf surface(centroids.X.size(), 0);
    double horiz = 0.0, t = 0.0, trueAnomaly = 0.0, semiMajorAxis = 0.0, meanMotion = 0.0;
    for (const Orbits::CircularOrbit& orbit: orbits) {
        semiMajorAxis = orbit.semiMajorAxis();
        meanMotion = sqrt(Earth::mu / (semiMajorAxis*semiMajorAxis*semiMajorAxis));
        horiz = horizon(semiMajorAxis / Earth::radius, alpha);

        cosInclination = cos(orbit.inclination);
        sinInclination = sin(orbit.inclination);
        t = 0.0;
        while (t < settings.timeDuration + settings.deltaT/2.0) {
            trueAnomaly = meanMotion*t;
            cosNode = cos(orbit.ascendingNode - Earth::angularVelocity*t);
            sinNode = sin(orbit.ascendingNode - Earth::angularVelocity*t);
            cosAnomaly = cos(trueAnomaly + orbit.initialPhase);
            sinAnomaly = sin(trueAnomaly + orbit.initialPhase);

            px =                cosNode*cosAnomaly - cosInclination*sinAnomaly*sinNode;
            py = cosInclination*cosNode*sinAnomaly +                cosAnomaly*sinNode;
            pz = sinInclination*sinAnomaly;
            for (size_t i = 0; i < surface.size(); i++)
                surface[i] |= ((centroids.X[i]*px + centroids.Y[i]*py + centroids.Z[i]*pz) > horiz);
            t = t + settings.deltaT;
        }
    }
    return surface;
}

double Surface::computeTime(const Grid::Centroids& centroids, const Orbits::Constellation& orbits, const Settings::Sets& settings)
{
    const double alpha = MathStuff::degreesToRad(settings.coneAngle) / 2.0;
    double cos_node = 0.0, sin_node = 0.0, cos_i = 0.0, cos_AnomalyPlusPhase = 0.0, sin_i = 0.0,
           sin_AnomalyPlusPhase = 0.0, px = 0.0, py = 0.0,         pz = 0.0;

    Surface::Surf    surface(centroids.X.size(), 0);
    std::vector<double> time(centroids.X.size(), 0.0); // Time, Dr. Freeman? Is it really that time again?
    double max_time = 0.0;

    double trueAnomaly = 0.0, t = 0.0, horiz = 0.0, semiMajorAxis = 0.0, meanAngularVelocity = 0.0;
    while (t < settings.timeDuration + settings.deltaT/2.0) {
        for (const Orbits::CircularOrbit& orbit: orbits) {
            semiMajorAxis = orbit.semiMajorAxis();
            meanAngularVelocity = sqrt(Earth::mu / (semiMajorAxis*semiMajorAxis*semiMajorAxis));
            trueAnomaly = meanAngularVelocity*t;
            horiz = horizon(semiMajorAxis / Earth::radius, alpha);

            cos_i = cos(orbit.inclination);
            sin_i = sin(orbit.inclination);
            cos_node = cos(orbit.ascendingNode - Earth::angularVelocity*t);
            sin_node = sin(orbit.ascendingNode - Earth::angularVelocity*t);
            cos_AnomalyPlusPhase = cos(trueAnomaly + orbit.initialPhase);
            sin_AnomalyPlusPhase = sin(trueAnomaly + orbit.initialPhase);

            px =       cos_node*cos_AnomalyPlusPhase - cos_i*sin_AnomalyPlusPhase*sin_node;
            py = cos_i*cos_node*sin_AnomalyPlusPhase +       cos_AnomalyPlusPhase*sin_node;
            pz = sin_i*sin_AnomalyPlusPhase;
            for (size_t i = 0; i < surface.size(); i++)
                surface[i] |= ((centroids.X[i]*px + centroids.Y[i]*py + centroids.Z[i]*pz) > horiz);
        }
        for (size_t i = 0; i < surface.size(); i++) {
            if (surface[i]) {
                if (time[i] > max_time)
                    max_time = time[i];
                time[i] = 0.0;
            }
            else time[i] += settings.deltaT;
            surface[i] = 0;
        }
        t = t + settings.deltaT;
    }
    for (size_t i = 0; i < surface.size(); i++) {
        if (time[i] > max_time)
            max_time = time[i];
    }
    return max_time;
}

Surface::Times Surface::computeTimeFull(const Grid::Centroids& centroids, const Orbits::Constellation& orbits, const Settings::Sets& settings)
{
    const double alpha = MathStuff::degreesToRad(settings.coneAngle) / 2.0;
    double cos_node = 0.0, sin_node = 0.0, cos_i = 0.0, cos_AnomalyPlusPhase = 0.0, sin_i = 0.0,
           sin_AnomalyPlusPhase = 0.0, px = 0.0, py = 0.0, pz = 0.0;

    Surface::Surf surface(centroids.X.size(), 0);
    std::vector<double>     time(centroids.X.size(), 0.0); // Time, Dr. Freeman? Is it really that time again?
    std::vector<double> max_time(centroids.X.size(), 0.0);
    double trueAnomaly = 0.0, t = 0.0, horiz = 0.0, semiMajorAxis = 0.0, meanAngularVelocity = 0.0;

    while (t < settings.timeDuration + settings.deltaT/2.0) {
        for (const Orbits::CircularOrbit& orbit: orbits) {
            semiMajorAxis = orbit.semiMajorAxis();
            meanAngularVelocity = sqrt(Earth::mu / (semiMajorAxis*semiMajorAxis*semiMajorAxis));
            trueAnomaly = meanAngularVelocity*t;
            horiz = horizon(semiMajorAxis / Earth::radius, alpha);

            cos_i = cos(orbit.inclination);
            sin_i = sin(orbit.inclination);
            cos_node = cos(orbit.ascendingNode - Earth::angularVelocity*t);
            sin_node = sin(orbit.ascendingNode - Earth::angularVelocity*t);
            cos_AnomalyPlusPhase = cos(trueAnomaly + orbit.initialPhase);
            sin_AnomalyPlusPhase = sin(trueAnomaly + orbit.initialPhase);

            px =       cos_node*cos_AnomalyPlusPhase - cos_i*sin_AnomalyPlusPhase*sin_node;
            py = cos_i*cos_node*sin_AnomalyPlusPhase +       cos_AnomalyPlusPhase*sin_node;
            pz = sin_i*sin_AnomalyPlusPhase;
            for (size_t i = 0; i < surface.size(); i++)
                surface[i] |= ((centroids.X[i]*px + centroids.Y[i]*py + centroids.Z[i]*pz) > horiz);
        }
        for (size_t i = 0; i < surface.size(); i++) {
            if (surface[i]) {
                if (time[i] > max_time[i])
                    max_time[i] = time[i];
                time[i] = 0.0;
            }
            else time[i] = time[i] + settings.deltaT;
            surface[i] = 0;
        }
        t = t + settings.deltaT;
    }
    for (size_t i = 0; i < surface.size(); i++)
        if (time[i] > max_time[i])
            max_time[i] = time[i];
    return max_time;
}

//SatelliteSurface::Surface SatelliteSurface::compute(const std::vector<Grid::Node>& centroids,
//                                                    const Orbits::EllipticalConstellation& orbits,
//                                                    const SatelliteSurface::Settings& settings)
//{
//    const double alpha = MathStuff::degreesToRad(settings.coneAngle) / 2.0;
//    double cos_omega = 0.0, sin_omega = 0.0,     cos_i = 0.0, c_anomalyPlusArgument = 0.0,
//               sin_i = 0.0, ci_comega = 0.0, ci_somega = 0.0, s_anomalyPlusArgument = 0.0;
//    double  centroid_x = 0.0,  centroid_y = 0.0,  centroid_z = 0.0,
//                    xC = 0.0,          zC = 0.0, z = 0.0;
//    SatelliteSurface::Surface surface(centroids.size(), 0);

//    double t = 0.0, H = 0.0, horiz = 0.0, trueAnomaly = 0.0, E = 0.0;
//    //QTextStream out(stdout);
//    //out << "delta t = " << settings.deltaT << "\n";
//    for (const Orbits::EllipticalOrbit& orbit: orbits) {
//        cos_i = cos(orbit.inclination);
//        sin_i = sin(orbit.inclination);
//        t = 0.0;
//        while (t < settings.timeDuration + settings.deltaT/2.0) {
//            E = Orbits::solveKeplersEquation(t, orbit, 1.0E-6);
//            trueAnomaly = Orbits::trueAnomaly(E, orbit.eccentricity());

//            H = Orbits::height(E, orbit) / Earth::radius;
//            horiz = horizon(H, alpha);
//            //out << "t = " << t << ", Hsin(a) = " << H*sin(alpha) << ", H = " << H << ", horiz = " << horiz << "\n";

//            cos_omega = cos(orbit.ascendingNode - Earth::angularVelocity*t);
//            sin_omega = sin(orbit.ascendingNode - Earth::angularVelocity*t);
//            ci_comega = cos_i*cos_omega, ci_somega = cos_i*sin_omega;
//            c_anomalyPlusArgument = cos(trueAnomaly + orbit.argumentOfPeriapsis);
//            s_anomalyPlusArgument = sin(trueAnomaly + orbit.argumentOfPeriapsis);
//            for (int i = 0; i < surface.size(); i++) {
//                if (surface[i] == 0) {
//                    centroid_x = centroids[i][0];
//                    centroid_y = centroids[i][1];
//                    centroid_z = centroids[i][2];
//                    xC = centroid_x*cos_omega + centroid_y*sin_omega;
//                    zC = centroid_y*ci_comega + centroid_z*sin_i - centroid_x*ci_somega;

//                    z = xC*c_anomalyPlusArgument + zC*s_anomalyPlusArgument;
//                    if (z > horiz)
//                        surface[i] = 1;
//                }
//            }
//            t = t + settings.deltaT;
//        }
//    }
//    return surface;
//}

//double SatelliteSurface::computeTime(const QVector<Grid::Node>& centroids,
//                                     const QVector<Orbits::EllipticalOrbit>& orbits,
//                                     const SatelliteSurface::Settings& settings)
//{
//    const double alpha = MathStuff::degreesToRad(settings.coneAngle) / 2.0;
//    double cos_omega = 0.0, sin_omega = 0.0,     cos_i = 0.0, c_gammaPlusPhase = 0.0,
//               sin_i = 0.0, ci_comega = 0.0, ci_somega = 0.0, s_gammaPlusPhase = 0.0;
//    double  centroid_x = 0.0,  centroid_y = 0.0,  centroid_z = 0.0,
//                    xC = 0.0,          zC = 0.0, z = 0.0;
//    double maxT = 0.0;
//    SatelliteSurface::Surface surface(centroids.size(), 0);
//    QVector<double> time(centroids.size(), 0.0); // Time, Dr. Freeman? Is it really that time again?

//    double t = 0.0, H = 0.0, horiz = 0.0;
//    double trueAnomaly = 0.0, E = 0.0;
//    //QTextStream out(stdout);
//    while (t < settings.timeDuration + settings.deltaT/2.0) {
//        for (const Orbits::EllipticalOrbit& orbit: orbits) {
//            E = Orbits::solveKeplersEquation(t, orbit, 1.0E-6);
//            trueAnomaly = Orbits::trueAnomaly(E, orbit.eccentricity());

//            H = Orbits::height(E, orbit) / Earth::radius;
//            horiz = horizon(H, alpha);

//            cos_i = cos(orbit.inclination);
//            sin_i = sin(orbit.inclination);
//            cos_omega = cos(orbit.ascendingNode - Earth::angularVelocity*t);
//            sin_omega = sin(orbit.ascendingNode - Earth::angularVelocity*t);
//            ci_comega = cos_i*cos_omega, ci_somega = cos_i*sin_omega;
//            c_gammaPlusPhase = cos(trueAnomaly + orbit.argumentOfPeriapsis);
//            s_gammaPlusPhase = sin(trueAnomaly + orbit.argumentOfPeriapsis);
//            for (int i = 0; i < surface.size(); i++) {
//                if (surface[i] == '\0') {
//                    centroid_x = centroids[i][0];
//                    centroid_y = centroids[i][1];
//                    centroid_z = centroids[i][2];
//                    xC = centroid_x*cos_omega + centroid_y*sin_omega;
//                    zC = centroid_y*ci_comega + centroid_z*sin_i - centroid_x*ci_somega;

//                    z = xC*c_gammaPlusPhase + zC*s_gammaPlusPhase;
//                    if (z > horiz)
//                        surface[i] = 1;
//                }
//            }
//        }
//        for (int i = 0; i < surface.size(); i++) {
//            if (surface[i] == '\1') {
//                if (time[i] > maxT)
//                    maxT = time[i];

//                time[i] = 0.0;
//            }
//            else time[i] = time[i] + settings.deltaT;
//            surface[i] = 0;
//        }
//        t = t + settings.deltaT;
//    }
//    return maxT;
//}

double Surface::sumArea(const Grid::Areas& areas, const Surface::Surf& surface)
{
    double area = 0.0;
    for (size_t i = 0; i < surface.size(); i++)
        if (surface[i])
            area = area + areas[i];
    return area;
}

void Surface::writeToTextFile(const Surface::Surf& surface, const std::string& filepath)
{
    std::ofstream out(filepath);
    for (const char surfElement: surface)
        out << (char)(surfElement + '0') << " ";
}
