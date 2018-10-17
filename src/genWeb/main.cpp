#include "grid.h"
#include "satellitesurface.h"
#include "orbits.h"
//#include "simplex.h"
#include "earth.h"
#include "messages.h"
#include "mathstuff.h"
#include "examples.h"
//#include "distances.h"

#include <QTime>
#include <QFile>
#include <QTextStream>
#include <utility>

// double mean(const std::vector<double>& v)
// {
//     double sum = 0.0;
//     for (const double elem: v)
//         sum = sum + elem;
//     return sum/(1.0*v.size());
// }

void out_histogram()
{
    // Чтение сетки
    //bool ok = false;
    //Grid::SphereGrid sphereGrid = Grid::readGrid({"gridCentr-small.txt", "gridAreas-small.txt"}, &ok);
    Grid::Centroids centroids = Grid::readCentroids("gridCentr-small.txt");
//    if (!ok)
//        return;

    // Чтение и печать орбит
    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    Orbits::printCircularOrbits(orbits);

    // Вычисление области покрытия и сохранение результата в файл для визуализации в Mathematica.
    Settings::Sets settings = Settings::readSettings("settings.ini");
    Settings::printAlgorithmSettings(settings);
    //settings.deltaT = 15.0;
    //SatelliteSurface::Surface surface = SatelliteSurface::compute(sphereGrid.centroids, orbits, parameters);

    QFile histoFile("hist_time.txt");
    histoFile.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outHist(&histoFile);

    auto timeData = Surface::computeTimeFull(centroids, orbits, settings);

    outHist.setRealNumberPrecision(6);
    for (size_t i = 0; i < timeData.size(); i++) {
        outHist << centroids.X[i] << " " << centroids.Y[i] << " " << centroids.Z[i] << " "
                << timeData[i] << "\n";
    }
    //    for (size_t i = 0; i < time_data.size(); i++) {
//        outHist << "{" << centroids.X[i] << "," << centroids.Y[i] << ","
//                << centroids.Z[i] << ","
//                << time_data[i] << /*( (i < time_data.size()-1) ? "," : "" ) <<*/ "}";
//        outHist << ((i < time_data.size() - 1) ? "," : "");
//    }
//    for (size_t i = 0; i < time_data.size(); i++) {
//        outHist << time_data[i] << /*( (i < time_data.size()-1) ? "," : "" ) <<*/ "\n";
//        //outHist << ((i < time_data.size() - 1) ? "," : "");
//    }
        //outHist << time_data[i] << ( (i < time_data.size()-1) ? "," : "" );
    //outHist << "}";
    histoFile.close();
}

int main(int, char **)
{
    Examples::swarm_optimisation();
    // Examples::example_surface_computation_circular();
    // Examples::example_surface_computation_elliptical();
    
    return 0;
}
