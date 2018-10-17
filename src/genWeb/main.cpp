#include "grid.h"
#include "satellitesurface.h"
#include "orbits.h"
#include "earth.h"
#include "messages.h"
#include "mathstuff.h"
#include "examples.h"

#include <QTime>
#include <fstream>
#include <utility>

void out_histogram()
{
    // Чтение сетки
    //Grid::SphereGrid sphereGrid = Grid::readGrid({"gridCentr-small.txt", "gridAreas-small.txt"}, &ok);
    Grid::Centroids centroids = Grid::readCentroids("gridCentr-small.txt");

    // Чтение и печать орбит
    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    Orbits::printCircularOrbits(orbits);

    // Вычисление области покрытия и сохранение результата в файл для визуализации в Mathematica.
    Settings::Sets settings = Settings::readSettings("settings.ini");
    Settings::printSettings(settings);
    //settings.deltaT = 15.0;
    //SatelliteSurface::Surface surface = SatelliteSurface::compute(sphereGrid.centroids, orbits, parameters);

    auto timeData = Surface::computeTimeFull(centroids, orbits, settings);

    std::ofstream outHist("hist_time.txt");
    for (size_t i = 0; i < timeData.size(); i++) {
        //outHist << centroids.X[i] << " " << centroids.Y[i] << " " << centroids.Z[i] << " "
        outHist << timeData[i] << "\n";
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
    outHist.close();
}

int main(int, char **)
{
    Examples::exampleGridGeneration(4);

    auto orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    Orbits::printCircularOrbits(orbits);

    // Examples::swarm_optimisation();
    // Examples::example_surface_computation_circular();
    // Examples::example_surface_computation_elliptical();
    
    return 0;
}
