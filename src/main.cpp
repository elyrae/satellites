#include "grid.h"
#include "satellitesurface.h"
#include "orbits.h"
#include "earth.h"
#include "messages.h"
#include "mathstuff.h"
#include "examples.h"

#include <fstream>
#include <iostream>

void writeTimeHistogram()
{
    // Чтение сетки
    Grid::Centroids centroids = Grid::readCentroids("gridCentroidsNew.txt");

    // Чтение и печать орбит
    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    Orbits::printCircularOrbits(orbits);

    // Вычисление области покрытия и сохранение результата в файл для визуализации в Mathematica.
    Settings::Sets settings = Settings::readSettings("settings.ini");
    Settings::printSettings(settings);

    auto timeData = Surface::computeTimeFull(centroids, orbits, settings);
    std::ofstream outHist("timeHistogram.txt");
    for (size_t i = 0; i < timeData.size(); i++) {
        //outHist << centroids.X[i] << " " << centroids.Y[i] << " " << centroids.Z[i] << " "
        outHist << timeData[i] << "\n";
    }
}

int main(int, char **)
{
    // Examples::exampleGridGeneration(4);

    // writeTimeHistogram();
    
    Examples::SwarmOptimisation();

    // Examples::example_surface_computation_circular();
    // Examples::example_surface_computation_elliptical();
    
    return 0;
}
