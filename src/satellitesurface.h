#ifndef SATELLITESURFACE_H
#define SATELLITESURFACE_H

#include "grid.h"
#include "orbits.h"
#include "settings.h"

#include <string>

namespace Surface {
    // Каждому треугольнику сетки ставим в соответствие флаг "покрыт связью"
    // Использование аж целого инта связано с заметно большей скоростью счета
    // Подумаю еще, как сделать получше
    using Surf = std::vector<int>;
    using Times = std::vector<double>;

    // нахождение области покрытия группировкой
    Surf compute(const Grid::Centroids& centroids,
                 const Orbits::Constellation& orbits,
                 const Settings::Sets& parameters = Settings::defaultParameters);
//    Surf compute(const std::vector<Grid::Node>& centroids,
//                 const Orbits::EllipticalConstellation& orbits,
//                 const Settings& settings);

    // поиск максимального времени ожидания
    double computeTime(const Grid::Centroids& centroids,
                       const Orbits::Constellation& orbits,
                       const Settings::Sets& parameters = Settings::defaultParameters);
//    double computeTime(const std::vector<Grid::Point>& centroids,
//                       const Orbits::Constellation& orbits,
//                       const Settings::Sets& parameters = Settings::defaultParameters);

    // максимальное время ожидания для каждой точки сетки
    Times computeTimeFull(const Grid::Centroids& centroids,
                          const Orbits::Constellation& orbits,
                          const Settings::Sets& parameters = Settings::defaultParameters);

    double sumArea(const Grid::Areas& areas, const Surf& surface);

    void writeToTextFile(const Surf& surface, const std::string& filepath);
}

#endif // SATELLITESURFACE_H
