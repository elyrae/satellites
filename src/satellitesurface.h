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

    inline double horizon(const double H, const double alpha);
    inline double horizon(const double H, const double alpha, const double delta);

    void fill_orbits(const std::vector<Orbits::Constellation> &configurations, 
                     const Settings::Sets &settings, 
                     std::vector<float> &x, std::vector<float> &y, std::vector<float> &z, std::vector<float> &h);

    // нахождение области покрытия группировкой
    Surf compute(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                 const Settings::Sets &parameters = Settings::default_parameters);

    // поиск максимального времени ожидания
    double compute_time(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                        const Settings::Sets &parameters = Settings::default_parameters);

    double compute_timeOMP(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                           const Settings::Sets &parameters = Settings::default_parameters);

    // максимальное время ожидания для каждой точки сетки
    Times compute_time_full(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                            const Settings::Sets &parameters = Settings::default_parameters);

    double sum_area(const Grid::Areas &areas, const Surf &surface);
    void write(const Surf &surface, const std::string &filepath);
}

#endif // SATELLITESURFACE_H
