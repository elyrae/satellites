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

    class Timegrid {
        Grid::Centroids centroids;
        Grid::Areas areas;
        Settings::Sets sets;

        Orbits::Constellation orb;

    public:
        Timegrid(const size_t cons_size, const std::string &centroids_filepath, 
                                         const std::string &areas_filepath,
                                         const std::string &settings_filepath);
        double compute_max_time();
        double compute_covered_area();

        Orbits::Constellation &orbits();
    };

    inline double horizon(const double H, const double alpha);
    inline double horizon(const double H, const double alpha, const double delta);

    void fill_orbits(const std::vector<Orbits::Constellation> &configurations, 
                     const Settings::Sets &settings, 
                     std::vector<float> &x, std::vector<float> &y, std::vector<float> &z, std::vector<float> &h);

    // нахождение области покрытия группировкой
    Surf compute(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                 const Settings::Sets &parameters = Settings::defaults);

    // поиск максимального времени ожидания
    double compute_time(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                        const Settings::Sets &parameters = Settings::defaults);
    double compute_timeOMP(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                           const Settings::Sets &parameters = Settings::defaults);

    // максимальное время ожидания для каждой точки сетки
    Times compute_time_full(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                            const Settings::Sets &parameters = Settings::defaults);

    double sum_area(const Grid::Areas &areas, const Surf &surface);
    void write(const Surf &surface, const std::string &filepath);
}

#endif // SATELLITESURFACE_H
