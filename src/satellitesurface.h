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
        double area;

        Settings::Sets sets;
        Settings::Sets sets_for_area;
        Orbits::Constellation def_constellation;
    public:
        Timegrid(const size_t grid_iterations, const Orbits::Constellation &def_cons, const Settings::Sets &_sets);
        double compute_max_time(const Orbits::Constellation &orb)     const;
        double compute_covered_area(const Orbits::Constellation &orb) const;

        const Orbits::Constellation &default_constellation() const { return def_constellation; };
        const Settings::Sets &settings()                     const { return sets; };
        double full_area()                                   const { return area; };
    };

    void fill_orbits(const std::vector<Orbits::Constellation> &configurations, 
                     const Settings::Sets &settings, 
                     std::vector<float> &x, std::vector<float> &y, std::vector<float> &z, std::vector<float> &h);

    // нахождение области покрытия группировкой
    Surf compute(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                 const Settings::Sets &parameters = Settings::defaults);
    double compute_area(const Grid::Centroids &centroids, const Grid::Areas &areas, const Orbits::Constellation &orbits,
                        const Settings::Sets &parameters = Settings::defaults);

    // поиск максимального времени ожидания
    double compute_time(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                        const Settings::Sets &parameters = Settings::defaults);
    // double compute_timeOMP(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
    //                        const Settings::Sets &parameters = Settings::defaults);

    // максимальное время ожидания для каждой точки сетки
    Times compute_time_full(const Grid::Centroids &centroids, const Orbits::Constellation &orbits,
                            const Settings::Sets &parameters = Settings::defaults);

    void write(const Surf &surface, const std::string &filepath);
    void write(const Times &time,   const std::string &filepath);
}

#endif // SATELLITESURFACE_H
