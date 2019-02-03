#ifndef ORBITS_H
#define ORBITS_H

#include <string>
#include <vector>

// Конечно же, круговая орбита является частным случаем эллиптической, но
// такое разделение удобно, так как алгоритмы счета немного разные
enum class OrbitType {Circular, Elliptical};

namespace Orbits {
    struct CircularOrbit {
        double ascending_node; // долгота восходящего узла, рад
        double inclination;   // наклонение, рад
        double initial_phase;  // начальная фаза, рад

        double height; // высота орбиты над уровнем моря, м

        double semi_major_axis() const; // большая полуось = радиус Земли + высота орбиты
        double period()          const; // период орбиты, c
        double mean_motion()     const; // средняя угловая скорость, рад/с
    };

    struct EllipticalOrbit {
        double ascending_node;        // Долгота восходящего узла, рад
        double inclination;           // Наклонение, рад
        double argument_of_periapsis; // Аргумент перицентра, рад

        double apocenter;  // Апоцентр, м
        double pericenter; // Перицентр, м

        double mean_anomaly; // Средняя аномалия, рад

        double eccentricity()    const; // эксцентриситет
        double semi_major_axis() const; // большая полуось
        double period()          const; // период обращения

        double height(const double E) const; // высота спутника над уровнем моря от экс. аномалии
    };

    using Constellation  = std::vector<CircularOrbit>;
    using ConstellationUnion = std::vector<Constellation>;

    using EllipticalConstellation = std::vector<EllipticalOrbit>;

    double solve_keplers_equation(const double t, const EllipticalOrbit& orbit, const double eps = 1.0E-4);
    double true_anomaly(const double E, const double e);

    Constellation read_circular_orbits(const std::string& filepath);
    void print_circular_orbits(const Constellation& orbits);

    // EllipticalConstellation readEllipticalOrbits(const QString& filepath);
    // void printEllipticalOrbits(const EllipticalConstellation& orbits);
}

#endif // ORBITS_H
