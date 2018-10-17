#ifndef ORBITS_H
#define ORBITS_H

#include <string>
#include <vector>
//#include <QString>

// Конечно же, круговая орбита является частным случаем эллиптической, но
// такое разделение удобно, так как алгоритмы счета немного разные
enum class OrbitType {Circular, Elliptical};

namespace Orbits {
    struct CircularOrbit {
        double ascendingNode; // долгота восходящего узла, рад
        double inclination;   // наклонение, рад
        double initialPhase;  // начальная фаза, рад

        double height; // высота орбиты над уровнем моря, м

        double semiMajorAxis() const; // большая полуось = радиус Земли + высота орбиты
        double period()        const; // период, c
        double meanMotion()    const; // средняя угловая скорость, рад/с
    };

    struct EllipticalOrbit {
        double ascendingNode;       // Долгота восходящего узла, рад
        double inclination;         // Наклонение, рад
        double argumentOfPeriapsis; // Аргумент перицентра, рад

        double apocenter;           // Апоцентр, м
        double pericenter;          // Перицентр, м

        double meanAnomaly;         // Средняя аномалия, рад

        double eccentricity()  const; // эксцентриситет
        double semiMajorAxis() const; // большая полуось
        double period()        const; // период обращения

        double height(const double E) const; // высота спутника над уровнем моря от экс. аномалии
    };

    using Constellation = std::vector<CircularOrbit>;
    using EllipticalConstellation = std::vector<EllipticalOrbit>;

    double solveKeplersEquation(const double t, const EllipticalOrbit& orbit, const double eps = 1.0E-4);
    double trueAnomaly(const double E, const double e);

    // Constellation readCircularOrbits(const QString& filepath);
    Constellation readCircularOrbits(const std::string& filepath);
    void printCircularOrbits(const Constellation& orbits);

    // EllipticalConstellation readEllipticalOrbits(const QString& filepath);
    // void printEllipticalOrbits(const EllipticalConstellation& orbits);
}

#endif // ORBITS_H
