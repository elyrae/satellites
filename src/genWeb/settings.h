#ifndef SETTINGS_H
#define SETTINGS_H

#include <QString>

namespace Settings {
    struct Sets {
        double coneAngle;    // Угол раствора конуса связи при спутнике
        double timeDuration; // Время счета
        double deltaT;       // Шаг расчета
    };

    const Sets defaultParameters = { .coneAngle = 120.0,
                                         .timeDuration = 3600.0,
                                         .deltaT = 120.0 };

    Sets readSettings(const QString& iniFile);
    void printAlgorithmSettings(const Sets& params);
    bool isCorrectSets(const Sets& params);
}

#endif // SETTINGS_H
