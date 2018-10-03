#include "settings.h"
#include "messages.h"

#include <QString>
#include <QSettings>
#include <QTextStream>

Settings::Sets Settings::readSettings(const QString &iniFile)
{
    QSettings settings(iniFile, QSettings::IniFormat);
    Settings::Sets parameters;
    const Settings::Sets defaults = Settings::defaultParameters;

    parameters.coneAngle    = settings.value("Satellites/ConeAngle", defaults.coneAngle).toDouble();
    parameters.deltaT       = settings.value("Orbits/DeltaT", 60.0).toDouble();
    parameters.timeDuration = settings.value("Orbits/TimeDuration", 60.0).toDouble() * 60.0;

    return Settings::isCorrectSets(parameters) ? parameters : defaults;
}

void Settings::printAlgorithmSettings(const Settings::Sets& params)
{
    QTextStream out(stdout);
    out << Messages::settingsMessage.arg(params.coneAngle)
                                    .arg(params.timeDuration)
                                    .arg(params.deltaT);
}

bool Settings::isCorrectSets(const Settings::Sets &params)
{
    return ((0.0 < params.coneAngle) && (params.coneAngle < 180.0) &&
            (0.0 < params.deltaT) && (0.0 < params.timeDuration));
}
