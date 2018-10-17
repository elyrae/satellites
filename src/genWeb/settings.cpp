#include "settings.h"
#include "messages.h"
#include "INIReader.h"

#include <cstdio>
#include <iostream>

Settings::Sets Settings::readSettings(const std::string& iniFile)
{
    Settings::Sets parameters;
    INIReader reader(iniFile);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load '" << iniFile << "'\n";
        return Settings::defaultParameters;
    }

    parameters.coneAngle    = reader.GetReal("Satellites", "ConeAngle", 120.0);
    parameters.deltaT       = reader.GetReal("Orbits", "DeltaT", 60.0);
    parameters.timeDuration = reader.GetReal("Orbits", "TimeDuration", 60.0) * 60.0;
    return Settings::isCorrectSets(parameters) ? parameters : Settings::defaultParameters;
}

void Settings::printSettings(const Settings::Sets& params)
{
    printf(Messages::settingsMessage.c_str(), params.coneAngle, params.timeDuration, params.deltaT);
}

bool Settings::isCorrectSets(const Settings::Sets &params)
{
    return ((0.0 < params.coneAngle) && (params.coneAngle < 180.0) &&
            (0.0 < params.deltaT) && (0.0 < params.timeDuration));
}
