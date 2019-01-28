#include "settings.h"
#include "messages.h"
#include "INIReader.h"

#include <cstdio>
#include <iostream>

Settings::Sets Settings::read_settings(const std::string& ini_file)
{
    Settings::Sets parameters;
    INIReader reader(ini_file);
    if (reader.ParseError() < 0) {
        std::cout << "Can't load '" << ini_file << "'\n";
        return Settings::default_parameters;
    }

    parameters.cone_angle    = reader.GetReal("Satellites", "ConeAngle", 120.0);
    parameters.delta_t       = reader.GetReal("Orbits", "DeltaT", 60.0);
    parameters.time_duration = reader.GetReal("Orbits", "TimeDuration", 60.0) * 60.0;
    return Settings::is_correct_sets(parameters) ? parameters : Settings::default_parameters;
}

void Settings::print_settings(const Settings::Sets& params)
{
    printf(Messages::settings_message.c_str(), params.cone_angle, params.time_duration, params.delta_t);
}

bool Settings::is_correct_sets(const Settings::Sets &params)
{
    return ((0.0 < params.cone_angle) && (params.cone_angle < 180.0) &&
            (0.0 < params.delta_t) && (0.0 < params.time_duration));
}
