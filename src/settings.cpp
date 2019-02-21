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
        return Settings::defaults;
    }

    parameters.cone_angle      = reader.GetReal("Satellites", "ConeAngle",      120.0);
    parameters.elevation_angle = reader.GetReal("Satellites", "ElevationAngle",   0.0);
    
    parameters.time_duration = reader.GetReal("Orbits", "TimeDuration", 60.0) * 60.0;
    parameters.delta_t       = reader.GetReal("Orbits", "DeltaT", 60.0);
    return Settings::is_correct_sets(parameters) ? parameters : Settings::defaults;
}

void Settings::print_settings(const Settings::Sets& par)
{
    printf(Messages::settings_message.c_str(), par.cone_angle, par.elevation_angle, par.time_duration, par.delta_t);
}

bool Settings::is_correct_sets(const Settings::Sets &params)
{
    return ((0.0 < params.cone_angle)       && (params.cone_angle      < 180.0) &&
            !(params.elevation_angle < 0.0) && (params.elevation_angle <  90.0) &&
            (0.0 < params.delta_t)          && (0.0 < params.time_duration));
}
