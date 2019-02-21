#ifndef MESSAGES_H
#define MESSAGES_H

#include <string>

namespace Messages {
    const std::string circular_orbit_header  = "[Ascending node|Inclination|Initial phase|     Height]\n";
    const std::string circular_orbit_message = "[%10.2f deg|%7.2f deg|%9.2f deg|%8.2f km]\n";
    const std::string elliptical_orbit_message = "[%f deg|%f deg|%f deg|%f km|%f km|%f deg]\n";
    const std::string settings_message = "Algorithm parameters:\n"
                                         "    Satellites cone angle: %.2f degrees.\n"
                                         "    Satellites elevation angle: %.2f degrees.\n"
                                         "    Duration: %.2f seconds.\n"
                                         "    Time step: %.2f seconds.\n";
}

#endif // MESSAGES_H
