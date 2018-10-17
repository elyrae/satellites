#ifndef MESSAGES_H
#define MESSAGES_H

#include <string>

namespace Messages {
    // const std::string meshGeneratedMessage = "Mesh:\n"
    //                                      "    %1 iterations\n"
    //                                      "    %2 nodes\n"
    //                                      "    %3 cells.\n"
    //                                      "Timing: %4ms.\n\n";
    const std::string circularOrbitHeader  = "[Ascending node|Inclination|Initial phase|    Height]\n";
    const std::string circularOrbitMessage = "[%10.2f deg|%7.2f deg|%9.2f deg|%8.2f km]\n";
    const std::string ellipticalOrbitMessage = "[%f deg|%f deg|%f deg|%f km|%f km|%f deg]\n";
    const std::string settingsMessage = "Algorithm parameters:\n"
                                        "    Satellite cone angle: %f degrees.\n"
                                        "    Duration: %f seconds.\n"
                                        "    Time step: %f seconds.\n";
}

#endif // MESSAGES_H
