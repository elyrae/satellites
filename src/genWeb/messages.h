#ifndef MESSAGES_H
#define MESSAGES_H

#include <QString>

namespace Messages {
    const QString meshGeneratedMessage = "Mesh:\n"
                                         "    %1 iterations\n"
                                         "    %2 nodes\n"
                                         "    %3 cells.\n"
                                         "Timing: %4ms.\n\n";
    const QString circularOrbitMessage = "    [%1 deg|%2 deg|%3 deg|%4 km]\n";
    const QString ellipticalOrbitMessage = "[%1 deg|%2 deg|%3 deg|%4 km|%5 km|%6 deg|%7 min]\n";
    const QString settingsMessage = "Algorithm parameters:\n"
                                    "    Cone angle: %1 degrees.\n"
                                    "    Duration: %2 seconds.\n"
                                    "    DeltaT: %3 seconds.\n";
}

#endif // MESSAGES_H
