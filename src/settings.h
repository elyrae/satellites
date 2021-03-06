#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>

namespace Settings {
    struct Sets {
        double cone_angle;      // Угол раствора конуса связи при спутниках
        double elevation_angle; // Угол возвышения спутников над горизонтом

        double time_duration; // Время счета, с
        double delta_t;       // Шаг расчета, с
    };

    const Sets defaults = { .cone_angle = 120.0, .elevation_angle = 0.0, .time_duration = 3600.0, .delta_t = 120.0};

    Sets read_settings(const std::string &ini_file);
    void print_settings(const Sets &params);

    bool is_correct_sets(const Sets &params);
}

#endif // SETTINGS_H
