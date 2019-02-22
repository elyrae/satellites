#ifndef EXAMPLES_H
#define EXAMPLES_H

namespace Examples {
    void grid_generation(const int iterations);
    
    // Пример единичного вычисления максимального времени ожидания 
    void time_computation();

    // Вычисление макс. ожидания точек сетки по отдельности для
    // высот орбит в 1500 и 2000 км и возвышениями спутников в 0 и 10 градусов над горизонтом. 
    void time_histogram();

    // Сравнение макс. времени ожидания для разного размера сеток
    void compare_different_grids();
    
    // Что случится с временем ожидания, если спутники в группировке будут выходить из строя?
    void destroy_satellites();

    // Вычисление зоны покрытия спутниковой группировкой
    void surface_computation();

    // Оптимизация времени ожидания методом роя частиц.  
    void swarm_optimisation();
}

#endif // EXAMPLES_H
