#ifndef GRID_H
#define GRID_H

#include <array>
#include <vector>
// #include <utility>

namespace Grid {
    using Point = std::array<double, 3>;
    using Cell = std::array<int, 3>;

    using Nodes = std::vector<Point>;
    using Cells = std::vector<Cell>;
    using Areas = std::vector<double>;

    struct TriangularGrid {
        Nodes nodes;
        Cells cells;
    };

    struct Centroids {
        std::vector<double> X;
        std::vector<double> Y;
        std::vector<double> Z;
    };

    // struct SphereGrid {
    //     Centroids centroids;
    //     Areas areas;
    //     double area;
    // };

    TriangularGrid generate(const size_t iterations);

    Centroids centroids(const TriangularGrid &grid);
    Areas areas(const TriangularGrid &grid);

    void write_nodes(const Nodes &nodes, const std::string &file);
    void write_cells(const Cells &cells, const std::string &file);
    void write_centroids(const Centroids &centroids, const std::string &file);
    void write_areas(const Areas &areas, const std::string &file);

    Nodes     read_nodes(const std::string& file);
    Cells     read_cells(const std::string& file);
    Centroids read_centroids(const std::string& file);
    Areas     read_areas(const std::string& file);
}

#endif // GRID_H
