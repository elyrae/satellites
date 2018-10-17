#ifndef GRID_H
#define GRID_H

#include <array>
#include <vector>
#include <utility>

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

    struct SphereGrid {
        std::vector<Point> centroids;
        std::vector<double> areas;
        double area;
    };

    struct Centroids {
        std::vector<double> X;
        std::vector<double> Y;
        std::vector<double> Z;
    };

    TriangularGrid generate(const int iterations);

    std::vector<Point> centroids(const Grid::TriangularGrid& grid);
    Areas areas(const Grid::TriangularGrid& grid);

    void writeNodes(const std::vector<Point>& nodes, const std::string& file);
    void writeCells(const std::vector<Cell>& cells, const std::string& file);
    void writeCentroids(const std::vector<Point>& centroids, const std::string& file);
    void writeAreas(const std::vector<double>& areas, const std::string& file);

    Nodes readNodes(const std::string& file);
    Cells readCells(const std::string& file);
    Centroids readCentroids(const std::string& file);
    Areas readAreas( const std::string& file);

}

#endif // GRID_H
