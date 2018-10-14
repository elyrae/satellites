#ifndef GRID_H
#define GRID_H

#include <array>
#include <vector>
#include <utility>
#include <QString>

namespace Grid {
    const int reserveLimit = 100000;

    using Point = std::array<double, 3>;
    using Cell = std::array<int, 3>;
    using Areas = std::vector<double>;

    struct TriangularGrid {
        std::vector<Point> nodes;
        std::vector<Cell> cells;
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
    std::vector<double> areas(const Grid::TriangularGrid& grid);

    void writeSphereGrid(const SphereGrid& sphereGrid, const std::pair<QString, QString>& files);
    void writeTriangularGrid(const Grid::TriangularGrid& triangGrid,
                             const std::pair<QString, QString>& files);
    SphereGrid readGrid(const std::pair<QString, QString>& files, bool *ok);
    Centroids  readCentroids(const QString& file, bool *ok);
    Areas      readAreas(const QString& file, bool *ok);
}

#endif // GRID_H
