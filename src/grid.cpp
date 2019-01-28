#include "grid.h"

#include <fstream>
#include <cmath>

const int out_precision = 12;

const Grid::Cells initial_cells =
{{1, 11, 7}, {1,  7, 6}, {1,  6, 10}, {1, 10, 3},
 {1, 3, 11}, {4,  8, 0}, {5,  4,  0}, {9, 5,  0},
 {2, 9, 0 }, {8,  2, 0}, {11, 9,  7}, {7, 2,  6},
 {6, 8, 10}, {10, 4, 3}, {3,  5, 11}, {4, 10, 8},
 {5, 3, 4 }, {9, 11, 5}, {2,  7,  9}, {8, 6,  2}};

const Grid::Nodes initial_nodes =
{{0                  ,                   0, -0.9510565162951536},
 {0                  ,                   0,  0.9510565162951536},
 {-0.8506508083520399,                   0, -0.42532540417602  },
 {0.85065080835203990,                   0,  0.42532540417602  },
 {0.68819096023558680,                -0.5, -0.42532540417602  },
 {0.68819096023558680,                 0.5, -0.42532540417602  },
 {-0.6881909602355868,                -0.5,  0.42532540417602  },
 {-0.6881909602355868,                 0.5,  0.42532540417602  },
 {-0.2628655560595668, -0.8090169943749474, -0.42532540417602  },
 {-0.2628655560595668,  0.8090169943749474, -0.42532540417602  },
 {0.26286555605956680, -0.8090169943749474,  0.42532540417602  },
 {0.26286555605956680,  0.8090169943749474,  0.42532540417602  }};

// Середина отрезка
Grid::Point center(const Grid::Point &first, const Grid::Point &second)
{
    return {(first[0] + second[0]) / 2.0, (first[1] + second[1]) / 2.0, (first[2] + second[2]) / 2.0};
}

// Нормирование всех точек на единицу
void normalize(std::vector<Grid::Point> &points)
{
    for(Grid::Point &p: points)
        for (int i = 0; i < 3; i++)
            p[i] /= sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

void normalize(Grid::Centroids &p)
{
    double norm = 0.0;
    for (size_t i = 0; i < p.X.size(); i++) {
        norm = sqrt(p.X[i]*p.X[i] + p.Y[i]*p.Y[i] + p.Z[i]*p.Z[i]);

        p.X[i] /= norm;
        p.Y[i] /= norm;
        p.Z[i] /= norm;
    }
}

// Создание треугольной сетки с количеством итераций iterations
Grid::TriangularGrid Grid::generate(const int iterations)
{
    Grid::TriangularGrid grid;
    grid.nodes = initial_nodes;
    grid.cells = initial_cells;

    Grid::Cells new_cells;

    int k = 0;
    for (int it = 1; it <= iterations; it++)
    {
        new_cells.clear();
        for (const Grid::Cell &cell: grid.cells) {
            k = grid.nodes.size();

            new_cells.push_back({cell[0], k,       k + 2  });
            new_cells.push_back({k,       cell[1], k + 1  });
            new_cells.push_back({k + 2,   k + 1,   cell[2]});
            new_cells.push_back({k,       k + 1,   k + 2  });

            grid.nodes.push_back(center(grid.nodes[cell[0]], grid.nodes[cell[1]]));
            grid.nodes.push_back(center(grid.nodes[cell[1]], grid.nodes[cell[2]]));
            grid.nodes.push_back(center(grid.nodes[cell[0]], grid.nodes[cell[2]]));
        }
        grid.cells = new_cells;
    }
    normalize(grid.nodes);
    return grid;
}

// Получение центральных точек всех треугольных граней.
Grid::Centroids Grid::centroids(const Grid::TriangularGrid &grid)
{
    Grid::Centroids centroids;
    centroids.X.reserve(grid.cells.size());
    centroids.Y.reserve(grid.cells.size());
    centroids.Z.reserve(grid.cells.size());

    for (const Grid::Cell &cell: grid.cells) {
        centroids.X.push_back((grid.nodes[cell[0]][0] + grid.nodes[cell[1]][0] + grid.nodes[cell[2]][0]) / 3.0);
        centroids.Y.push_back((grid.nodes[cell[0]][1] + grid.nodes[cell[1]][1] + grid.nodes[cell[2]][1]) / 3.0);
        centroids.Z.push_back((grid.nodes[cell[0]][2] + grid.nodes[cell[1]][2] + grid.nodes[cell[2]][2]) / 3.0);
    }
    normalize(centroids);

    return centroids;
}

// Вычисление площадей всех треугольных граней сетки
Grid::Areas Grid::areas(const Grid::TriangularGrid &grid)
{
    double ab = 0.0, bc = 0.0, ac = 0.0, p = 0.0;
    Grid::Point a, b, c;
    Grid::Areas areas;
    areas.reserve(grid.cells.size());

    for (const Grid::Cell &cell: grid.cells) {
        ab = 0.0; bc = 0.0; ac = 0.0;
        a = grid.nodes[cell[0]]; b = grid.nodes[cell[1]]; c = grid.nodes[cell[2]];
        for (int i = 0; i < 3; i++) {
            ab += (b[i] - a[i])*(b[i] - a[i]);
            bc += (c[i] - b[i])*(c[i] - b[i]);
            ac += (c[i] - a[i])*(c[i] - a[i]);
        }
        ab = std::sqrt(ab); 
        bc = std::sqrt(bc); 
        ac = std::sqrt(ac);

        p = (ab + bc + ac) / 2.0;
        areas.push_back(std::sqrt( p*(p - ab)*(p - bc)*(p - ac) ));
    }
    return areas;
}

void Grid::write_centroids(const Grid::Centroids &centroids, const std::string &file)
{
    std::ofstream out(file);
    out.precision(out_precision);
    for (size_t i = 0; i < centroids.X.size(); i++)
        out << centroids.X[i] << " " << centroids.Y[i] << " " << centroids.Z[i] << "\n";
}

void Grid::write_areas(const Grid::Areas &areas, const std::string &file)
{
    std::ofstream out(file);
    out.precision(out_precision);
    for (const double area: areas)
        out << area << "\n";
}

void Grid::write_nodes(const Grid::Nodes &nodes, const std::string &file)
{
    std::ofstream out(file);
    out.precision(out_precision);
    for (const Grid::Point &node: nodes)
        out << node[0] << " " << node[1] << " " << node[2] << "\n";
}

void Grid::write_cells(const Grid::Cells &cells, const std::string &file)
{
    std::ofstream out(file);
    for (const Grid::Cell &cell: cells)
        out << cell[0] << " " << cell[1] << " " << cell[2] << "\n";
}

Grid::Nodes Grid::read_nodes(const std::string &file)
{
    std::ifstream in(file);
    Grid::Nodes nodes;

    double x = 0.0, y = 0.0, z = 0.0;
    while (in.good()) {
        in >> x;
        in >> y;
        in >> z;
            
        nodes.push_back({x, y, z});
    }
    return nodes;
}

Grid::Cells Grid::read_cells(const std::string &file)
{
    std::ifstream in(file);
    Grid::Cells cells;

    int m = 0, n = 0, k = 0;
    while (in.good()) {
        in >> m;
        in >> n;
        in >> k;
            
        cells.push_back({m, n, k});
    }
    return cells;
}

Grid::Centroids Grid::read_centroids(const std::string &file)
{
    std::ifstream in(file);
    Grid::Centroids centroids;

    double x = 0.0, y = 0.0, z = 0.0;
    while (in.good()) {
        in >> x;
        in >> y;
        in >> z;
            
        centroids.X.push_back(x);
        centroids.Y.push_back(y);
        centroids.Z.push_back(z);
    }
    return centroids;
}

Grid::Areas Grid::read_areas(const std::string &file)
{
    std::ifstream in(file);
    Grid::Areas areas;

    double area = 0.0;
    while (in.good()) {
        in >> area;
        areas.push_back(area);
    }
    return areas;
}
