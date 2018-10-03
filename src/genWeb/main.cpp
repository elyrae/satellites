#include "grid.h"
#include "satellitesurface.h"
#include "orbits.h"
//#include "simplex.h"
#include "earth.h"
#include "messages.h"
#include "mathstuff.h"
#include "examples.h"
//#include "distances.h"

#include <QTime>
#include <QFile>
#include <QTextStream>
#include <utility>

double mean(const std::vector<double>& v)
{
    double sum = 0.0;
    for (const double elem: v)
        sum = sum + elem;
    return sum/(1.0*v.size());
}

//void plot2d_mathematica()
//{
//    // Чтение сетки
//    bool ok = false;
//    Grid::SphereGrid sphereGrid = Grid::readGrid({"gridCentr-small.txt", "gridAreas-small.txt"}, &ok);
//    if (!ok)
//        return;

//    // Чтение и печать орбит
//    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
//    Orbits::printCircularOrbits(orbits);

//    // Вычисление области покрытия и сохранение результата в файл для визуализации в Mathematica.
//    Settings::Sets parameters = Settings::readSettings("settings.ini");
//    //Settings::Sets parametersArea = parameters;
//    //SatelliteSurface::Surface surface = SatelliteSurface::compute(sphereGrid.centroids, orbits, parameters);

//    const QString pointF = "{{%1,%2},%3}";
////    double psi_1 = -M_PI*0.0, psi_2 = -M_PI*0.0, delta = M_PI / 40.0;
////    const double bound_1 = M_PI*2.0, bound_2 = M_PI*2.0;
//            double psi_1 = 0.0*M_PI, psi_2   = 0.0*M_PI, delta = M_PI / 20.0;
//    const double bound_1 = 2.0*M_PI, bound_2 = 2.0*M_PI;

//    QTextStream out(stdout);
//    QFile plotFile("plot2D.txt");
//    plotFile.open(QIODevice::WriteOnly | QIODevice::Text);
//    QTextStream outPlot(&plotFile);
//    QTime timer;

//    outPlot << "{";
//    timer.start();

//    std::vector<double> times;

//    double area = 0.0, f = 0.0;
//    while (psi_1 < bound_1 /*2.0*M_PI + delta/2.0*/) {
//        psi_2 = 0.0;
//        while (psi_2 < bound_2 /*2.0*M_PI + delta/2.0*/) {
////            orbits[1].ascendingNode = psi_1;
////            orbits[2].ascendingNode = psi_2;
////            orbits[1].initialPhase = psi_1;
////            orbits[2].initialPhase = psi_2;
////            orbits[1].ascendingNode = psi_1;
////            orbits[2].ascendingNode = psi_2;

////            orbits[0].ascendingNode = psi_1;
////            orbits[0].initialPhase  = psi_2;

////            orbits[0].inclination = psi_1;
////            orbits[1].inclination = psi_1;
////            orbits[2].inclination = psi_1;

////            orbits[3].inclination = psi_2;
////            orbits[4].inclination = psi_2;
////            orbits[5].inclination = psi_2;
//            orbits[0].ascendingNode = psi_1;
//            orbits[1].ascendingNode = psi_1;
//            orbits[2].ascendingNode = psi_1;

//            orbits[3].ascendingNode = psi_2;
//            orbits[4].ascendingNode = psi_2;
//            orbits[5].ascendingNode = psi_2;
//            timer.restart();
//            //area = SatelliteSurface::computeArea(sphereGrid, orbits, parametersArea) / ((4.0*M_PI));
//            f = Surface::computeTime(sphereGrid.centroids, orbits, parameters);
//            //+ 1.0E6*(1.0 - area)*(1.0 - area);

//            outPlot << pointF.arg(psi_1).arg(psi_2)
//                   //.arg(SatelliteSurface::computeArea(sphereGrid, orbits, parameters) / ((4.0*M_PI)));
//                   //.arg(SatelliteSurface::computeTime(sphereGrid.centroids, orbits, parameters));
//                   .arg(f);
//                   //.arg(Dist::computeTime(sphereGrid.centroids, orbits));
//                   //      .arg(SatelliteSurface::computeArea(sphereGrid.centroids, sphereGrid.,
//                   //                                         orbits, parameters) / (4.0*M_PI));
//            outPlot << ",";
//            //outPlot << ((psi_2 + delta < bound_2) ? "," : "");
//            times.push_back(timer.elapsed()*1.0);
//            out << timer.elapsed() << "(" << mean(times) << ")" << ", ";

//            out.flush();
//            psi_2 += delta;
//        }
//        //out << "\n"; out.flush();
//        psi_1 += delta;
//    }
//    outPlot << "}";
//    plotFile.close();
//    out << "\n";
//}

//void plot2d_mathematica_area()
//{
//    // Чтение сетки
//    bool ok = false;
//    Grid::SphereGrid sphereGrid = Grid::readGrid({"gridCentroids.txt", "gridAreas.txt"}, &ok);
//    if (!ok)
//        return;

//    // Чтение и печать орбит
//    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
//    Orbits::printCircularOrbits(orbits);

//    // Вычисление области покрытия и сохранение результата в файл для визуализации в Mathematica.
//    Settings::Sets parameters = Settings::readSettings("settings.ini");
//    //SatelliteSurface::Surface surface = SatelliteSurface::compute(sphereGrid.centroids, orbits, parameters);

//    const QString pointF = "{{%1,%2},%3},";
////    double psi_1 = -M_PI*0.0, psi_2 = -M_PI*0.0, delta = M_PI / 40.0;
////    const double bound_1 = M_PI*2.0, bound_2 = M_PI*2.0;
//    double psi_1 = 0.0, psi_2 = 0.0, delta = M_PI / 20.0;
//    const double bound_1 = 2.0*M_PI, bound_2 = 2.0*M_PI;

//    QTextStream out(stdout);
//    QFile plotFile("plot2D_area.txt");
//    plotFile.open(QIODevice::WriteOnly | QIODevice::Text);
//    QTextStream outPlot(&plotFile);
//    QTime timer;

//    outPlot << "{";
//    timer.start();
//    while (psi_1 < bound_1 /*2.0*M_PI + delta/2.0*/) {
//        psi_2 = 0.0;
//        while (psi_2 < bound_2 /*2.0*M_PI + delta/2.0*/) {
////            orbits[1].ascendingNode = psi_1;
////            orbits[2].ascendingNode = psi_2;
//            orbits[0].ascendingNode = psi_1;
//            orbits[1].ascendingNode = psi_1;
//            orbits[2].ascendingNode = psi_1;

//            orbits[3].ascendingNode = psi_2;
//            orbits[4].ascendingNode = psi_2;
//            orbits[5].ascendingNode = psi_2;
////            orbits[0].inclination = psi_1;
////            orbits[1].inclination = psi_1;
////            orbits[2].inclination = psi_1;

////            orbits[3].inclination = psi_2;
////            orbits[4].inclination = psi_2;
////            orbits[5].inclination = psi_2;
//            outPlot << pointF.arg(psi_1).arg(psi_2)
//                   .arg(Surface::computeArea(sphereGrid, orbits, parameters) / ((4.0*M_PI)));
//            out << timer.elapsed() << "ms ";
//            out.flush();
//            psi_2 += delta;
//        }
//        //out << "\n"; out.flush();
//        psi_1 += delta;
//    }
//    outPlot << "}";
//    plotFile.close();
//    out << "\n";
//}

//double test_int_sphere(const Orbits::Constellation& orbits)
//{
//    double t = 0;
//    double minim = Dist::nintegrateSphere(orbits, t), newint = 0.0;
//    while (t < 3600.0) {
//        t = t + 100.0;
//        newint = Dist::nintegrateSphere(orbits, t);
//        if (newint < minim)
//            minim = newint;
//    }
//    return minim;
////    out << "\n";
////    out << timer.elapsed() << "ms: " << i << ".\n";
//    //        out << "{" << t << "," << Dist::nintegrateSphere(orbits, t) << "},";
//    //        out.flush();
//    //out << QString("%1\n").arg(Dist::nintegrateCommunicationArea(orbits, {1.0, 0.0, 0.0}, 0.0),
//    //                           -15, 'g', 10);
//}

//void plot_int()
//{
//    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
//    Orbits::printCircularOrbits(orbits);

//    const QString pointF = "{{%1,%2},%3},";
//    double psi_1 = -M_PI*0.0, psi_2 = -M_PI*0.0, delta = M_PI / 20.0;
//    const double bound_1 = M_PI*2.0, bound_2 = M_PI*2.0;

//    QFile plotFile("plot2D_int.txt");
//    plotFile.open(QIODevice::WriteOnly | QIODevice::Text);
//    QTextStream outPlot(&plotFile);
//    QTextStream out(stdout);
//    outPlot << "{";
//    const double T = 120.0*60.0;
//    while (psi_1 < bound_1 /*2.0*M_PI + delta/2.0*/) {
//        psi_2 = -M_PI*0.0;
//        while (psi_2 < bound_2 /*2.0*M_PI + delta/2.0*/) {
//            orbits[1].ascendingNode = psi_1;
//            orbits[2].ascendingNode = psi_2;
////            orbits[0].ascendingNode = psi_1;
////            orbits[1].ascendingNode = psi_1;
////            orbits[2].ascendingNode = psi_1;

////            orbits[3].ascendingNode = psi_2;
////            orbits[4].ascendingNode = psi_2;
////            orbits[5].ascendingNode = psi_2;

////            orbits[1].initialPhase = psi_1;
////            orbits[2].initialPhase = psi_2;
////            orbits[1].ascendingNode = psi_1;
////            orbits[2].ascendingNode = psi_2;
////            orbits[0].ascendingNode = psi_1;
////            orbits[0].initialPhase  = psi_2;
////            orbits[0].inclination = psi_1;
////            orbits[1].inclination = psi_1;
////            orbits[2].inclination = psi_1;

////            orbits[3].inclination = psi_2;
////            orbits[4].inclination = psi_2;
////            orbits[5].inclination = psi_2;
////            orbits[0].ascendingNode = psi_1;
////            orbits[1].ascendingNode = psi_1;
////            //orbits[2].ascendingNode = psi_1;

////            orbits[2].ascendingNode = psi_2;
////            orbits[3].ascendingNode = psi_2;
//            //orbits[5].ascendingNode = psi_2;
//            outPlot << pointF.arg(psi_1).arg(psi_2).arg(test_int_sphere(orbits));
//            out << "+";
//            out.flush();
////            if (!((psi_1 + delta > bound_1) && (psi_2 + delta > bound_2)))
////                outPlot << ",";
//            psi_2 += delta;
//        }
//        //out << "\n"; out.flush();
//        psi_1 += delta;
//    }
//    //char c;
//    //outPlot >> c;
//    outPlot << "}";
//    plotFile.close();
//}



//void test_int_sphere()
//{
//    const Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
//    Orbits::printCircularOrbits(orbits);

//    QTextStream out(stdout);

//    double t = 0;
//    out << "{";
//    QTime timer;
//    timer.start();
//    double minim = Dist::nintegrateSphere(orbits, t);
//    while (t < 3600.0) {
//        t = t + 100.0;
////        out << "{" << t << "," << Dist::nintegrateSphere(orbits, t) << "},";
////        out.flush();
//        t = t + 100.0;
//    }
//    out << "\n";
//    out << timer.elapsed() << "ms: " << i << ".\n";
//    //out << QString("%1\n").arg(Dist::nintegrateCommunicationArea(orbits, {1.0, 0.0, 0.0}, 0.0),
//    //                           -15, 'g', 10);
//}

//int iterations(int argc, char *argv[])
//{
//    int it = 0;
//    if (argc > 1)
//        it = QString(argv[1]).toInt();
//    if (!((0 < it) && (it < 7)))
//        it = 5;
//    return it;
//}

//void new_time()
//{
//    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
//    Orbits::printCircularOrbits(orbits);

//    bool ok = false;
//    Grid::SphereGrid sphereGrid = Grid::readSphereGrid({"bigGridCentroids.txt", "bigGridAreas.txt"}, &ok);
//    if (!ok)
//        return;

//    QTextStream out(stdout);
//    out << Dist::computeTime(sphereGrid.centroids, orbits);
//}

//void gen()
//{
//    Grid::TriangularGrid grid = Grid::generate(1);
//    Grid::writeSphereGrid()
//}

void out_histogram()
{
    // Чтение сетки
    bool ok = false;
    Grid::SphereGrid sphereGrid = Grid::readGrid({"bigGridCentroids.txt", "bigGridAreas.txt"}, &ok);
    if (!ok)
        return;

    // Чтение и печать орбит
    Orbits::Constellation orbits = Orbits::readCircularOrbits("circularOrbits.txt");
    Orbits::printCircularOrbits(orbits);

    // Вычисление области покрытия и сохранение результата в файл для визуализации в Mathematica.
    Settings::Sets settings = Settings::readSettings("settings.ini");
    settings.deltaT = 15.0;
    //SatelliteSurface::Surface surface = SatelliteSurface::compute(sphereGrid.centroids, orbits, parameters);

    QFile histoFile("hist_time.txt");
    histoFile.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream outHist(&histoFile);

    auto time_data = Surface::computeTimeFull(sphereGrid.centroids, orbits, settings);

    outHist << "{";
    for (size_t i = 0; i < time_data.size(); i++) {
        outHist << "{" << sphereGrid.centroids[i][0] << "," << sphereGrid.centroids[i][1] << ","
                << sphereGrid.centroids[i][2] << ","
                << time_data[i] << /*( (i < time_data.size()-1) ? "," : "" ) <<*/ "}";
        outHist << ((i < time_data.size() - 1) ? "," : "");
    }
        //outHist << time_data[i] << ( (i < time_data.size()-1) ? "," : "" );
    outHist << "}";
    histoFile.close();
}

int main(int, char **)
{
    //out_histogram();
    //new_time();

    Examples::example_grid_generation();

    //test_int_sphere();
    //plot_int();
    //plot2d_mathematica();
    //plot2d_mathematica_area();
    //Examples::swarm_optimisation();
    Examples::example_surface_computation_circular();
    //Examples::example_optimization();
    //Examples::example_surface_computation_elliptical();
    //example_surface_computation_circular();
    return 0;
}
