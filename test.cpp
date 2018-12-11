/**
 * Build with the command
 * g++ --std=c++11 -I /path/to/your/eigen/libary -O3 -o testdt test.cpp
 */

#include "delaunay-trianglation.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <random>

using namespace Eigen;

/*
 * visualization command
 gnuplot -p -e "set size ratio -1; plot 'edges.dat' u 1:2:(\$3-\$1):(\$4-\$2) w vector notitle, 'points.dat' u 1:2:0 w labels offset 0, char 1 notitle"
 */
void test(int n)
{
    int seed = std::time(nullptr);

//    std::cout << "seed: " << seed << std::endl;

    std::srand(seed);
    Matrix<double, 2, Dynamic> points(2, n);
    points.setRandom();

    std::ofstream("points.dat") << points.transpose();

    Delaunay<double> dt(points);
    dt.triangulate();

    Matrix<int, 2, Dynamic> edges = dt.get_edges();

    int num_edges = edges.cols();
    MatrixXd out(4, num_edges);

    for (int i = 0; i < num_edges; ++i) {
        out.col(i).head(2) = points.col(edges(0, i));
        out.col(i).tail(2) = points.col(edges(1, i));
    }

    std::ofstream("edges.dat") << out.transpose();

}

int main(int argc, char ** argv)
{
    int n = 10;
    if (argc == 2) {
        try {
            n = std::stoi(argv[1]);
        }
        catch (std::exception & e) {
            std::cout << "invalid parameter" << std::endl;
            return 1;
        }

        if (n < 3) {
            std::cout << "number of points is too small" << std::endl;
            return 0;
        }
    }

    test(n);
    return 0;
}
