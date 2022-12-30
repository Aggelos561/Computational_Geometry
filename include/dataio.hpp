#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;

// Namespace dataio is used to read the file input and write result files

namespace dataio{

    bool getParameters(std::string& nameOfDirectory, std::string& outputFile, bool& preprocessEnabled, int argc, char** argv);

    std::vector<Point> readPoints(const std::string& );

    std::vector<std::pair<int, std::string>> findFiles(const std::string& );

    void writeToOutputFile(const std::string&, const std::vector<std::pair<double, double>>&, const std::vector<double>&, const std::vector<double>&, int, bool);
}