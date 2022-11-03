#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;

// Namespace dataio is used to read the file input and write result files

namespace dataio{

    bool getParameters(std::string&, std::string&, std::string&, int&, std::string&, int, char**);

    std::vector<Point> readPoints(const std::string& );

    void createResultsFile(const std::vector<Segment_2>&, const ft&, const std::chrono::milliseconds&, const ft& ,const std::string&, const std::string&, const int&, const std::string&);
}