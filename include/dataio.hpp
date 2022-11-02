#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;


namespace dataio{

    std::vector<Point> readPoints();

    void createResultsFile(const std::vector<Segment_2>&, const ft&, const std::chrono::milliseconds&, const ft& );
}