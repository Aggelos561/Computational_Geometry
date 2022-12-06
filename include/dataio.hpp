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

    bool getParameters(std::string& nameOfFile, std::string& outputFile, std::string& algorithm, std::string& algorithm_initial,std::string& initial, int& edge_selection_int,double& threshold,int& annealing,int& L, int argc, char** argv);

    std::vector<Point> readPoints(const std::string& );

    void createResultsFile(const std::vector<Segment_2> &polygLine, const ft& area, const ft& areaBefore, const ft& ratio,const ft& ratioBefore, const std::chrono::milliseconds& polygonizationDuration, const std::string& output, const std::string& algorithm, const int& edgeSelection, const std::string& initialization);
}