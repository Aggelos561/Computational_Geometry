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

    std::string getRunType(int, char**);

    bool getExecutionParameters(std::string&, std::string& , std::string& , std::string& ,std::string&, int&, int&, double&, std::string&, int&, int, char** , int&);

    bool getScoreParameters(std::string&, std::string&, bool&, int, char**);

    std::vector<Point> readPoints(const std::string& );

    std::vector<std::pair<int, std::string>> findFiles(const std::string& );

    void createResultsFile(const std::vector<Segment_2> &, const ft&, const ft&, const ft&, const ft&, const std::chrono::milliseconds&, const std::string&, const std::string&, const int&, const std::string&, const std::string&);

    void writeToOutputFile(const std::string&, const std::vector<std::pair<double, double>>&, const std::vector<double>&, const std::vector<double>&, int, bool);
}