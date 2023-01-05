#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <string>
#include <vector>
#include "preprocessor.hpp"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;

// Namespace dataio is used to read the file input and write result files

namespace showCasedAlgos{

    void  runAlgorithm(const std::string&, const std::vector<Point>&, std::vector<std::pair<double, double>>&, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, const std::chrono::milliseconds, int&, const std::pair<int, std::string>&, Preprocessor&);

    void partialWrite(const std::vector<std::pair<int, std::string>>&, std::vector<std::pair<double, double>>& , std::vector<std::vector<double>>& , std::vector<std::vector<double>>& , const std::pair<int, std::string>& ,  bool& , int , const std::string& , int , bool );

    void initPreprocess(std::vector<Point>, Preprocessor&, const std::vector<std::pair<int, std::string>>&);

}