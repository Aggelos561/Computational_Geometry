#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;


class Preprocessor{

    private:
        std::vector<Point> points;
        std::map<int, int> sizeToLocalL;
        std::map<int, int> sizeToGlobalL;
        std::map<int, int> sizeToSubdivisionL;
        int optimal_M;
        const float percentage = 0.10;
        
        int getSizeBucket(int );

    public:

        void defaultInput();

        void preprocessInput(const std::vector<Point>& );

        int getSimSubDiv_M(int );

        int getSimLocal_L(int );

        int getSimGlobal_L(int );

        int getSimSubDiv_L(int );

        std::vector<Point> generateSubPoints();

};