#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;


class Preprocessor{

    private:
        std::vector<Point> points;
        float percentage;
        std::vector<Point> subPoints;
        int optimalLocal_L = 1;
        int optimalSimLocal_L = 6000;
        int optimalSimGlobal_L = 5000;
        void generateSubPoints();


    public:
        Preprocessor(const std::vector<Point>& , float);

        void preProcessInput();

        int getLocal_L();

        int getSimLocal_L();

        int getSimGlobal_L();

        
};