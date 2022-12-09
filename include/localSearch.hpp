
#pragma once

#include "polygonization.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;


class localSearch: public Polygonization{

    private:

        ft optimisedArea;
        ft optimisedRatio;
        ft areaDiff;
        int L;
        double threshold;
        int mode; // 1 = min, 2 = max

        void findChanges(std::vector<Changes>&, std::vector<Point>&, const Segment_2&, const std::pair<Point, Point>&);

        void applyChanges(std::vector<Segment_2>&, std::vector<Changes>&);

        static bool sortAreaChanges(Changes&, Changes&);

    public:
        localSearch(const std::vector<Point>&, const std::vector<Segment_2>&, const ft&, const ft&, int, double, int);
        void start();
        ft getOptimisedArea();
        ft getOptimisedRatio();
};