#pragma once

#include "polygonization.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Splitters.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point>::type> Convex_hull_traits_2;
typedef CGAL::Epick::FT ft;
typedef CGAL::Epick Cartesian;
typedef CGAL::Search_traits_2<Cartesian> TreeTraits;
typedef CGAL::Kd_tree<TreeTraits> Kd_tree;
typedef Kd_tree::Tree Tree;
typedef CGAL::Fuzzy_iso_box<TreeTraits> FuzzyBox;



class simulatedAnnealing: public Polygonization{

    private:

        ft optimisedArea;
        ft optimisedRatio;
        int L;
        ft energy;
        int mode; // 1 = min, 2 = max
        int subDAlgo;

        static bool pointsYAscending(const Point&, const Point&);

        ft calcAreaDiff(const std::vector<Segment_2>&, const Point&, const Point&, const Point&, const Point&);

        ft energyCalc(const ft&, const ft&);
        
        transitionStep localTransition(std::vector<Segment_2>&, const Tree&);  

        transitionStep globalTransition(std::vector<Segment_2>&);

        void findGlobalChanges(std::vector<Changes>&, std::vector<Point>& , const Segment_2&, const std::pair<Point, Point>&);

        transitionStep applyGlobalChanges(std::vector<Segment_2>&, std::vector<Changes>&);

        void replace(const Segment_2&, const Segment_2&, const Segment_2&, const Segment_2&, const Segment_2&, std::vector<Segment_2>&, int, int, int);

        void KdTreeInit(const std::vector<Segment_2>&, Tree&);

        bool validityCheck(const Tree&, const std::vector<Segment_2>&, const Segment_2&, const Segment_2&, const Segment_2&);

        FuzzyBox getRectangeBox(std::vector<Point>&);

        ft metropolis(const ft&, const ft&);

        static bool lexOrderPoints(const Point&, const Point&);

    public:
        simulatedAnnealing(const std::vector<Point>&, const std::vector<Segment_2>&, const ft&, const ft&, int, int, int);
        void startAnnealing();
        void startSubdivision(std::vector<Point>&, int, const std::string&);
    
};