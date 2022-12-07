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
        int transitionMode; // 1 = local, 2 = global
        int subDAlgo; // 1 = incremental, 2 = convex hull
        bool segmentImmunity; // segment marked
        std::string initialization; // 1a, 1b, 2a, 2b
        int m;

        std::vector<Segment_2> untouchableVector;

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

        bool isGlobalValidPath(const std::vector<Segment_2>&, const std::vector<Point>&, const Segment_2&, const std::pair<Point, Point>&);

        FuzzyBox getRectangeBox(std::vector<Point>&);

        ft metropolis(const ft&, const ft&);

        static bool lexOrderPoints(const Point&, const Point&);

        void mergePolygons(std::vector<subTeam>&, std::vector<std::vector<Segment_2>>&);

        void subGlobalTransitions(std::vector<subTeam>&, std::vector<std::vector<Segment_2>>&);

        void subPolygonization(std::vector<subTeam>&, std::vector<std::vector<Segment_2>>&, int);
    
        void createSubsetPoints(std::vector<subTeam>&);

        std::pair<Segment_2, Segment_2> getMarkedSegments(const std::vector<Point>&);

        std::vector<Segment_2> getLowerHull(const std::vector<Point> &);

        bool spatialCondition(const std::vector<subTeam>&, int, int);

    public:

        simulatedAnnealing(const std::vector<Point>&, const std::vector<Segment_2>&, const ft&, const ft&, int, int, int);
        simulatedAnnealing(const std::vector<Point>&, int, int, int, const std::string&, int,int);
        simulatedAnnealing(const std::vector<Point>&, const std::vector<Segment_2>&, const ft&, const ft&, int, int, int, const std::vector<Segment_2>); 
        void startAnnealing();
        void startSubdivision();
        ft getArea();
        ft getRatio(ft CHArea);
};