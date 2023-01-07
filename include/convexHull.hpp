#pragma once

#include "polygonization.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;

// Class for convexHull algorithm
// Inherits from polygonization base class

class convexHull : public Polygonization{

    private:

        bool spatial_Subdivision;
        Segment_2 leftConnection;
        Segment_2 rightConnection;
        int polygonIndex;
        int maxPolygonIndex;
        std::vector<Segment_2> initialConvexHull;

        std::vector<Segment_2> getConvexHull(const std::vector<Point> &,std::vector<Point> &);
        
        void initialRun(const std::vector<Segment_2>& ,std::vector<Point>& ,std::vector<Segment_2>& );
       
        void initializeConvexHull(std::vector<Segment_2>&, const std::vector<Point>&, std::vector<Point>&);
        
        void insertBestPoint(const std::vector<pair>&, std::vector<Point>&,std::vector<Segment_2>&);
        
        bool findVisiblePoints(std::vector<visPoint>&, const Point&, const Segment_2&, const std::vector<Segment_2>& );

        Point findBestPoint(const std::vector<visPoint>&);

        static bool sortVisPointsAsc(const visPoint& a , const visPoint& b);
        
    public:
        convexHull(const std::vector<Point>&, int );
        convexHull(const std::vector<Point>&, int ,const Segment_2& , const Segment_2&, int, int);
        void start(const std::chrono::_V2::system_clock::time_point, const std::chrono::milliseconds);
};