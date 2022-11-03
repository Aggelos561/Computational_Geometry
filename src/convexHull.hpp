#pragma once

#include "polygonization.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;

class convexHull : public Polygonization{

    private:      
        std::vector<Segment_2> getConvexHull(std::vector<Point> &,std::vector<Point> &);
        
        void initialRun(std::vector<Segment_2>& ,std::vector<Point>& ,std::vector<Segment_2>& );
       
        void initializeConvexHull(std::vector<Segment_2>&, std::vector<Point>&, std::vector<Point>&);
        
        void insertBestPoint(std::vector<pair>&, std::vector<Point>&,std::vector<Segment_2>&);
        
        void findVisiblePoints(std::vector<visPoint>&, Point&, Segment_2&, std::vector<Segment_2>& );

        Point findBestPoint(std::vector<visPoint>& ,std::vector<Point>&, std::vector<Segment_2>&, Segment_2&);
        

    
    public:
        convexHull(std::vector<Point>&,int );
        void start();
};