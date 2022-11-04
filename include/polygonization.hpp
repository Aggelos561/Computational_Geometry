#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;


typedef struct visible_point {
  Point cor;
  double distance;
} visPoint;

typedef struct pairPointSeg {
  Point cor;
  Segment_2 seg;
} pair;


// Base class polygonization.
// Used as base class for Incremental and convexHull classes

class Polygonization{

    protected:
        int edgeSelection;
        std::vector<Point> points;
        std::vector<Point> remainingPoints;
        std::vector<Segment_2> polygLine;
        std::vector<Point> polygLinePoints;
        ft totalArea;
        ft ratio;

        
        std::vector<Point> getPolyLinePoints(const std::vector<Segment_2>&);
    
        void deleteSegment(std::vector<Segment_2>&, const Segment_2&);  
        
        bool forceInsertPoint(std::vector<Segment_2> &, const Point &);

        void expandPolygonLine(std::vector<Segment_2>&, const Segment_2&, const Point&);   

        ft calcRatio(const std::vector<Segment_2>&, const ft&);
    
    public:
        Polygonization(const std::vector<Point>&, int);
        
        const ft& getArea();
        const ft& getRatio();
        const std::vector<Segment_2>& getPolygonLine();

};


