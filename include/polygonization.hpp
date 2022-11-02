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



class Polygonization{

    protected:
        std::vector<Point> points;
        std::vector<Point> remainingPoints;
        std::vector<Segment_2> polygLine;
        std::vector<Point> polygLinePoints;
        ft totalArea;
        ft ratio;

        
        std::vector<Point> getPolyLinePoints(std::vector<Segment_2>&);
    
        void deleteSegment(std::vector<Segment_2>&, Segment_2&);  
        
        bool forceInsertPoint(std::vector<Segment_2> &, Point &);

        void expandPolygonLine(std::vector<Segment_2>&, Segment_2&, Point&);   

        ft calcRatio(std::vector<Segment_2>&, ft&);
    
    public:
        Polygonization(std::vector<Point>&);
        
        const ft& getArea();
        const ft& getRatio();
        const std::vector<Segment_2>& getPolygonLine();

};


