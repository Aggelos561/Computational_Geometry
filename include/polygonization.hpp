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

    private:
        std::vector<Point> points;
        std::vector<Point> remainingPoints;
        std::vector<Segment_2> polygLine;
        std::vector<Point> polygLinePoints;
        ft totalArea;
        ft ratio;
        std::chrono::milliseconds duration;

        static bool sortYAsc(Point &, Point &);
        static bool sortYDesc(Point &, Point &);
        void sortPoints(std::vector<Point>&);

        std::vector<Point> getPolyLinePoints(std::vector<Segment_2>&);
        
        std::vector<Segment_2> getConvexHull(std::vector<Point>&);
        std::vector<Segment_2> getConvexHull(std::vector<Point> &,std::vector<Point> &);
        
        void initialRun(std::vector<Segment_2>& ,std::vector<Point>& ,std::vector<Segment_2>& );
        
        void initializeTriangle(std::vector<Segment_2>&, std::vector<Point>&, std::vector<Point>&);
        void initializeConvexHull(std::vector<Segment_2>&, std::vector<Point>&, std::vector<Point>&);
        void deleteSegment(std::vector<Segment_2>&, Segment_2&);
        
        ft getTriangleArea(Segment_2&, Point&);
        
        Segment_2 chooseVisibleSegment(std::vector<Segment_2>&, Point&, ft&);
        bool forceInsertPoint(std::vector<Segment_2> &, Point &);

        void expandPolygonLine(std::vector<Segment_2>&, Segment_2&, Point&);
        void insertBestPoint(std::vector<pair>&, std::vector<Point>&,std::vector<Segment_2>& , int );

        std::vector<Segment_2> getRedSegments(std::vector<Segment_2>&, std::vector<Segment_2>&, Point&);
        std::vector<Segment_2> findVisibleSegments(std::vector<Segment_2>&, std::vector<Segment_2>&, Point&);
        void findVisiblePoints(std::vector<visPoint>&, Point&, Segment_2&, std::vector<Segment_2>& );

        Point findBestPoint(std::vector<visPoint>& ,std::vector<Point>&, std::vector<Segment_2>&, Segment_2&);
        

        ft calcRatio(std::vector<Segment_2>&, ft&);
    
    public:
        Polygonization(std::vector<Point>&);
        
        const ft& getArea();
        const ft& getRatio();
        const std::chrono::milliseconds& getDuration();
        const std::vector<Segment_2>& getPolygonLine();

        void incremental();
        void convexHullAlg();
};


