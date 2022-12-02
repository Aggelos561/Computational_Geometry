#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <utility>
#include <vector>

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Splitters.h>
#include <CGAL/Fuzzy_iso_box.h>

typedef CGAL::Epick Cartesian;
typedef CGAL::Search_traits_2<Cartesian> TreeTraits;
typedef CGAL::Kd_tree<TreeTraits> Kd_tree;
typedef Kd_tree::Tree Tree;
typedef CGAL::Fuzzy_iso_box<TreeTraits> FuzzyBox;

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

typedef struct changes {
  std::pair<Point, Point> pairPointsSeq;
  std::vector<Point> path;
  ft areaDiff;
  Segment_2 segToRemove;
} Changes;

typedef struct transitionStep {
  std::vector<Segment_2> newPolygLine;
  ft areaDiff;
  bool simple;
} transitionStep;

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
        ft areaDiff;
        int L;
        double threshold;
        ft energy;
        
        std::vector<Point> getPolyLinePoints(const std::vector<Segment_2>&);
    
        void deleteSegment(std::vector<Segment_2>&, const Segment_2&);  
        
        bool forceInsertPoint(std::vector<Segment_2> &, const Point &);

        void expandPolygonLine(std::vector<Segment_2>&, const Segment_2&, const Point&);

        ft calcArea(const std::vector<Segment_2> &);

        ft calcPointsArea(std::vector<Point>);

        ft calcRatio(const std::vector<Segment_2>&, const ft&);

        std::vector<Segment_2> getConvexHull(const std::vector<Point>&);  

        std::vector<Point> getPathK(std::vector<Segment_2>&, int, int, std::pair<Point, Point>&);

        ft calculateDeletedArea(std::vector<Point>, const Segment_2&);

        ft calculateAddedArea(std::vector<Point>, const std::pair<Point, Point>&);

        void findChanges(std::vector<Changes>&, std::vector<Point>&, const Segment_2&, const std::pair<Point, Point>&);

        void applyChanges(std::vector<Segment_2>&, std::vector<Changes>&);

        static bool sortAreaChanges(Changes&, Changes&);

        static bool pointsYAscending(const Point&, const Point&);

        bool applyBlueRemoval(std::vector<Segment_2>&, Changes&);

        bool applyKPathRemoval(std::vector<Segment_2>&, Changes&);

        bool checkPolygonSimplicity(std::vector<Segment_2>&);

        bool isValidPath(const std::vector<Segment_2>&, const std::vector<Point>&, const Segment_2&, const std::pair<Point, Point>&);

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
        Polygonization(const std::vector<Point>&, int);
        
        const ft& getArea();
        const ft& getRatio();
        const std::vector<Segment_2>& getPolygonLine();
        void localSearch(std::vector<Segment_2>&);
        void simulatedAnnealing(std::vector<Segment_2>&);
        void spatialSubdivision(std::vector<Point>&, int, const std::string&);
};


