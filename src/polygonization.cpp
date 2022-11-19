#include "../include/dataio.hpp"
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <iostream>
#include <vector>
#include "../include/polygonization.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point>::type> Convex_hull_traits_2;
typedef CGAL::Epick::FT ft;


// Contructor for base polygonization class
Polygonization::Polygonization(const std::vector<Point>& points, int edgeSelection){
	this->points = points;
	this->remainingPoints = points;
  this->edgeSelection = edgeSelection;
  this->totalArea = 0;
  this->ratio = 0;
}

// Return total area of polygon
const ft& Polygonization::getArea(){
	return totalArea;
}

// Return ratio for results file
const ft& Polygonization::getRatio(){
	return ratio;
}


// Return polygon line vector
const std::vector<Segment_2>& Polygonization::getPolygonLine(){
	return polygLine;
}



// Get all points from a vector of segments
std::vector<Point> Polygonization::getPolyLinePoints(const std::vector<Segment_2> &polygLine) {

  std::vector<Point> linePoints;

  for (int i = 0; i < polygLine.size(); i++) {
    linePoints.push_back(polygLine[i].point(0));
  }

  return linePoints;
}



// Deletes a spesific segment value from a vector of segments
void Polygonization::deleteSegment(std::vector<Segment_2> &polygLine, const Segment_2 &visibleSegment){
  polygLine.erase(std::remove(polygLine.begin(), polygLine.end(), visibleSegment), polygLine.end());
}



// if new point already on polygon line the force insert into polygon line because no visible segment exists
bool Polygonization::forceInsertPoint(std::vector<Segment_2> &polygLine, const Point &nextPoint) {

  bool forceInserted = false;

  for (int i = 0; i < polygLine.size(); i++) {
    if (polygLine[i].has_on(nextPoint)) {
      polygLine.insert(polygLine.begin() + i + 1, Segment_2(polygLine[i].source(), nextPoint));
      polygLine.insert(polygLine.begin() + i + 2, Segment_2(nextPoint, polygLine[i].target()));
      polygLine.erase(polygLine.begin() + i);
      forceInserted = true;
      break;
    }
  }

  // Check if point was force inserted
  if (forceInserted) {
    return true;
  }
	return false;
}



// Insert the two new segments in the right place in polygon line for the new point
void Polygonization::expandPolygonLine(std::vector<Segment_2> &polygLine, const Segment_2 &visibleSegment, const Point &nextPoint){

  int index = 0;
  for (int i = 0; i < polygLine.size(); i++) {
    if (polygLine[i].point(1) == visibleSegment.point(0)){
      
      polygLine.insert(polygLine.begin() + index + 1,Segment_2(visibleSegment.point(0), nextPoint));
      polygLine.insert(polygLine.begin() + index + 2, Segment_2(nextPoint, visibleSegment.point(1)));
      break;
      
    } 
    else if (polygLine[i].point(1) == visibleSegment.point(1)) {
      polygLine.insert(polygLine.begin() + index + 1,Segment_2(visibleSegment.point(1), nextPoint));
      polygLine.insert(polygLine.begin() + index + 2, Segment_2(nextPoint, visibleSegment.point(0)));
      break;
    }
    index++;
  }

}



// Calculating area / convex hull area ratio for results
ft Polygonization::calcRatio(const std::vector<Segment_2> &convexHull, const ft &area){

  std::vector<Point> convexHullPoints = getPolyLinePoints(convexHull);

  // Calculating convex hull area and divide with total polygon area
  ft convexHullArea = CGAL::polygon_area_2(convexHullPoints.begin(), convexHullPoints.end(), Convex_hull_traits_2(CGAL::make_property_map(convexHullPoints)));

  return area/convexHullArea;
}

std::vector<Point> Polygonization::getPath(std::vector<Segment_2>& polygLine ,std::vector<Point>& pointsOfPolygon,int k){
  std::vector<Point> path;
  for(int m = 0,a = 0; m < pointsOfPolygon.size(), a < k; m++,a++){
    path.push_back(pointsOfPolygon[m]);
  }
}

void Polygonization::optimizeLocalSearch(std::vector<Segment_2>& polygLine){
  std::vector<Changes> vecOfChanges;
  ft Da;
  do{
    //while ∆A ≥ threshold do
    //for every edge e ∈ S do
    for (int i = 0; i < polygLine.size(); i++) {
        for(int k = 1; k <= this->L; k++){
          for(int j = 0; j < k ; j++){
            std::vector<Point> pointsOfPolygon = getPolyLinePoints(polygLine);
            // for(int m = 0,a = 0; m < pointsOfPolygon.size(), a < k; m++,a++)
            //   path.push_back(pointsOfPolygon[m]);
          }
            if V moving to e increases area and retains simplicity then
              list T ← [e, V ]
            end if
          end for
        end for
        Apply all changes in T to S
        Keep best solution S′; ∆A ← Area(S′) − Area(S)
        end while
      }
    }
  }while(Da >= threshold)
}