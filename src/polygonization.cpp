#include "../include/dataio.hpp"
#include <CGAL/Bbox_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/enum.h>
#include <CGAL/intersections.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>
#include "../include/polygonization.hpp"

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Splitters.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Kd_tree_rectangle.h>

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


// Contructor for base polygonization class
Polygonization::Polygonization(const std::vector<Point>& points, int edgeSelection){
	this->points = points;
	this->remainingPoints = points;
  this->edgeSelection = edgeSelection;
  this->totalArea = 0;
  this->ratio = 0;
}


Polygonization::Polygonization(const std::vector<Point>& points, const std::vector<Segment_2>& polygLine, const ft& area, const ft& ratio){
	this->points = points;
  this->polygLine = polygLine;
  this->edgeSelection = 0;
  this->totalArea = area;
  this->ratio = ratio;
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


ft Polygonization::calcArea(const std::vector<Segment_2>& polygLine) {
  polygLinePoints = this->getPolyLinePoints(polygLine);
  ft totalArea = CGAL::polygon_area_2(polygLinePoints.begin(), this->polygLinePoints.end(), Convex_hull_traits_2(CGAL::make_property_map(polygLinePoints)));

  return totalArea;
}

ft Polygonization::getCHArea(){
  return this->CHArea;
}

// Calculating area / convex hull area ratio for results
ft Polygonization::calcRatio(const std::vector<Segment_2> &convexHull, const ft &area){

  std::vector<Point> convexHullPoints = getPolyLinePoints(convexHull);

  // Calculating convex hull area and divide with total polygon area
  ft convexHullArea = CGAL::polygon_area_2(convexHullPoints.begin(), convexHullPoints.end(), Convex_hull_traits_2(CGAL::make_property_map(convexHullPoints)));
  this->CHArea = convexHullArea;
  return area/convexHullArea;
}



// Find Complex Hull Points and rerurn a vector of segments
std::vector<Segment_2> Polygonization::getConvexHull(const std::vector<Point> &polygLinePoints) {

  std::vector<Point> convexHullPoints;

  CGAL::convex_hull_points_2(polygLinePoints.begin(), polygLinePoints.end(), std::back_inserter(convexHullPoints));

  Polygon_2 convexHullPolygon = Polygon_2();

  for (const Point& p : convexHullPoints) {
    convexHullPolygon.push_back(p);
  }

  std::vector<Segment_2> convexHullSegments;

  // Create A Vector That Stores Convex Hull segments
  for (const Segment_2& seg : convexHullPolygon.edges()) {
    convexHullSegments.push_back(seg);
  }

  return convexHullSegments;
}


// Get K path based on k number
std::vector<Point> Polygonization::getPathK(std::vector<Segment_2>& polygLine, int segmentIndex, int k, std::pair<Point, Point>& kPairSequence){
  
  std::vector<Point> k_path;

  int segIndex;
  for (segIndex = segmentIndex; segIndex < polygLine.size(); segIndex++){
    k_path.push_back(polygLine[segIndex].source());

    if (segIndex == (polygLine.size() - 1)){
      segIndex = -1;
    }

    k--;

    if (k == 0)
      break;
  }

  if (segmentIndex - 1 < 0){
    kPairSequence.first = polygLine[polygLine.size() - 1].source();
  }
  else {
    kPairSequence.first = polygLine[segmentIndex - 1].source();
  }
  
  kPairSequence.second = polygLine[segIndex < 0 ? 0 : segIndex].target();
  

 return k_path;

}



bool Polygonization::isValidPath(const std::vector<Segment_2>& polygLine, const std::vector<Point>& kPoints, const Segment_2& e, const std::pair<Point, Point>& kPairSequence){

  Point source = e.source();
  Point target = e.target();

  for (const Point& point : kPoints){
    if ((point == source) || (point == target))
      return false;
  }

  return true;
}


ft Polygonization::calculateDeletedArea(std::vector<Point> points, const Segment_2& e){

  points.push_back(e.target());
  points.push_back(e.source());
 

  return calcPointsArea(points);

}



ft Polygonization::calculateAddedArea(std::vector<Point> points, const std::pair<Point, Point>& pairSequencePoints){

  points.push_back(pairSequencePoints.second);
  points.push_back(pairSequencePoints.first);

  return calcPointsArea(points);

}




bool Polygonization::applyBlueRemoval(std::vector<Segment_2>& polygLine, Changes& change){

  bool segmentFound = false;
  for (int i = 0; i < polygLine.size(); i++){

    if (polygLine[i] == change.segToRemove){

      segmentFound = true;
      Segment_2 blueSeg = polygLine[i];

      if (change.path.size() == 1) {

        polygLine.insert(polygLine.begin() + (i + 1) % polygLine.size(), Segment_2(polygLine[i].source(), change.path[0]));
        polygLine.insert(polygLine.begin() + (i + 2) % polygLine.size(), Segment_2(change.path[0], polygLine[i].target()));

      }
      else if (change.path.size() == 2){
        
        for (int j = 0; j < change.path.size(); j++){
          
          if (j == 0){
            polygLine.insert(polygLine.begin() + (i + 1) % polygLine.size(), Segment_2(polygLine[i].source(), change.path[0]));
          }
          else if (j == change.path.size() - 1){
            polygLine.insert(polygLine.begin() + (i + j + 1) % polygLine.size(), Segment_2(change.path[j - 1], change.path[j]));
            polygLine.insert(polygLine.begin() + (i + j + 2) % polygLine.size(), Segment_2(change.path[change.path.size() - 1], polygLine[i].target()));
            
          }
        
        }

      }
      else{

        for (int j = 0; j < change.path.size(); j++){
          
          if (j == 0){
            polygLine.insert(polygLine.begin() + (i + 1) % polygLine.size(), Segment_2(polygLine[i].source(), change.path[0]));
          }
          else if (j == change.path.size() - 1){        
            polygLine.insert(polygLine.begin() + (i + j + 1) % polygLine.size(), Segment_2(change.path[j - 1], change.path[j]));
            polygLine.insert(polygLine.begin() + (i + j + 2) % polygLine.size(), Segment_2(change.path[change.path.size() - 1], polygLine[i].target()));
        
          }
          else{
            polygLine.insert(polygLine.begin() + (i + j + 1) % polygLine.size(), Segment_2(change.path[j - 1], change.path[j]));
          }
        
        }

      }

      polygLine.erase(polygLine.begin() + i);
      break;
    }
  }

  return segmentFound;

}


//  Remove K path from above
bool Polygonization::applyKPathRemoval(std::vector<Segment_2>& polygLine, Changes& change){
  
  int startIndex = 0;

  for (int i = 0; i < polygLine.size(); i++){
    if (polygLine[i].target() == change.path[0]){
      startIndex = i;
    }
  }

  int k = change.path.size() + 1;

  while(k != 0){
    
    polygLine.erase(polygLine.begin() + startIndex);
    k--;

    if (startIndex >= polygLine.size() - 1)
      startIndex = 0;
  }
  
  int insertionIndex = 0;
  for (insertionIndex = 0; insertionIndex < polygLine.size(); insertionIndex++){

    if (polygLine[insertionIndex].target() == change.pairPointsSeq.first){

      polygLine.insert(polygLine.begin() + (insertionIndex + 1) % polygLine.size(), Segment_2(change.pairPointsSeq.first, change.pairPointsSeq.second));      
      break;
      
    }
  
  }

  Polygon_2 polygon = Polygon_2();

  for (Segment_2& segment : polygLine){
    polygon.push_back(segment.source());
  }

  for (Point& p : change.path){

    K traits;

    CGAL::Bounded_side result = CGAL::bounded_side_2(polygon.begin(), polygon.end(), p, traits);

    if (result != CGAL::ON_BOUNDED_SIDE){
      return true;
    }

  }

  return false;

}


// Area based on ccw or cw points
ft Polygonization::calcPointsArea(std::vector<Point> polygPoints) {

  ft totalArea = abs(CGAL::polygon_area_2(polygPoints.begin(), polygPoints.end(), Convex_hull_traits_2(CGAL::make_property_map(polygPoints))));

  return totalArea;
}


//  debugging method for checking simplicity
bool Polygonization::checkPolygonSimplicity(std::vector<Segment_2>& polygLine){

  Polygon_2 polygon = Polygon_2();

  for (Segment_2& segment : polygLine){
    polygon.push_back(segment.source());
  }


  for (int i = 0; i < polygLine.size(); i++){

    if (polygLine[i % polygLine.size()].target() != polygLine[(i+1) % polygLine.size()].source()){
      return false;
    }
  }
  return polygon.is_simple();

}

// compare segments having same source target or reverse
bool Polygonization::looseSegCompare(const Segment_2& seg1, const Segment_2& seg2){

  if (seg1.source() == seg2.source() && seg1.target() == seg2.target())
    return true;
  
  if (seg1.source() == seg2.target() && seg1.target() == seg2.source())
    return true;
  
  return false;

}

// For local search and global transition only. Checking spesific segments intersection
bool Polygonization::isSimple(const Changes& change, const std::vector<Segment_2>& newPolygonLine){
  
  Segment_2 addedSegmentleft(change.segToRemove.source(), change.path[0]);
  Segment_2 addedSegmentRight(change.path[change.path.size() - 1], change.segToRemove.target());
  Segment_2 upperAddedSegment(change.pairPointsSeq.first, change.pairPointsSeq.second);

  auto result1 = CGAL::intersection(upperAddedSegment, addedSegmentleft);
  auto result2 = CGAL::intersection(upperAddedSegment, addedSegmentRight);

  if (result1 || result2){
    return false;
  }

  std::vector<Segment_2> possibleInterSegs = {addedSegmentleft, addedSegmentRight, upperAddedSegment};


  for (const Segment_2& newSegment : possibleInterSegs){
    for (const Segment_2& newPolygonSegment : newPolygonLine){

      if (newPolygonSegment == newSegment)
        continue;

      const auto result = CGAL::intersection(newSegment, newPolygonSegment);

      if (result){
        if (const Point *p = boost::get<Point>(&*result)) {

          if (*p == newPolygonSegment.source() || *p == newPolygonSegment.target()) {
            continue;
          } 
          else {
            return false;
          }
        } 
        else {
        return false;
        }
      }
    }
  }

  return true;
}