#include "../include/dataio.hpp"
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
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


ft Polygonization::calcArea(const std::vector<Segment_2>& polygLine) {
  polygLinePoints = this->getPolyLinePoints(polygLine);
  ft totalArea = CGAL::polygon_area_2(polygLinePoints.begin(), this->polygLinePoints.end(), Convex_hull_traits_2(CGAL::make_property_map(polygLinePoints)));

  return totalArea;
}

ft Polygonization::calcPointsArea(std::vector<Point>& polygPoints) {
  ft totalArea = abs(CGAL::polygon_area_2(polygPoints.begin(), polygPoints.end(), Convex_hull_traits_2(CGAL::make_property_map(polygPoints))));

  return totalArea;
}


// Calculating area / convex hull area ratio for results
ft Polygonization::calcRatio(const std::vector<Segment_2> &convexHull, const ft &area){

  std::vector<Point> convexHullPoints = getPolyLinePoints(convexHull);

  // Calculating convex hull area and divide with total polygon area
  ft convexHullArea = CGAL::polygon_area_2(convexHullPoints.begin(), convexHullPoints.end(), Convex_hull_traits_2(CGAL::make_property_map(convexHullPoints)));

  return area/convexHullArea;
}

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
  
  kPairSequence.second = polygLine[segIndex].target();
  

 return k_path;

}

bool Polygonization::isValidPath(const std::vector<Point>& kPoints, const Segment_2& e){

  Point source = e.source();
  Point target = e.target();

  for (const Point& point : kPoints){
    if ((point == source) || (point == target))
      return false;
  }

  return true;
}


ft Polygonization::calculateDeletedArea(std::vector<Point> points, const Segment_2& e){

  points.push_back(e.source());
  points.push_back(e.target());

  return calcPointsArea(points);

}


ft Polygonization::calculateAddedArea(std::vector<Point> points, const std::pair<Point, Point>& pairSequencePoints){

  points.push_back(pairSequencePoints.second);
  points.push_back(pairSequencePoints.first);
  
  return calcPointsArea(points);

}


void Polygonization::findChanges(std::vector<Changes>& possibleChanges, std::vector<Point>& points , const Segment_2& e, const std::pair<Point, Point>& pairSequencePoints){

  ft removedArea = calculateDeletedArea(points, e);
  
  ft addedArea = calculateAddedArea(points, pairSequencePoints);

  if (addedArea > removedArea){
    Changes change = Changes();

    change.pairPointsSeq = pairSequencePoints;
    change.path = points;
    change.segToRemove = e;
    change.areaDiff = addedArea - removedArea;

    possibleChanges.push_back(change);
  }

}

bool Polygonization::sortAreaChanges(Changes& a, Changes& b){
  
  if (a.areaDiff > b.areaDiff)
    return true;

  return false;
  
}

void Polygonization::applyBlueRemoval(std::vector<Segment_2>& polygLine, Changes& change){

  for (int i = 0; i < polygLine.size(); i++){

    if (polygLine[i] == change.segToRemove){

      Segment_2 blueSeg = polygLine[i];
      polygLine.erase(polygLine.begin() + i);

      if (change.path.size() == 1) {
        if (i + 1 > polygLine.size() - 1){
          polygLine.insert(polygLine.begin()  + 1, Segment_2(polygLine[i].source(), change.path[0]));
        }
        else{
        polygLine.insert(polygLine.begin() + i + 1, Segment_2(polygLine[i].source(), change.path[0]));
        }

        if (i + 2 > polygLine.size() - 1){
         polygLine.insert(polygLine.begin() + i + 1, Segment_2(change.path[0], polygLine[i].target()));
        }
        
        else{
          polygLine.insert(polygLine.begin() + i + 2, Segment_2(change.path[0], polygLine[i].target()));
        }
        
        polygLine.erase(polygLine.begin() + i);
      }
      else {
        
        for (int j = 0; j < change.path.size(); j++){
          
          if (j == 0){
            polygLine.insert(polygLine.begin() + i + j + 1, Segment_2(polygLine[i].source(), change.path[0]));
          }
          else if (j == change.path.size() - 1){
            polygLine.insert(polygLine.begin() + i + j + 1, Segment_2(change.path[change.path.size() - 1], polygLine[i].target()));
          }
          else{
            polygLine.insert(polygLine.begin() + i + j + 1, Segment_2(change.path[j - 1], change.path[j]));
          }
        
        }

        polygLine.erase(polygLine.begin() + i);
      }
      break;
    }
  }
}


void Polygonization::applyKPathRemoval(std::vector<Segment_2>& polygLine, Changes& change){

  int insertionIndex = 0;
  for (insertionIndex = 0; insertionIndex < polygLine.size(); insertionIndex++){
    
    if (polygLine[insertionIndex].source() == change.pairPointsSeq.first){
      polygLine.insert(polygLine.begin() + insertionIndex + 1, Segment_2(change.pairPointsSeq.first, change.pairPointsSeq.second));
      break;
    }
  
  }

  for (int i = insertionIndex; i < polygLine.size(); i++){

    if (polygLine[i].source() == change.path[0]){

      for (int j = 0; j < change.path.size(); j++){
        polygLine.erase(polygLine.begin() + i);
      }

      break;
    }

    if (i == polygLine.size() - 1){
      i = -1;
    }
  }

}


bool Polygonization::checkPolygonSimplicity(std::vector<Segment_2>& polygLine){

  Polygon_2 polygon = Polygon_2();

  for (Segment_2& segment : polygLine){
    polygon.push_back(segment.source());
  }

  return polygon.is_simple();

}

void Polygonization::applyChanges(std::vector<Segment_2>& polygLine, std::vector<Changes>& possibleChanges){

  std::sort(possibleChanges.begin(), possibleChanges.end(), this->sortAreaChanges);

  std::vector<Segment_2> changedPolygLine = polygLine, prevPolygon = polygLine;

  for (Changes& change : possibleChanges){
    
    std::cout << "Test_0" << std::endl;
    applyBlueRemoval(changedPolygLine, change);
    std::cout << "Test_1" << std::endl;

    applyKPathRemoval(changedPolygLine, change);
    std::cout << "Test_2" << std::endl;
    
    for (Segment_2& segment : changedPolygLine){
      std::cout << "Changed --> " << segment << std::endl;
    }

    if (!checkPolygonSimplicity(changedPolygLine)){
      changedPolygLine = prevPolygon;
    }
    else{

      for (Segment_2& segment : changedPolygLine){
        std::cout << "SIMPLE --> " << segment << std::endl;
      }

      prevPolygon = changedPolygLine;
    }
  }

  for (Segment_2& segment : polygLine){
    std::cout << "Previous --> " << segment << std::endl;
  }
  std::cout << "-------------------------------------" << std::endl;

  for (Segment_2& segment : changedPolygLine){
    std::cout << "New --> " << segment << std::endl;
  }
  std::cout << "=====================================" << std::endl;

}


void Polygonization::localSearch(std::vector<Segment_2>& polygLine){

  std::vector<Changes> possibleChanges;
  ft Da = -1;
  threshold = 1;
  this->L = 1;

  do{

    for (Segment_2& e : polygLine) {
      
      int segmentIndex = 0;
      std::cout << " Next segment is " << e << std::endl;
        
      for (Segment_2& v : polygLine){

        for(int k = 1; k <= this->L; k++){

          std::pair<Point, Point> kPairSequence;
          std::vector<Point> kPoints = getPathK(polygLine, segmentIndex, k, kPairSequence);  

          if (!isValidPath(kPoints, e)){
              break;
          }

          findChanges(possibleChanges, kPoints, e, kPairSequence);

          for (const Point& p : kPoints){
            std::cout << "Point --> " << p << std::endl;
          }

        }
        segmentIndex++;
      }

      applyChanges(polygLine, possibleChanges);
    }
  } while(Da >= threshold);

  // for (const Changes& change: possibleChanges){
  //   std::cout << "Blue --> " << change.segToRemove << " Pair --> (" << change.pairPointsSeq.first <<", " << change.pairPointsSeq.second << ")" << std::endl ;
  // }

}


//while ∆A ≥ threshold do
//for every edge e ∈ S do
// if V moving to e increases area and retains simplicity then
// list T ← [e, V ]
// end if
// end for
// end for
// Apply all changes in T to S
// Keep best solution S′; ∆A ← Area(S′) − Area(S)
// end while