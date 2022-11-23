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

ft Polygonization::calcPointsArea(std::vector<Point> polygPoints) {

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
  
  kPairSequence.second = polygLine[segIndex < 0 ? 0 : segIndex].target();
  

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

  

  // std::cout << "Deleted Area = " << (long int) removedArea << " , Added Area = " << (long int) addedArea << ", Area Diff = " << (long int) addedArea - removedArea << std::endl;
  Changes change = Changes();

  change.pairPointsSeq = pairSequencePoints;
  change.path = points;
  change.segToRemove = e;
  change.areaDiff = addedArea - removedArea;

  possibleChanges.push_back(change);


}

bool Polygonization::sortAreaChanges(Changes& a, Changes& b){
  
  if (a.areaDiff > b.areaDiff)
    return true;

  return false;
  
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
            polygLine.insert(polygLine.begin() + (i + j + 1) % polygLine.size(), Segment_2(change.path[change.path.size() - 1], polygLine[i].target()));
            
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

  if (!segmentFound){
    return false;
  }
  return true;

}


void Polygonization::applyKPathRemoval(std::vector<Segment_2>& polygLine, Changes& change){

  // int count = 1;
  // for (Segment_2& segment : polygLine){
  //   std::cout << count << ". KPATH0 ~> " << segment << std::endl;
  //   count++;
  // }
  
  int startIndex = 0;
  for (int i = 0; i < polygLine.size(); i++){
    if (polygLine[i].target() == change.path[0]){
      startIndex = i;
    }
  }

  int k = change.path.size() + 1;

  // std::cout << "a"<< std::endl;

  // std::cout << "index => " << startIndex << " max index => " << polygLine.size() - 1 << std::endl; 
  while(k != 0){
    
    polygLine.erase(polygLine.begin() + startIndex);
    k--;

    if (startIndex >= polygLine.size() - 1)
      startIndex = 0;
  }
  

  // count = 1;
  // std::cout << std::endl;
  // for (Segment_2& segment : polygLine){
  //   std::cout << count << ". KPATH1 ~> " << segment << std::endl;
  //   count++;
  // }

  int insertionIndex = 0;
  for (insertionIndex = 0; insertionIndex < polygLine.size(); insertionIndex++){

    if (polygLine[insertionIndex].target() == change.pairPointsSeq.first){

      polygLine.insert(polygLine.begin() + (insertionIndex + 1) % polygLine.size(), Segment_2(change.pairPointsSeq.first, change.pairPointsSeq.second));      
      break;
      
    }
  
  }

  // count = 1;
  // std::cout << std::endl;
  // for (Segment_2& segment : polygLine){
  //   std::cout << count << ". KPATH2 ~> " << segment << std::endl;
  //   count++;
  // }
}


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

void Polygonization::applyChanges(std::vector<Segment_2>& polygLine, std::vector<Changes>& possibleChanges){

  std::sort(possibleChanges.begin(), possibleChanges.end(), this->sortAreaChanges);

  std::vector<Segment_2> changedPolygLine = polygLine, prevPolygon = polygLine;

  for (Changes& change : possibleChanges){

    applyKPathRemoval(changedPolygLine, change);

    bool blueRemoved = applyBlueRemoval(changedPolygLine, change);

    if (!checkPolygonSimplicity(changedPolygLine) || !blueRemoved){
      changedPolygLine = prevPolygon;
    }
    else{
      
      if (calcArea(changedPolygLine) < calcArea(polygLine)){

        prevPolygon = changedPolygLine;
        polygLine = prevPolygon;

        ft newArea = calcArea(polygLine);


        std::cout << "New Polygon Area --> " << (long int) newArea << std::endl;
        std::cout << "Change Kpoint --> " << change.path[0] << std::endl;
        std::cout << "Blue Segment --> " << change.segToRemove << std::endl;

        for (Segment_2& segment : changedPolygLine){
          std::cout << "SIMPLE --> " << segment << std::endl;
        }

        std::cout << std::endl;

        break;
      }
    }
  }


}


void Polygonization::localSearch(std::vector<Segment_2>& polygLine){

  ft Da = 1000;
  threshold = 1;
  this->L = 3;

  // polygLine = {Segment_2(Point(40, 20), Point(30, 30)), Segment_2(Point(30, 30), Point(22, 24)), Segment_2(Point(22, 24), Point(19, 22)),
  //               Segment_2(Point(19, 22), Point(15, 26)), Segment_2(Point(15, 26), Point(0, 50)), Segment_2(Point(0, 50), Point(-2, 25)),
  //               Segment_2(Point(-2, 25), Point(16, 4)), Segment_2(Point(16, 4), Point(25, 3)), Segment_2(Point(25, 3), Point(33, 13)),
  //               Segment_2(Point(33, 13), Point(40, 20))};

  std::vector<Segment_2> test = polygLine;

  do{

    std::vector<Changes> possibleChanges;
    
    for (Segment_2& e : polygLine) {
 
      int segmentIndex = 0;
      // std::cout << " Next segment is " << e << std::endl;

      for (Segment_2& v : polygLine){

        for(int k = 1; k <= this->L; k++){

          std::pair<Point, Point> kPairSequence;
          std::vector<Point> kPoints = getPathK(polygLine, segmentIndex, k, kPairSequence);

          if (!isValidPath(kPoints, e))
            break;

          findChanges(possibleChanges, kPoints, e, kPairSequence);

        }
        segmentIndex++;
      }

    }

    applyChanges(polygLine, possibleChanges);

    threshold++;

  } while(Da >= threshold);


  for (Segment_2& e : test) {

    std::cout << e << std::endl;

  }

  std::cout << std::endl; 
  std::cout << "Area Before ==> " << (long int) calcArea(test) << std::endl;
  std::cout << "Area After ==> " << (long int) calcArea(polygLine) << std::endl;

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