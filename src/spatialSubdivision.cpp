#include "../include/dataio.hpp"
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/enum.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>
#include "../include/polygonization.hpp"
#include "../include/simulatedAnnealing.hpp"
#include "../include/incremental.hpp"
#include "../include/convexHull.hpp"



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point>::type> Convex_hull_traits_2;
typedef CGAL::Epick::FT ft;



// Spatial Condition
// Needs to be monote
// collinear checking so ch algorithm wont break segments
bool simulatedAnnealing::spatialCondition(const std::vector<subTeam>& subPolPoints, int condIndex, int currentIndex){

  bool conditionSpat = (subPolPoints[condIndex].markedSegments.second.source().y() < subPolPoints[condIndex].markedSegments.second.target().y() && 
  subPolPoints[condIndex].markedSegments.second.source().x() < subPolPoints[condIndex].markedSegments.second.target().x()) && 
  (subPolPoints[currentIndex].markedSegments.first.source().y() > subPolPoints[currentIndex].markedSegments.first.target().y() 
  && subPolPoints[currentIndex].markedSegments.first.source().x() < subPolPoints[currentIndex].markedSegments.first.target().x());

  std::vector<Point> remainingPointsLeft = subPolPoints[condIndex].points;
  std::vector<Point> remainingPointsRight = subPolPoints[currentIndex].points;

  remainingPointsLeft.erase(std::remove(remainingPointsLeft.begin(), remainingPointsLeft.end(), subPolPoints[condIndex].markedSegments.second.source()), remainingPointsLeft.end());
  remainingPointsLeft.erase(std::remove(remainingPointsLeft.begin(), remainingPointsLeft.end(), subPolPoints[condIndex].markedSegments.second.target()), remainingPointsLeft.end());

  remainingPointsRight.erase(std::remove(remainingPointsRight.begin(), remainingPointsRight.end(), subPolPoints[currentIndex].markedSegments.first.source()), remainingPointsRight.end());
  remainingPointsRight.erase(std::remove(remainingPointsRight.begin(), remainingPointsRight.end(), subPolPoints[currentIndex].markedSegments.first.target()), remainingPointsRight.end());


  for (const Point& p : remainingPointsLeft){
    if (CGAL::collinear(subPolPoints[condIndex].markedSegments.second.source(), subPolPoints[condIndex].markedSegments.second.target(), p)){
      return false;
    }
  }

  for (const Point& p : remainingPointsRight){
    if (CGAL::collinear(subPolPoints[condIndex].markedSegments.first.source(), subPolPoints[condIndex].markedSegments.first.target(), p)){
      return false;
    }
  }

  return conditionSpat;
}


// Get lower hull for spesific set of points
std::vector<Segment_2> simulatedAnnealing::getLowerHull(const std::vector<Point> &polygLinePoints) {

  std::vector<Point> lowerHullPoints;

  CGAL::lower_hull_points_2(polygLinePoints.begin(), polygLinePoints.end(), std::back_inserter(lowerHullPoints));

  Polygon_2 lowerHullLine = Polygon_2();

  for (const Point& p : lowerHullPoints) {
    lowerHullLine.push_back(p);
  }

  std::vector<Segment_2> lowerHullSegments;

  // Create A Vector That Stores Convex Hull segments
  for (const Segment_2& seg : lowerHullLine.edges()) {
    lowerHullSegments.push_back(seg);
  }

  return lowerHullSegments;
}




// Merging all subpolygons into one final polygon
void simulatedAnnealing::mergePolygons(std::vector<subTeam>& subPolPoints, std::vector<std::vector<Segment_2>>& allPolygLines){

  for (int i = 0; i < allPolygLines.size(); i++){
    
    int indexStart;
    int indexEnd;

    for (int j = 0; j < allPolygLines[i].size(); j++){
    
      if (allPolygLines[i][j].source() == subPolPoints[i].markedSegments.first.source()){
        indexStart = j;
      }
      if (allPolygLines[i][j].target() == subPolPoints[i].markedSegments.second.target()){
          indexEnd = j;
      }
    }

    
    while (indexStart != indexEnd){
       polygLine.push_back(allPolygLines[i][indexStart]);
       indexStart = indexStart + 1 > allPolygLines[i].size() - 1 ? 0 : indexStart + 1;
    }

    if (i < allPolygLines.size() - 1){
      polygLine.push_back(allPolygLines[i][indexStart]);
    }

  }

  for (int i = allPolygLines.size() - 1; i >= 0; i--){
    
    int indexStart;
    int indexEnd;
    
    for (int j = 0; j < allPolygLines[i].size(); j++){

      if (allPolygLines[i][j].source() == subPolPoints[i].markedSegments.second.target()){
        indexStart = j;
      }
    }

    if (i == allPolygLines.size() - 1){
      polygLine.push_back(allPolygLines[i][indexStart - 1 >= 0 ? indexStart - 1 : allPolygLines[i].size() - 1]);
    }

    for (int j = 0; j < allPolygLines[i].size(); j++){

      if (allPolygLines[i][j].source() == polygLine[polygLine.size() - 1].target()){
        indexStart = j;
      }
      if (allPolygLines[i][j].target() == subPolPoints[i].markedSegments.first.source()){
         indexEnd = j;
      }
    }

    while (indexStart != indexEnd){
       polygLine.push_back(allPolygLines[i][indexStart]);
       indexStart = indexStart + 1 > allPolygLines[i].size() - 1 ? 0 : indexStart + 1;
    }
    polygLine.push_back(allPolygLines[i][indexStart]);
  }

  for (int i = 0; i < polygLine.size(); i++){
    for (int j = 0; j < subPolPoints.size() - 1; j++){

      if (subPolPoints[j].markedSegments.second == polygLine[i]){
        int leftIndex = i;
        int rightIndex = (i + 1) % polygLine.size();
        int indexDeleted = (i + 1) % polygLine.size();

        Segment_2 segmentInsert = Segment_2(polygLine[leftIndex].source(), polygLine[rightIndex].target());
        Segment_2 segForTriangle =  Segment_2(polygLine[rightIndex].target(), polygLine[leftIndex].source());

        std::vector<Segment_2> triangleVec = {polygLine[leftIndex], polygLine[rightIndex], segForTriangle};
        optimisedArea += std::abs(calcArea(triangleVec));

        polygLine.insert(polygLine.begin() + i, segmentInsert);
        polygLine.erase(polygLine.begin() + indexDeleted);
        polygLine.erase(polygLine.begin() + indexDeleted);
        
      }
    }
  }

}

// get marked segments for spesific subset of points
std::pair<Segment_2, Segment_2> simulatedAnnealing::getMarkedSegments(const std::vector<Point>& subPoints){

  std::pair<Segment_2, Segment_2> markedSegments;

  std::vector<Segment_2> lowerSegments =  getLowerHull(subPoints);

  markedSegments.first = lowerSegments[0];
  
  markedSegments.second = Segment_2(lowerSegments[lowerSegments.size() - 1 >= 0 ? lowerSegments.size() - 1 : 0].source(), subPoints[subPoints.size() - 1]);

  return markedSegments;

}

// Create subset point teams based on spatial search condition
void simulatedAnnealing::createSubsetPoints(std::vector<subTeam>& subPolPoints){

  int m = this->m, pIndex = 0, threshold = 6, n = points.size(), k = std::ceil((double)(n - 1)/(double)(m - 1));

  // Sorting lex order
  std::sort(points.begin(), points.end(), lexOrderPoints);

 
  subTeam firstSubSet;
  subPolPoints.push_back(firstSubSet);
  
  int pointIndex = 0;

  for (int i = 0; i < m; i++){
    subPolPoints[0].points.push_back(points[pointIndex]);
    pointIndex++;
  }
  pointIndex--;
  subPolPoints[0].markedSegments = getMarkedSegments(subPolPoints[0].points);

  
  for (int i = 1; i < k; i++){
    
    subTeam subTeam;
    subPolPoints.push_back(subTeam);

    for (int j = 0; j < m; j++){
      subPolPoints[i].points.push_back(points[pointIndex]);
      pointIndex++;
      
      if (pointIndex == n)
        break;

    }
    
    subPolPoints[i].markedSegments = getMarkedSegments(subPolPoints[i].points);
    
    int condIndex = i - 1;

  
    while (subPolPoints[i].points.size() > 0 && !spatialCondition(subPolPoints, condIndex, i)){

      subPolPoints[condIndex].points.push_back(subPolPoints[i].points[1]);
      subPolPoints[i].points.erase(subPolPoints[i].points.begin());
      
      if (pointIndex < n){
        subPolPoints[i].points.push_back(points[pointIndex]);
        pointIndex++;
      }

      subPolPoints[condIndex].markedSegments = getMarkedSegments(subPolPoints[condIndex].points);
      subPolPoints[i].markedSegments = getMarkedSegments(subPolPoints[i].points);
      
      if (pointIndex >= n)
        break;
    
    }

    if (subPolPoints[i].points.size() < threshold){
      for (int mergeIndex = 1; mergeIndex < subPolPoints[i].points.size(); mergeIndex++){
        subPolPoints[condIndex].points.push_back(subPolPoints[i].points[mergeIndex]);
      }
      subPolPoints[condIndex].markedSegments = getMarkedSegments(subPolPoints[condIndex].points);
      subPolPoints.erase(subPolPoints.begin() + i);
      break;
    }

    if (pointIndex >= n)
      break;

    pointIndex--;
  }

}

// Polygonization for every sub point
void simulatedAnnealing::subPolygonization(std::vector<subTeam>& subPolPoints, std::vector<std::vector<Segment_2>>& allPolygLines, int edge_selection){

  for (int i = 0; i < subPolPoints.size(); i++){

    if (this->subDAlgo == 1){
      try {
        
        bool markedSegmentEffect = false;
        if (initialization == "1a" && i != 0){
          markedSegmentEffect = true;
        }
        else if (initialization == "1b" && i != subPolPoints.size() - 1){
          markedSegmentEffect = true;
        }

        Incremental pol = initialization == "1a" ? Incremental(subPolPoints[i].points, edge_selection, "1a", subPolPoints[i].markedSegments.first.target(), markedSegmentEffect) : Incremental(subPolPoints[i].points, edge_selection, "1b", subPolPoints[i].markedSegments.second.source(), markedSegmentEffect);
        pol.start();

        bool leftSegmentFound = false, rightSegmentFound = false;
        for (const Segment_2& seg : pol.getPolygonLine()){
          if (looseSegCompare(seg, subPolPoints[i].markedSegments.first)){
            leftSegmentFound = true;
          }
          else if (looseSegCompare(seg, subPolPoints[i].markedSegments.second)){
            rightSegmentFound = true;
          }
        }

        if (i != 0 && i != subPolPoints.size() - 1 && (!leftSegmentFound || !rightSegmentFound)){
          throw incerementalFailure("Marked Segment Not Found");
        }
        if (i == 0 && !rightSegmentFound){
          throw incerementalFailure("Marked Segment Not Found");
        }
        if (i == subPolPoints.size() - 1 && !leftSegmentFound){
          throw incerementalFailure("Marked Segment Not Found");
        }

        allPolygLines[i] = pol.getPolygonLine();
      } 
      catch (incerementalFailure incFail) {
        
        convexHull pol = convexHull(subPolPoints[i].points, edge_selection, subPolPoints[i].markedSegments.first, subPolPoints[i].markedSegments.second, i, subPolPoints.size() - 1);
        pol.start();
        allPolygLines[i] = pol.getPolygonLine();
      }

    }
    else if (this->subDAlgo == 2){

      convexHull pol = convexHull(subPolPoints[i].points, edge_selection, subPolPoints[i].markedSegments.first, subPolPoints[i].markedSegments.second, i, subPolPoints.size() - 1);
      pol.start();
      allPolygLines[i] = pol.getPolygonLine();
    }

  }

}


// Global Transitions for every sub polygon
void simulatedAnnealing::subGlobalTransitions(std::vector<subTeam>& subPolPoints, std::vector<std::vector<Segment_2>>& allPolygLines){
  
  for (int i = 0; i < allPolygLines.size(); i++){
    
    std::cout << " i = " << i << std::endl;
    std::vector<Segment_2> markedSegments;

    if (i != 0)
      markedSegments.push_back(subPolPoints[i].markedSegments.first);
    
    if (i != allPolygLines.size() - 1)
      markedSegments.push_back(subPolPoints[i].markedSegments.second);

    simulatedAnnealing annealing = simulatedAnnealing(points , allPolygLines[i], calcArea(allPolygLines[i]), calcRatio(allPolygLines[i], calcArea(allPolygLines[i])), 1000, 2, 2, markedSegments);
    annealing.startAnnealing();

    allPolygLines[i] = annealing.getPolygonLine();
  }

}

void simulatedAnnealing::startSubdivision(){

  // Create subset teams
  std::vector<subTeam> subPolPoints;
  createSubsetPoints(subPolPoints);

  // Polygonization for every subset
  std::vector<std::vector<Segment_2>> allPolygLines(subPolPoints.size());
  subPolygonization(subPolPoints, allPolygLines, edgeSelection);

  // // Global transitions for every sub polygon
  subGlobalTransitions(subPolPoints, allPolygLines);

  // Merge all subset polygons
  mergePolygons(subPolPoints, allPolygLines);

  // Local transitions for final polygon
  simulatedAnnealing annealing = simulatedAnnealing(points , polygLine, calcArea(polygLine), calcRatio(polygLine, calcArea(polygLine)), 2, 2, 1);
  annealing.startAnnealing();
  polygLine = annealing.getPolygonLine();


  for (const Segment_2& segment : polygLine){
      std::cout << "CONSTRUCTED ==> " << segment << std::endl;
  }

  optimisedArea = calcArea(polygLine);

  std::cout << " Real Area After ==> " <<(long int) optimisedArea << std::endl;
  std::cout << "Segments Size: " << polygLine.size() << std::endl;
  std::cout << "Polygon Simplicity: " << checkPolygonSimplicity(polygLine) << std::endl;
}