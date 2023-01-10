#include "../include/convexHull.hpp"
#include "../include/polygonization.hpp"
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <iostream>
#include <vector>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point>::type> Convex_hull_traits_2;
typedef CGAL::Epick::FT ft;

// Constructor for convex hull class
convexHull::convexHull(const std::vector<Point> &points, int edgeSelection) : Polygonization(points, edgeSelection) {
}

convexHull::convexHull(const std::vector<Point> &points, int edgeSelection, const Segment_2& leftConnection, const Segment_2& rightConnection, int polygonIndex, int maxPolygonIndex) : Polygonization(points, edgeSelection) {
  this->leftConnection = leftConnection;
  this->rightConnection = rightConnection;
  this->spatial_Subdivision = true;
  this->polygonIndex = polygonIndex;
  this->maxPolygonIndex = maxPolygonIndex;
}



// convex hull algorithm MAIN method
void convexHull::start(const std::chrono::_V2::system_clock::time_point startTime, const std::chrono::milliseconds cutOff){

  srand(time(NULL));

  //Store remaining points
  remainingPoints = points;

  // Store init convex hull into polygon line
  initializeConvexHull(polygLine, points, remainingPoints);

  // Get polygon lone points
 	polygLinePoints = getPolyLinePoints(polygLine);

  // Get current convex hull segments
  std::vector<Segment_2> currConvexHullSegments = getConvexHull(polygLinePoints, remainingPoints);


  // Run to find the first next new point
  initialRun(currConvexHullSegments, remainingPoints, polygLine);

  // if new point on convex hull then force insert it into polygon line
  for (int i = 0; i < remainingPoints.size(); i++){

    if (forceInsertPoint(polygLine, remainingPoints[i])) {
      remainingPoints.erase(remainingPoints.begin() + i);
      i--;
    }
  
  }

  while(remainingPoints.size() > 0){

    if (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime) > cutOff){
      throw cutOffAbort("Cut off time exceeded");
    }

    std::vector<pair> bestPoints;

    //For every polyon line segment find visible and shortest point
    for(int i = 0; i < polygLine.size(); i++){

      if(this->spatial_Subdivision){
        
        if(polygonIndex != 0 && looseSegCompare(polygLine[i], this->leftConnection)){
          continue;
        }
        else if (polygonIndex != maxPolygonIndex && looseSegCompare(polygLine[i], this->rightConnection)){
          continue;
        }
      }
      

      std::vector<visPoint> allPossiblePoints;
      //Find visible points on this segment
      for(int j = 0; j < remainingPoints.size(); j++){
        double distance = squared_distance(remainingPoints[j], polygLine[i]);
        Point insPoint = remainingPoints[j];
        visPoint inserted;
        inserted.cor = insPoint;
        inserted.distance = distance;
        allPossiblePoints.push_back(inserted);
      }

      sort(allPossiblePoints.begin(), allPossiblePoints.end(), sortVisPointsAsc);

      std::vector<visPoint> visPoints;
      
      findVisiblePoints(visPoints, allPossiblePoints[0].cor, polygLine[i], polygLine);
   
      if(visPoints.size() == 0){
        continue;
      }

      Point bestPoint = visPoints[0].cor;
      
      Segment_2 bestSeg = polygLine[i];
      pair best = {bestPoint, bestSeg};
      bestPoints.push_back(best);
    }
    //Select best segmnent,point for insertion in polygon
    insertBestPoint(bestPoints, remainingPoints, polygLine);
  }

  totalArea = calcArea(polygLine);
  ratio = calcRatio(initialConvexHull, totalArea);
}


// Calculate convex hull
std::vector<Segment_2> convexHull::getConvexHull(const std::vector<Point> &polygLinePoints, std::vector<Point> &remainingPoints) {

  std::vector<Point> convexHullPoints;

  // Find Complex Hull Points
  CGAL::convex_hull_2(polygLinePoints.begin(), polygLinePoints.end(), std::back_inserter(convexHullPoints));

  Polygon_2 convexHullPolygon = Polygon_2();

  for (auto it = convexHullPoints.begin(); it != convexHullPoints.end(); ++it) {
    convexHullPolygon.push_back(*it);
    remainingPoints.erase(std::remove(remainingPoints.begin(), remainingPoints.end(), *it), remainingPoints.end());
  }

  std::vector<Segment_2> convexHullSegments;

  // Create A Vector That Stores Convex Hull segments
  for (const Segment_2& seg : convexHullPolygon.edges()) {
    convexHullSegments.push_back(seg);
  }

  return convexHullSegments;
}

// Get first convex hull and return vector of segments
void convexHull::initializeConvexHull(std::vector<Segment_2> &polygLine, const std::vector<Point> &points, std::vector<Point> &remainingPoints) {

  std::vector<Point> convexPoints;

  for (int i = 0; i < points.size(); i++)
    convexPoints.push_back(points[i]);

  polygLine = getConvexHull(convexPoints, remainingPoints);
  this->initialConvexHull = polygLine;

}

void convexHull::initialRun(const std::vector<Segment_2> &currConvexHullSegments, std::vector<Point> &remainingPoints, std::vector<Segment_2> &polygLine) {
  
  bool done = false;

  for (int i = 0; i < currConvexHullSegments.size(); i++) {
    
    std::vector<visPoint> visPoints;

    if (spatial_Subdivision && (polygonIndex != 0 && looseSegCompare(currConvexHullSegments[i], leftConnection)) || (polygonIndex != maxPolygonIndex && looseSegCompare(currConvexHullSegments[i], rightConnection)))
      continue;
    
    if (done) 
      return;
    
    for (int j = 0; j < remainingPoints.size(); j++) {
      Segment_2 segmentsArray[] = {
          Segment_2(remainingPoints[j], currConvexHullSegments[i].point(0)),
          Segment_2(remainingPoints[j], currConvexHullSegments[i].point(1)),
          Segment_2(remainingPoints[j], Point((currConvexHullSegments[i].point(0).x() + currConvexHullSegments[i].point(1).x()) / 2, ((currConvexHullSegments[i].point(0).y() + currConvexHullSegments[i].point(1).y()) / 2)))};
      
			bool flag = false;

      for (const Segment_2& polygSeg : polygLine) {
       
        if (currConvexHullSegments[i] == polygSeg)
          continue;

        for (int k = 0; k < 3; k++) {
          const auto result = intersection(segmentsArray[k], polygSeg);

          if (result) {

            if (const Point *p = boost::get<Point>(&*result)) {
              if (*p == currConvexHullSegments[i].point(0) || *p == currConvexHullSegments[i].point(1)) {
                continue;
              } else {
                flag = true;
                break;
              }
            } else {
              flag = true;
              break;
            }

          }
        }
        
        if (flag) {
          break;
        }
      }

      if (!flag) {
        double distance = squared_distance(remainingPoints[j], currConvexHullSegments[i]);
        Point insPoint = remainingPoints[j];
        visPoint inserted;
        inserted.cor = insPoint;
        inserted.distance = distance;
        visPoints.push_back(inserted);
        done = true;
      }

    }

    if (!visPoints.size())
      continue;
    
    int index = 0;
    Point bestPoint = visPoints[0].cor;
    double bestDist = visPoints[0].distance;
    
    for (int k = 0; k < visPoints.size(); k++) {
      if (bestDist > visPoints[i].distance) {
        bestPoint = visPoints[i].cor;
        bestDist = visPoints[i].distance;
      }
    }
    
    for (int m = 0; m < remainingPoints.size(); m++) {
      
      if (remainingPoints[m] == bestPoint) {
        remainingPoints.erase(remainingPoints.begin() + m);
        Segment_2 visSeg = currConvexHullSegments[i];
        deleteSegment(polygLine, visSeg);
        expandPolygonLine(polygLine, visSeg, bestPoint);
        break;
      }
    }
  }
}


// find visible points based on a vector of segments polygon line
bool convexHull::findVisiblePoints(std::vector<visPoint> &visPoints, const Point &remainingPoint, const Segment_2 &seg, const std::vector<Segment_2> &polygLine) {
  
  Segment_2 segmentsArray[] = {
      Segment_2(remainingPoint, seg.point(0)),
      Segment_2(remainingPoint, seg.point(1)),
      Segment_2(remainingPoint, Point((seg.point(0).x() + seg.point(1).x()) / 2, ((seg.point(0).y() + seg.point(1).y()) / 2)))};
  
	bool flag = false;
  for (const Segment_2& polygSeg : polygLine) {

    for (int k = 0; k < 3; k++) {
      const auto result = intersection(segmentsArray[k], polygSeg);

      if (result) {
        if (const Point *p = boost::get<Point>(&*result)) {
          if (*p == seg.point(0) || *p == seg.point(1) || *p == Point((seg.point(0).x() + seg.point(1).x()) / 2,((seg.point(0).y() + seg.point(1).y()) / 2))) {
            continue;
          } else {
            flag = true;
            break;
          }
        } else {
          flag = true;
          break;
        }
      }
    }

    if (flag) {
      break;
    }
  }

  if (!flag) {
    double distance = squared_distance(remainingPoint, seg);
    Point insPoint = remainingPoint;
    visPoint inserted;
    inserted.cor = insPoint;
    inserted.distance = distance;
    visPoints.push_back(inserted);
    return true;
  }
  else{
    visPoints.clear();
    return false;
  }

}


// Find closest point
Point convexHull::findBestPoint(const std::vector<visPoint> &visPoints) {
  int index = 0;
  Point bestPoint = visPoints[0].cor;
  double bestDist = visPoints[0].distance;

  for (int k = 0; k < visPoints.size(); k++) {
    if (bestDist > visPoints[k].distance) {
      bestPoint = visPoints[k].cor;
      bestDist = visPoints[k].distance;
    }
  }
  return bestPoint;
}

void convexHull::insertBestPoint(const std::vector<pair> &bestPoints, std::vector<Point> &remainingPoints, std::vector<Segment_2> &polygLine) {
  
	if (this->edgeSelection != 1) {
  
    if (bestPoints.size() == 0){
      throw polygonizationFailure("Convex Hull Failure");
    }
  
    std::vector<Segment_2> testPolyg = polygLine;
    deleteSegment(testPolyg, bestPoints[0].seg);

    expandPolygonLine(testPolyg, bestPoints[0].seg, bestPoints[0].cor);
    
    std::vector<Point> polygLinePoints = getPolyLinePoints(testPolyg);
    
    ft chosenArea = CGAL::polygon_area_2(polygLinePoints.begin(), polygLinePoints.end(), Convex_hull_traits_2(CGAL::make_property_map(polygLinePoints)));
    
    int index = 0;
    pair bestPair;
    bestPair.cor = bestPoints[0].cor;
    bestPair.seg = bestPoints[0].seg;
    
    for (int i = 0; i < bestPoints.size(); i++) {
      
      std::vector<Segment_2> testPolyg = polygLine;
      deleteSegment(testPolyg, bestPoints[i].seg);
      expandPolygonLine(testPolyg, bestPoints[i].seg, bestPoints[i].cor);
      std::vector<Point> polygLinePoints = getPolyLinePoints(testPolyg);
      
      ft polArea = CGAL::polygon_area_2(polygLinePoints.begin(), polygLinePoints.end(), Convex_hull_traits_2(CGAL::make_property_map(polygLinePoints)));
      
      if (this->edgeSelection == 2 && polArea < chosenArea) {
        index = i;
        chosenArea = polArea;
        bestPair.cor = bestPoints[i].cor;
        bestPair.seg = bestPoints[i].seg;
      } 
      else if (this->edgeSelection == 3 && polArea > chosenArea) {
        index = i;
        chosenArea = polArea;
        bestPair.cor = bestPoints[i].cor;
        bestPair.seg = bestPoints[i].seg;
      }
      
    }
    
    for (int m = 0; m < remainingPoints.size(); m++) {
      if (remainingPoints[m] == bestPair.cor) {

        remainingPoints.erase(remainingPoints.begin() + m);
        deleteSegment(polygLine, bestPair.seg);
        expandPolygonLine(polygLine, bestPair.seg, bestPair.cor);
        break;
        
      }
    }
  } else {
    
    if (bestPoints.size() == 0){
      throw polygonizationFailure("Convex Hull Failure");
    }

    int rindex = rand() % bestPoints.size();
    
    pair bestPair;

    bestPair.cor = bestPoints[rindex].cor;
    bestPair.seg = bestPoints[rindex].seg;
    
    for (int m = 0; m < remainingPoints.size(); m++) {
      if (remainingPoints[m] == bestPair.cor) {
        remainingPoints.erase(remainingPoints.begin() + m);
        deleteSegment(polygLine, bestPair.seg);
        expandPolygonLine(polygLine, bestPair.seg, bestPair.cor);
        break;
      }
    }
  }
}

// Sort Y Ascending
bool convexHull::sortVisPointsAsc(const visPoint& a , const visPoint& b) {
  if (a.distance < b.distance)
    return true;
  
  return false;
}