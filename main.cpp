#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/enum.h>
#include <CGAL/intersections.h>
#include <CGAL/property_map.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <iterator>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point>::type> Convex_hull_traits_2;

std::vector<Point> readPoints() {

  std::vector<Point> points{
      Point(5, 1),   Point(8, 2),   Point(23, 4),  Point(0, 6),   Point(50, 6),
      Point(10, 4),  Point(-10, 2), Point(-20, 3)};
  return points;
  
}

void sortPoints(std::vector<Point> &points) {
  std::sort(points.begin(), points.end(), std::greater<Point>());
}


// Get all points from a vector of segments
std::vector<Point> getPolyLinePoints(std::vector<Segment_2> &polygLine) {

  std::vector<Point> linePoints;

  for (int i = 0; i < polygLine.size(); i++) {
    linePoints.push_back(polygLine[i].point(0));
  }

  return linePoints;
}

std::vector<Segment_2> getConvexHull(std::vector<Point> &polygLinePoints) {

  std::vector<Point> convexHullPoints;

  // Find Complex Hull Points
  CGAL::convex_hull_2(polygLinePoints.begin(), polygLinePoints.end(), std::back_inserter(convexHullPoints));

  Polygon_2 convexHullPolygon = Polygon_2();

  for (auto it = convexHullPoints.begin(); it != convexHullPoints.end(); ++it) {
    convexHullPolygon.push_back(*it);
  }

  std::vector<Segment_2> convexHullSegments;

  // Create A Vector That Stores Convex Hull segments
  for (Segment_2 seg : convexHullPolygon.edges()) {
    convexHullSegments.push_back(seg);
  }

  return convexHullSegments;
}

//Compares current convex hull with a convex hull that has the new point inserted
// Returns a vector with all the red segments inserted
std::vector<Segment_2> getRedSegments(std::vector<Segment_2> &currConvexHullSegments, std::vector<Segment_2> &nextConvexHullSegments) {

  std::vector<Segment_2> convexRedSegments;

  for (Segment_2 seg : currConvexHullSegments) {

    int segmentCounter = std::count(nextConvexHullSegments.begin(), nextConvexHullSegments.end(), seg);

    // If the new convex hull does not have this segment then this segment is RED
    if (!segmentCounter){
      convexRedSegments.push_back(seg);
    }
  }

  return convexRedSegments;
}

// Initialize triangle 
void initializeTriangle(std::vector<Segment_2> &polygLine, std::vector<Point> &points, std::vector<Point> &remainingPoints) {

  std::vector<Point> trianglePoints;

  for (int i = 0; i < 3; i++){
    trianglePoints.push_back(points[i]);
    remainingPoints.erase(remainingPoints.begin());
  }

  polygLine = getConvexHull(trianglePoints);
}

Segment_2 pickRandomRedSegment(std::vector<Segment_2> &visibleSegment) {
  // Pick Random red segment
  int randomRedIndex = rand() % visibleSegment.size();

  return visibleSegment[randomRedIndex];
}

void deletePolygonLineSegment(std::vector<Segment_2> &polygLine, Segment_2 &visibleSegment){
   polygLine.erase(std::remove(polygLine.begin(), polygLine.end(), visibleSegment), polygLine.end());
}

void expandPolygonLine(std::vector<Segment_2> &polygLine, Segment_2 &visibleSegment, Point &nextPoint){
  // Insert the two new segments in the right place in polygon line
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

Segment_2 findVisibleSegment(std::vector<Segment_2> &polygLine, Segment_2 &convexSegment, Point &nextPoint) {

  int segmentExists = std::count(polygLine.begin(), polygLine.end(), convexSegment);

  // Convex hull segment == polygon line segment
  if (segmentExists) {
    return convexSegment;
  }

  // Convex hull segment not on polygon line segment
  Segment_2 segment = convexSegment;
  std::cout << "convex segment is " << segment << std::endl;
  int index = 0;
  
  for (int i = 0; i < polygLine.size(); i++){
    if (polygLine[i].point(0) == segment.point(0)){
      index = i;
      std::cout <<" IM AT -> " << polygLine[i] << std::endl;
      break;
    }
  }
  
  std::vector<Segment_2> visibleSegments;

  for (int i = index; i < polygLine.size(); i++){
    
    if (polygLine[i].point(1) != segment.point(1)){
      visibleSegments.push_back(polygLine[i]);
    }
    else{
      visibleSegments.push_back(polygLine[i]);
      break;
    }
  }

  for (Segment_2 vSeg : visibleSegments){
    std::cout << "Visible ==> " << vSeg << std::endl; 
  }

  for (Segment_2 vSeg : visibleSegments) {

    Segment_2 segmentsArray[] = {Segment_2(vSeg.point(0), nextPoint),
                                 Segment_2(nextPoint, vSeg.point(1)),
                                 Segment_2(nextPoint, Point(((vSeg.point(0).x() + vSeg.point(1).x()) / 2), ((vSeg.point(0).y() + vSeg.point(1).y()) /2)))
                                };

    for (Segment_2 polygSeg : polygLine){

      if (vSeg == polygSeg) continue;

        for (int k = 0; k < 3; k++) {
          const auto result = intersection(segmentsArray[k], polygSeg);

          if (result) {

            const Point *p = boost::get<Point>(&*result);

            if (*p == polygSeg.point(0) || *p == polygSeg.point(1)) {
              continue;
            } else {
              visibleSegments.erase(std::remove(visibleSegments.begin(), visibleSegments.end(), vSeg), visibleSegments.end());
            }
          }

        }
    }

  }
  srand(time(NULL));
  int randomIndex = rand() % (visibleSegments.size() - 1 - 0 + 1) + 0;
  return visibleSegments[randomIndex];
}

void createResultsFile(std::vector<Segment_2> &polygLine){
  
  std::ofstream outdata;              
  
  int i;                         

  outdata.open("results.txt"); 
  
  if (!outdata) {               
    std::cerr << "Error: file could not be opened" << std::endl;
    exit(1);
  }

  for (i = 0; i < polygLine.size(); ++i)
    outdata << polygLine[i] << std::endl;
  
  outdata.close();

}

int main() {

  // Read Points
  std::vector<Point> points = readPoints();

  // Get Total Points From Read
  int total_points = points.size();

  // Sort points in Descending Order
  sortPoints(points);

  std::vector<Point> remainingPoints = points;

  // Polygon Line Vector
  std::vector<Segment_2> polygLine;

  // Create A Starting Triangle
  initializeTriangle(polygLine, points, remainingPoints);

  // Main Loop Calculate Convex Hull
  // Find Red segments
  // Find visible segments
  // Expand the polygon line

  for (const Point &p : remainingPoints) {
    std::cout << "Remaining --> " << p << std::endl;
  }

  // Loop until there are no remaining points left
  while (remainingPoints.size() > 0) {

    std::cout << std::endl;

    // Get Polygon Line Points To Compute Convex Hull
    std::vector<Point> polygLinePoints = getPolyLinePoints(polygLine);

    // Calculate current convex hull
    std::vector<Segment_2> currConvexHullSegments = getConvexHull(polygLinePoints);

    // Get next point and delete it from remaining
    Point nextPoint = remainingPoints[0];
    std::cout << "New Point: " << nextPoint << std::endl;
    remainingPoints.erase(remainingPoints.begin());

    // Push next point
    polygLinePoints.push_back(nextPoint);

    // Calculate convex hull with the new point inside
    std::vector<Segment_2> nextConvexHullSegments = getConvexHull(polygLinePoints);

    for (int i = 0; i < nextConvexHullSegments.size(); i++) {
      std::cout << "New Convex Segment: " << nextConvexHullSegments[i] << std::endl;
    }

    // To get red segments compare current convex hull and convex hull with the
    // new point if the new novex hull DOES NOT have any segments from the
    // current convex hull then these segments are RED
    std::vector<Segment_2> redSegments = getRedSegments(currConvexHullSegments, nextConvexHullSegments);

    for (int i = 0; i < redSegments.size(); i++) {
      std::cout << "Red Segment: " << redSegments[i] << std::endl;
    }

    Segment_2 visibleConvexSegment = pickRandomRedSegment(redSegments);
    
    Segment_2 visibleSegment = findVisibleSegment(polygLine, visibleConvexSegment, nextPoint);

    deletePolygonLineSegment(polygLine, visibleSegment);
    
    expandPolygonLine(polygLine, visibleSegment, nextPoint);

    // Print all segments of the current polygon
    std::cout << std::endl;
    for (Segment_2 line : polygLine) {
      std::cout << "Line --> " << line << std::endl;
    }

    std::cout << std::endl;
  }

  createResultsFile(polygLine);

  return 0;
}