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

  std::vector<Point> points{Point(2, 1), Point(1, 2), Point(3, 4), Point(2, 6),
                            Point(3, 6), Point(0, 4), Point(-1, 2)};
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
std::vector<Segment_2> getRedSegments(std::vector<Segment_2>& currConvexHullSegments, std::vector<Segment_2>& nextConvexHullSegments) {

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

Segment_2 pickRandomRedSegment(std::vector<Segment_2> &redSegments) {
  // Pick Random red segment
  int randomRedIndex = rand() % redSegments.size();

  return redSegments[randomRedIndex];
}

void deletePolygonLineSegment(std::vector<Segment_2> &polygLine, Segment_2 &redSegment){
   polygLine.erase(std::remove(polygLine.begin(), polygLine.end(), redSegment), polygLine.end());
}

void expandPolygonLine(std::vector<Segment_2> &polygLine, Segment_2 &redSegment, Point &nextPoint){
    // Insert the two new segments in the right place in polygon line
    int index = 0;
    for (int i = 0; i < polygLine.size(); i++) {
      if (polygLine[i].point(1) == redSegment.point(0)){
        
        polygLine.insert(polygLine.begin() + index + 1,Segment_2(redSegment.point(0), nextPoint));
        polygLine.insert(polygLine.begin() + index + 2, Segment_2(nextPoint, redSegment.point(1)));
        break;
        
      } 
      else if (polygLine[i].point(1) == redSegment.point(1)) {
        polygLine.insert(polygLine.begin() + index + 1,Segment_2(redSegment.point(1), nextPoint));
        polygLine.insert(polygLine.begin() + index + 2, Segment_2(nextPoint, redSegment.point(0)));
        break;
      }
      index++;
    }
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

    // Delete previous segment and connect the two new segments with the new point
    Segment_2 redSegment = pickRandomRedSegment(redSegments);

    deletePolygonLineSegment(polygLine, redSegment);
    
    expandPolygonLine(polygLine, redSegment, nextPoint);

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