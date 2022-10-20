#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <CGAL/property_map.h>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<
    K, CGAL::Pointer_property_map<Point>::type>
    Convex_hull_traits_2;

std::vector<Point> readPoints() {

  std::vector<Point> points{Point(2, 1), Point(1, 2), Point(3, 4), Point(2, 6),
                            Point(3, 6)};
  return points;
}

void sortPoints(std::vector<Point> &points) {
  std::sort(points.begin(), points.end(), std::greater<Point>());
}

void initializeTriangle(std::vector<Segment_2> &polygLine,
                        std::vector<Point> &points,
                        std::vector<Point> &remainingPoints) {

  for (int i = 0; i < 3; i++) {
    polygLine.push_back(Segment_2(points[i], points[i + 1]));
    remainingPoints.erase(remainingPoints.begin());
  }
}

std::vector<Point> getPolyLinePoints(std::vector<Segment_2> &polygLine) {

  std::vector<Point> linePoints;

  for (Segment_2 seg : polygLine) {

    linePoints.push_back(seg.point(0));
  }

  return linePoints;
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

  // Must be 3 (Segments Of Triangle)
  int polygLineLength = polygLine.size();

  // Main Loop Calculate Convex Hull
  // Find Red segments
  // Find visible segments
  // Expand the polygon line

  for (const Point &p : remainingPoints) {
    std::cout << "Remaining --> " << p << std::endl;
  }

  while (remainingPoints.size() > 0) {

    std::cout << std::endl;

    // Get Polygon Line Points To Compute Complex Hull
    std::vector<Point> polygLinePoints = getPolyLinePoints(polygLine);

    // Vector Instance To Store Convex Hull Result
    std::vector<std::size_t> indices(polygLinePoints.size()), convexHullPoints;
    std::iota(indices.begin(), indices.end(), 0);

    // Find Complex Hull
    CGAL::convex_hull_2(
        indices.begin(), indices.end(), std::back_inserter(convexHullPoints),
        Convex_hull_traits_2(CGAL::make_property_map(polygLinePoints)));

    for (int i = 0; i < convexHullPoints.size(); i++) {
      std::cout << "points[" << i << "] = " << convexHullPoints[i] << std::endl;
    }

    Point nextPoint = remainingPoints[0];
    std::cout << "New Point: " << nextPoint << std::endl;

    // const auto result = CGAL::intersection();

    remainingPoints.erase(remainingPoints.begin());

    // Push 2 New Segments To Create the Expand
    // The Polygon Line Based On The New Point

    // polygLine.push_back(Segment_2());
    // polygLine.push_back(Segment_2());

    polygLineLength = polygLine.size();
  }

  return 0;
}