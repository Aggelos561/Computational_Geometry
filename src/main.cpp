#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <chrono>
#include <iostream>
#include <vector>
#include "../include/polygonization.hpp"
#include "../include/convexHull.hpp"
#include "../include/dataio.hpp"
#include "../include/incremental.hpp"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;



int main() {
  
  std::vector<Point> points = dataio::readPoints();

  convexHull pol = convexHull(points);

  auto start = std::chrono::high_resolution_clock::now();

  pol.start();
  
  auto stop = std::chrono::high_resolution_clock::now();

  std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  dataio::createResultsFile(pol.getPolygonLine(), pol.getArea(), duration, pol.getRatio());

  return 0;
}