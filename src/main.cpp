#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <chrono>
#include <iostream>
#include <vector>
#include "../include/dataio.hpp"
#include "../include/polygonization.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;



int main() {
  
  std::vector<Point> points = dataio::readPoints();

  Polygonization pol = Polygonization(points);

  pol.convexHullAlg();

  dataio::createResultsFile(pol.getPolygonLine(), pol.getArea(), pol.getDuration(), pol.getRatio());

  return 0;
}