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



int main(int argc, char** argv) {
  
  std::string nameOfFile;
  std::string outputFile;
  std::string algorithm;
  int edge_selection;
  std::string initialization;

  // Reading and parsing parameters
  bool parsed = dataio::getParameters(nameOfFile, outputFile, algorithm, edge_selection, initialization, argc, argv);

  if (!parsed){
    std::cout << "Input parameters invalid" << std::endl;
    return 1;
  }

  // Reading points from the input file
  std::vector<Point> points = dataio::readPoints(nameOfFile);
  // Calling incremental algorithm
  if(algorithm == "incremental"){

    Polygonization pol = Polygonization(points, edge_selection);

    pol.spatialSubdivision(points, edge_selection, initialization);
    // Incremental pol = Incremental(points, edge_selection, initialization);

    // auto start = std::chrono::high_resolution_clock::now();

    // Incremental algorithm begins
    // pol.start();

    // auto stop = std::chrono::high_resolution_clock::now();

    //calculating duration of algorithm
    // std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    // std::vector<Segment_2> polygonLine = pol.getPolygonLine();

    // pol.localSearch(polygonLine);
    
    // pol.simulatedAnnealing(polygonLine);


    // Write data into result file
    // dataio::createResultsFile(pol.getPolygonLine(), pol.getArea(), duration, pol.getRatio(), outputFile, algorithm, edge_selection, initialization);
  }

  // Calling convex hull algorithm
  else if(algorithm == "convex_hull"){
    convexHull pol = convexHull(points,edge_selection);

    auto start = std::chrono::high_resolution_clock::now();

    // Convex hull algorithm begins
    pol.start();

    auto stop = std::chrono::high_resolution_clock::now();

    // Calulating duration of algorithm
    std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    // Write data into result file

    dataio::createResultsFile(pol.getPolygonLine(), pol.getArea(), duration, pol.getRatio(), outputFile, algorithm, edge_selection, initialization);
  }

  return 0;
}