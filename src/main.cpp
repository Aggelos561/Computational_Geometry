#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <chrono>
#include <iostream>
#include <vector>
#include "../include/polygonization.hpp"
#include "../include/convexHull.hpp"
#include "../include/dataio.hpp"
#include "../include/incremental.hpp"
#include "../include/localSearch.hpp"
#include "../include/simulatedAnnealing.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;



int main(int argc, char** argv) {
  
  std::string nameOfFile;
  std::string outputFile;
  std::string algorithm_initial;
  std::string algorithm;
  int annealing;
  int edge_selection;
  std::string initialization;
  double threshold;
  int L;


  // Reading and parsing parameters
  //bool parsed = dataio::getParameters(nameOfFile, outputFile, algorithm, edge_selection, initialization, argc, argv);
  bool parsed = dataio::getParameters(nameOfFile, outputFile, algorithm, algorithm_initial, initialization, edge_selection, threshold,  annealing, L, argc, argv);
  std::cout << "nameOfFile " << nameOfFile << " " << "outputFile " << outputFile << " algorithm " << algorithm << " algorith_initial " << algorithm_initial << " annealing " << annealing << " edge_selection " << edge_selection << " initialization " << initialization << " threshold " << threshold << " L " << L << std::endl;
  if (!parsed){
    std::cout << "Input parameters invalid" << std::endl;
    return 1;
  }

  // Reading points from the input file
  std::vector<Point> points = dataio::readPoints(nameOfFile);
  // Calling incremental algorithm
  if(algorithm == "incremental"){
    Incremental pol = Incremental(points, edge_selection, initialization);

    pol.start();

    auto start = std::chrono::high_resolution_clock::now();


    simulatedAnnealing sim = simulatedAnnealing(points, pol.getPolygonLine(), pol.getArea(), pol.getRatio(), L, edge_selection, annealing); // l mode transition
   
    sim.startAnnealing();

    auto stop = std::chrono::high_resolution_clock::now();
  }

  // Calling convex hull algorithm
  else if(algorithm == "convex_hull"){
    convexHull pol = convexHull(points,edge_selection);

    auto start = std::chrono::high_resolution_clock::now();

    // Convex hull algorithm begins
    pol.start();

    simulatedAnnealing sim = simulatedAnnealing(points, pol.getPolygonLine(), pol.getArea(), pol.getRatio(), L, edge_selection, annealing); // l mode transition
   
    sim.startAnnealing();

    auto stop = std::chrono::high_resolution_clock::now();

    // Calulating duration of algorithm
    std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    // Write data into result file

    dataio::createResultsFile(pol.getPolygonLine(), pol.getArea(), duration, pol.getRatio(), outputFile, algorithm, edge_selection, initialization);
  }

  return 0;
}