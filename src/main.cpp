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

int main(int argc, char **argv)
{

  std::string nameOfFile;
  std::string outputFile;
  std::string algorithm_initial;
  std::string algorithm;
  int m;
  int annealing;
  int edge_selection;
  std::string initialization;
  double threshold;
  int L;

  // Reading and parsing parameters
  // bool parsed = dataio::getParameters(nameOfFile, outputFile, algorithm, edge_selection, initialization, argc, argv);
  bool parsed = dataio::getParameters(nameOfFile, outputFile, algorithm, algorithm_initial, initialization, edge_selection, threshold, annealing, L, argc, argv,m);
  std::cout << "nameOfFile " << nameOfFile << " "
            << "outputFile " << outputFile << " algorithm " << algorithm << " algorith_initial " << algorithm_initial << " annealing " << annealing << " edge_selection " << edge_selection << " initialization " << initialization << " threshold " << threshold << " L " << L << std::endl;
  if (!parsed)
  {
    std::cout << "Input parameters invalid" << std::endl;
    return 1;
  }

  // Reading points from the input file
  std::vector<Point> points = dataio::readPoints(nameOfFile);
  // Calling incremental algorithm
  if (algorithm_initial == "incremental")
  {
    Incremental pol = Incremental(points, edge_selection, initialization);
    auto start = std::chrono::high_resolution_clock::now();

    pol.start();
    ft areaBefore = pol.getArea();
    ft ratioBefore = pol.getRatio();
    std::cout << "Before " << areaBefore << " " << ratioBefore << std::endl;
    std::chrono::milliseconds duration;
    ft areaNow,ratioNow;
    if (algorithm == "simulated_annealing")
    { int mode = 1;
      simulatedAnnealing sim = simulatedAnnealing(points, L, edge_selection, edge_selection - 1, initialization, mode,m);
      sim.startSubdivision();

      auto stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      areaNow = sim.getArea();
      ratioNow = sim.getRatio(pol.getCHArea());
    }
    else if (algorithm == "local_search")
    {
      localSearch local = localSearch(points, pol.getPolygonLine(), pol.getArea(), pol.getRatio(), L, threshold, edge_selection - 1);
      local.start();
      auto stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      areaNow = local.getArea();
      ratioNow = local.getRatio(pol.getCHArea());
    }
    dataio::createResultsFile(pol.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization);
  }
  // Calling convex hull algorithm
  else if (algorithm_initial == "convex_hull")
  {
    convexHull pol = convexHull(points, edge_selection);

    auto start = std::chrono::high_resolution_clock::now();

    // Convex hull algorithm begins
    pol.start();
    ft areaBefore = pol.getArea();
    ft ratioBefore = pol.getRatio();
    std::cout << "Before " << areaBefore << " " << ratioBefore << std::endl;
    std::chrono::milliseconds duration;
    ft areaNow,ratioNow;
    if (algorithm == "simulated_annealing")
    { 
      int mode = 2;
      simulatedAnnealing sim = simulatedAnnealing(points, L, edge_selection, edge_selection - 1, initialization, mode,m);
      sim.startSubdivision();


      auto stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      areaNow = sim.getArea();
      ratioNow = sim.getRatio(pol.getCHArea());
    }
    else if (algorithm == "local_search")
    {
      localSearch local = localSearch(points, pol.getPolygonLine(), pol.getArea(), pol.getRatio(), L, threshold, edge_selection - 1);
      local.start();
      auto stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      areaNow = local.getArea();
      ratioNow = local.getRatio(pol.getCHArea());
    }
    
    dataio::createResultsFile(pol.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization);
  }

  return 0;
}