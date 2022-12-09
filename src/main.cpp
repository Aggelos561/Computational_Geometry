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
  std::string annealing;
  int edge_selection;
  std::string initialization;
  double threshold;
  int L;
  std::string max_min;
  
  // Reading and parsing parameters
  // bool parsed = dataio::getParameters(nameOfFile, outputFile, algorithm, edge_selection, initialization, argc, argv);
  bool parsed = dataio::getParameters(nameOfFile, outputFile, algorithm, algorithm_initial, initialization, edge_selection, threshold, annealing, L, argc, argv,m);
  if(edge_selection == 2){
    max_min = "min";
  }
  else{
    max_min = "max";
  }
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
    
    if (algorithm == "simulated_annealing")
    { int mode = 1;
      if(annealing == "local"){
        Incremental pol = Incremental(points, edge_selection, initialization);
        auto start = std::chrono::high_resolution_clock::now();

        pol.start();
        ft areaBefore = pol.getArea();
        ft ratioBefore = pol.getRatio();
      
        simulatedAnnealing sim = simulatedAnnealing(points, pol.getPolygonLine(), areaBefore, ratioBefore, L, edge_selection - 1,1);
        sim.startAnnealing();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization,max_min);
      }
      else if(annealing == "global"){
        Incremental pol = Incremental(points, edge_selection, initialization);
        auto start = std::chrono::high_resolution_clock::now();

        pol.start();
        ft areaBefore = pol.getArea();
        ft ratioBefore = pol.getRatio();
        simulatedAnnealing sim = simulatedAnnealing(points, pol.getPolygonLine(), areaBefore, ratioBefore, L, edge_selection - 1,2);
        sim.startAnnealing();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization,max_min);
      }
      else{
        auto start = std::chrono::high_resolution_clock::now();
        simulatedAnnealing sim = simulatedAnnealing(points, L, edge_selection, edge_selection - 1, initialization, m,1);
        sim.startSubdivision();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        ft areaBefore = sim.getArea();
        ft ratioBefore = sim.getRatio();
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization,max_min);
      }
    }
    else if (algorithm == "local_search")
    { 
      Incremental pol = Incremental(points, edge_selection, initialization);
      auto start = std::chrono::high_resolution_clock::now();

      pol.start();
      ft areaBefore = pol.getArea();
      ft ratioBefore = pol.getRatio();
      std::cout << "Before " << areaBefore << " " << ratioBefore << std::endl;
      std::chrono::milliseconds duration;
      ft areaNow,ratioNow;
      localSearch local = localSearch(points, pol.getPolygonLine(), pol.getArea(), pol.getRatio(), L, threshold, edge_selection - 1);
      local.start();
      auto stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      areaNow = local.getOptimisedArea();
      ratioNow = local.getOptimisedRatio();
      dataio::createResultsFile(local.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization,max_min);
    }
    
  }
  // Calling convex hull algorithm
  else if (algorithm_initial == "convex_hull")
  {
    std::chrono::milliseconds duration;
    ft areaNow,ratioNow;
    if (algorithm == "simulated_annealing")
    { 
      if(annealing == "local"){
        convexHull pol = convexHull(points, edge_selection);
        auto start = std::chrono::high_resolution_clock::now();

        pol.start();
        ft areaBefore = pol.getArea();
        ft ratioBefore = pol.getRatio();
      
        simulatedAnnealing sim = simulatedAnnealing(points, pol.getPolygonLine(), areaBefore, ratioBefore, L, edge_selection - 1,1);
        sim.startAnnealing();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization,max_min);
      }
      else if(annealing == "global"){
        convexHull pol = convexHull(points, edge_selection);
        auto start = std::chrono::high_resolution_clock::now();

        pol.start();
        ft areaBefore = pol.getArea();
        ft ratioBefore = pol.getRatio();
        simulatedAnnealing sim = simulatedAnnealing(points, pol.getPolygonLine(), areaBefore, ratioBefore, L, edge_selection - 1,2);
        sim.startAnnealing();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization,max_min);
      }
      else{
        simulatedAnnealing sim = simulatedAnnealing(points, L, edge_selection, edge_selection - 1, initialization, m,2);
        auto start = std::chrono::high_resolution_clock::now();
        sim.startSubdivision();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        ft areaBefore = sim.getArea();
        ft ratioBefore = sim.getRatio();
        std::cout << ratioBefore << " " << ratioNow << std::endl;
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization,max_min);
      }


      // auto stop = std::chrono::high_resolution_clock::now();
      // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      // areaNow = sim.getArea();
      // ratioNow = sim.getRatio(pol.getCHArea());
    }
    else if (algorithm == "local_search")
    {
      convexHull pol = convexHull(points, edge_selection);

      auto start = std::chrono::high_resolution_clock::now();

      // Convex hull algorithm begins
      pol.start();
      ft areaBefore = pol.getArea();
      ft ratioBefore = pol.getRatio();
      localSearch local = localSearch(points, pol.getPolygonLine(), pol.getArea(), pol.getRatio(), L, threshold, edge_selection - 1);
      local.start();
      
      auto stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      areaNow = local.getOptimisedArea();
      ratioNow = local.getOptimisedRatio();
      dataio::createResultsFile(local.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization,max_min);
    }
  }

  return 0;
}