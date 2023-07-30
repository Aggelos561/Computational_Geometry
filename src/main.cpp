#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include "../include/polygonization.hpp"
#include "../include/convexHull.hpp"
#include "../include/dataio.hpp"
#include "../include/incremental.hpp"
#include "../include/localSearch.hpp"
#include "../include/simulatedAnnealing.hpp"
#include "../include/preprocessor.hpp"
#include "../include/showCasedAlgos.hpp"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Epick::FT ft;


void runExecution(const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, const int , const int, const int, const int, const int);
void runScores(const std::string&, const std::string&, const bool&);


int main(int argc, char **argv){

  std::string nameOfDirectory, nameOfFile, outputFile, algorithm, algorithm_initial, initialization, annealing;
  int edge_selection_int, polygon_edge_selection, L, m;

  double threshold;
  bool preprocessEnabled, parsed;

  std::string executionType = dataio::getRunType(argc, argv);

  if (executionType == "EXECUTION"){
    // Reading parameters from command line
    parsed = dataio::getExecutionParameters(nameOfFile, outputFile, algorithm, algorithm_initial, initialization, polygon_edge_selection, edge_selection_int, threshold, annealing, L, argc, argv, m);
  }

  else{
    // Reading parameters from command line
    parsed = dataio::getScoreParameters(nameOfDirectory, outputFile, preprocessEnabled, argc, argv);
  }

  // Checking if parameters are parsed correctly
  if (!parsed){
    std::cout << "Invalid Input Parameters" << std::endl;
    return -1;
  }

  if (executionType == "EXECUTION"){
    runExecution(nameOfFile, outputFile, algorithm_initial, algorithm, annealing, initialization, edge_selection_int, polygon_edge_selection, L, m, threshold);
  }

  else{
    runScores(nameOfDirectory, outputFile, preprocessEnabled);
  }

  return 0;
}


void runExecution(const std::string& nameOfFile, const std::string& outputFile, const std::string& algorithm_initial, const std::string& algorithm, const std::string& annealing, const std::string& initialization, const int edge_selection, const int polygonEdgeSelection, const int L, const int m, const int threshold){

  std::string max_min = edge_selection == 2 ? "min" : "max";

  // Reading points from the input file
  std::vector<Point> points = dataio::readPoints(nameOfFile);
  // Calling incremental algorithm

  if (algorithm_initial == "incremental"){
    
    if (algorithm == "simulated_annealing"){ 

      if(annealing == "local"){
        Incremental pol = Incremental(points, polygonEdgeSelection, initialization);
        
        pol.start();
        ft areaBefore = pol.getArea();
        ft ratioBefore = pol.getRatio();
        
        auto start = std::chrono::high_resolution_clock::now();
        simulatedAnnealing sim = simulatedAnnealing(points, pol.getPolygonLine(), areaBefore, ratioBefore, L, edge_selection - 1, 1);
        sim.startAnnealing();
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration, outputFile, algorithm, edge_selection, initialization,max_min);
      }
      else if(annealing == "global"){
        Incremental pol = Incremental(points, polygonEdgeSelection, initialization);

        pol.start();
        ft areaBefore = pol.getArea();
        ft ratioBefore = pol.getRatio();

        auto start = std::chrono::high_resolution_clock::now();
        simulatedAnnealing sim = simulatedAnnealing(points, pol.getPolygonLine(), areaBefore, ratioBefore, L, edge_selection - 1, 2);
        sim.startAnnealing();
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration, outputFile, algorithm, edge_selection, initialization,max_min);
      }
      else{
        
        auto start = std::chrono::high_resolution_clock::now();
        
        simulatedAnnealing sim = simulatedAnnealing(points, L, polygonEdgeSelection, edge_selection - 1, initialization, m, 1);
        sim.startSubdivision();
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        ft areaBefore = sim.getArea();
        ft ratioBefore = sim.getRatio();

        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization,max_min);
      }
    }
    else if (algorithm == "local_search"){

        Incremental pol = Incremental(points, polygonEdgeSelection, initialization);

        pol.start();
        ft areaBefore = pol.getArea();
        ft ratioBefore = pol.getRatio();

        auto start = std::chrono::high_resolution_clock::now();
        localSearch local = localSearch(points, pol.getPolygonLine(), pol.getArea(), pol.getRatio(), L, threshold, edge_selection - 1);
        local.start();
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        ft areaNow = local.getOptimisedArea();
        ft ratioNow = local.getOptimisedRatio();
        
        dataio::createResultsFile(local.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization, max_min);
      }
      
    }
  // Calling convex hull algorithm
  else if (algorithm_initial == "convex_hull"){

    std::chrono::milliseconds duration;
    ft areaNow,ratioNow;

    if (algorithm == "simulated_annealing"){

      if(annealing == "local"){

        convexHull pol = convexHull(points, polygonEdgeSelection);

        pol.start();

        ft areaBefore = pol.getArea();
        ft ratioBefore = pol.getRatio();

        auto start = std::chrono::high_resolution_clock::now();
        simulatedAnnealing sim = simulatedAnnealing(points, pol.getPolygonLine(), areaBefore, ratioBefore, L, edge_selection - 1, 1);
        sim.startAnnealing();
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization,max_min);
      
      }
      else if(annealing == "global"){
        
        convexHull pol = convexHull(points, polygonEdgeSelection);

        pol.start();
        ft areaBefore = pol.getArea();
        ft ratioBefore = pol.getRatio();
        
        auto start = std::chrono::high_resolution_clock::now();
        simulatedAnnealing sim = simulatedAnnealing(points, pol.getPolygonLine(), areaBefore, ratioBefore, L, edge_selection - 1, 2);
        sim.startAnnealing();
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();
        
        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization, max_min);
      }
      else{

        simulatedAnnealing sim = simulatedAnnealing(points, L, polygonEdgeSelection, edge_selection - 1, initialization, m, 2);
        auto start = std::chrono::high_resolution_clock::now();

        sim.startSubdivision();

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        ft areaBefore = sim.getArea();
        ft ratioBefore = sim.getRatio();
        
        ft areaNow = sim.getOptimisedArea();
        ft ratioNow = sim.getOptimisedRatio();

        dataio::createResultsFile(sim.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization, max_min);
      }

    }
    else if (algorithm == "local_search"){
      convexHull pol = convexHull(points, polygonEdgeSelection);

      // Convex hull algorithm begins
      pol.start();
      ft areaBefore = pol.getArea();
      ft ratioBefore = pol.getRatio();

      auto start = std::chrono::high_resolution_clock::now();
      localSearch local = localSearch(points, pol.getPolygonLine(), pol.getArea(), pol.getRatio(), L, threshold, edge_selection - 1);
      local.start();
      
      auto stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      
      areaNow = local.getOptimisedArea();
      ratioNow = local.getOptimisedRatio();
      
      dataio::createResultsFile(local.getPolygonLine(), areaNow, areaBefore, ratioNow, ratioBefore, duration,outputFile, algorithm, edge_selection, initialization, max_min);
    }

  }
}



// Run scores
void runScores(const std::string& nameOfDirectory, const std::string& outputFile, const bool& preprocessEnabled){

  // Finds files and stores the name of every file and points size of the spesific dir
  std::vector<std::pair<int, std::string>> filesNPoints = dataio::findFiles(nameOfDirectory);
  int previousPointsSize = filesNPoints[0].first;

  std::vector<Point> points;
  
  int fileIndex = 0;
  int scoreIndex = 0;
  bool firstWrite = true;
  
  std::chrono::milliseconds cutOff;

  // Vectors for scores
  std::vector<std::pair<double, double>> scores(4);
  std::vector<std::vector<double>> boundMinScores(4);
  std::vector<std::vector<double>> boundMaxScores(4);

  // Contruct inctance preprocessor
  Preprocessor processor(filesNPoints, preprocessEnabled);
  
  // If enabled then do the preprocessing
  if (preprocessEnabled){
    showCasedAlgos::initPreprocess(points, processor, filesNPoints);
  }
  // Else just give the parameter L for every algorithm a spesific L parameter (depends on algorithm)
  else{
    processor.defaultInput();
  }

  // For every file in the given dir that was stores in the vector eariler
  // --> execute the algorithms
  // --> store scores
  // --> if the input of spesific size is completed then write the scores to ouput file
  // --> reset scores vectors

  for (const std::pair<int, std::string>& f : filesNPoints){

    // cut off
    cutOff = std::chrono::milliseconds(500 * f.first);

    fileIndex++;
    scoreIndex = 0;

    std::cout << "File: " << f.second << std::endl;
    
    showCasedAlgos::partialWrite(filesNPoints, scores, boundMinScores, boundMaxScores, f, firstWrite, previousPointsSize, outputFile, fileIndex, false);

    points = dataio::readPoints(f.second);

    showCasedAlgos::runAlgorithm("INC+GLOBAL+LOCAL", points, scores, boundMinScores, boundMaxScores, cutOff, scoreIndex, f, processor);
    std::cout << "Passed INC+GLOBAL+LOCAL"  << std::endl;

    showCasedAlgos::runAlgorithm("SUBDIVISION", points, scores, boundMinScores, boundMaxScores, cutOff, scoreIndex, f, processor);
    std::cout << "Passed SUBDIVISION"  << std::endl;

    showCasedAlgos::runAlgorithm("INC+LOCAL", points, scores, boundMinScores, boundMaxScores, cutOff, scoreIndex, f, processor);
    std::cout << "Passed INC+LOCAL"  << std::endl;

    showCasedAlgos::runAlgorithm("CONVEX+LOCAL", points, scores, boundMinScores, boundMaxScores, cutOff, scoreIndex, f, processor);
    std::cout << "Passed CONVEX+LOCAL"  << std::endl;

    if (fileIndex == filesNPoints.size())
      showCasedAlgos::partialWrite(filesNPoints, scores, boundMinScores, boundMaxScores, f, firstWrite, f.first, outputFile, fileIndex, true);

    previousPointsSize = f.first;

  }
  
}
