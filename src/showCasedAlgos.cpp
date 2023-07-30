#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <chrono>
#include <iostream>
#include <limits>
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
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;



// Runs spesific algorithm
// 4 algorithms are showcased for project 3
// --> incremental with global step and local step optimization
// --> Spatial Subdivision (Runs incremental polygonization and only if it fails its going to run convex hull)
// --> incremental with local step optimization
// --> convex hull with local step optimization
void showCasedAlgos::runAlgorithm(const std::string& algorithm, const std::vector<Point>& points, std::vector<std::pair<double, double>>& scores, std::vector<std::vector<double>>& boundMinScores, std::vector<std::vector<double>>& boundMaxScores, const std::chrono::milliseconds cutOff, int& scoreIndex, const std::pair<int, std::string>& f, Preprocessor& processor){

  // for mode 2 and 3 (min and max modes) execute the 4 algorithms based on algorithm parameter of function
  for (int mode = 2; mode <= 3; mode++){
    
    ft optimizedArea = 0, optimizedRatio = 0;
    std::chrono::milliseconds duration;

    try{

      if (algorithm == "INC+GLOBAL+LOCAL"){

        auto start = std::chrono::high_resolution_clock::now();

        Incremental incremental = Incremental(points, mode, "1a");
        incremental.start(start, cutOff, true);

        std::vector<Segment_2> polygonLine = incremental.getPolygonLine();

        simulatedAnnealing simulatedGlobal = simulatedAnnealing(points, incremental.getPolygonLine(), incremental.getArea(), incremental.getRatio(), processor.getSimGlobal_L(f.first), mode - 1, 2);
        simulatedGlobal.startAnnealing(start, cutOff, true);

        simulatedAnnealing simulatedLocal = simulatedAnnealing(points, simulatedGlobal.getPolygonLine(), simulatedGlobal.getOptimisedArea(), simulatedGlobal.getOptimisedRatio(), processor.getSimLocal_L(f.first), mode - 1, 1);
        simulatedLocal.startAnnealing(start, cutOff, true);

        auto stop = std::chrono::high_resolution_clock::now();

        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

        optimizedArea = simulatedLocal.getOptimisedArea();
        optimizedRatio = simulatedLocal.getOptimisedRatio();

      }
      else if (algorithm == "SUBDIVISION"){

        auto start = std::chrono::high_resolution_clock::now();
        
        simulatedAnnealing simulatedSubdivision = simulatedAnnealing(points, processor.getSimSubDiv_L(f.first), mode, mode - 1, "1a", processor.getSimSubDiv_M(f.first), 1);
        simulatedSubdivision.startSubdivision(start, cutOff, true);

        auto stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

        optimizedArea = simulatedSubdivision.getOptimisedArea();
        optimizedRatio = simulatedSubdivision.getOptimisedRatio();
      }
      else if (algorithm == "INC+LOCAL"){
        
        auto start = std::chrono::high_resolution_clock::now();
    
        Incremental incremental = Incremental(points, mode, "1a");
        incremental.start(start, cutOff, true);

        std::vector<Segment_2> polygonLine = incremental.getPolygonLine();

        simulatedAnnealing simulatedLocal = simulatedAnnealing(points, incremental.getPolygonLine(), incremental.getArea(), incremental.getRatio(), processor.getSimLocal_L(f.first), mode - 1, 1);
        simulatedLocal.startAnnealing(start, cutOff, true);

        auto stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

        optimizedArea = simulatedLocal.getOptimisedArea();
        optimizedRatio = simulatedLocal.getOptimisedRatio();
      
      }
      else if (algorithm == "CONVEX+LOCAL"){

        auto start = std::chrono::high_resolution_clock::now();

        try{
          
          convexHull convex = convexHull(points, mode);
          convex.start(start, cutOff, true);

          std::vector<Segment_2> polygonLine = convex.getPolygonLine();

          simulatedAnnealing simulatedLocal = simulatedAnnealing(points, convex.getPolygonLine(), convex.getArea(), convex.getRatio(), processor.getSimLocal_L(f.first), mode - 1, 1);
          simulatedLocal.startAnnealing(start, cutOff, true);

          auto stop = std::chrono::high_resolution_clock::now();
          duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

          optimizedArea = simulatedLocal.getOptimisedArea();
          optimizedRatio = simulatedLocal.getOptimisedRatio();
        }
        catch(polygonizationFailure convexFail){
          duration = std::chrono::milliseconds(500 * points.size() + 1);
        }

      }

    }
    // in case of cut off exceeding then the algorithm throws an exception so the scores are stored as 0 and 1
    catch(cutOffAbort abort){ 
      duration = std::chrono::milliseconds(500 * points.size() + 1);
    }

    // Store min, max scores and min, max bounds of them (depends on mode min max)
    if (mode == 2){
      scores[scoreIndex].first += duration > cutOff ? 1 : optimizedRatio;
      boundMinScores[scoreIndex].push_back(duration > cutOff ? 1 : optimizedRatio);
    }
    else{
      scores[scoreIndex].second += duration > cutOff ? 0 : optimizedRatio;
      boundMaxScores[scoreIndex].push_back(duration > cutOff ? 0 : optimizedRatio);
    }

  }

  scoreIndex++;
}


// Partial writing to output file (scores for the final board) every time a spesific size of points is completed
// Also resets scores and bound scores vectors
void showCasedAlgos::partialWrite(const std::vector<std::pair<int, std::string>>& filesNPoints, std::vector<std::pair<double, double>>& scores, std::vector<std::vector<double>>& boundMinScores, std::vector<std::vector<double>>& boundMaxScores, const std::pair<int, std::string>& f,  bool& firstWrite, int pointSize, const std::string& outputFile, int fileIndex, bool lastTime){

  if ((!lastTime && f.first != pointSize) || lastTime){
      
    std::vector<double> minScores, maxScores;

    for (int i = 0; i < boundMinScores.size(); i++){
      auto minScore = *std::max_element(boundMinScores[i].begin(), boundMinScores[i].end());
      minScores.push_back(minScore);
    }

    for (int i = 0; i < boundMaxScores.size(); i++){
      auto maxScore = *std::min_element(boundMaxScores[i].begin(), boundMaxScores[i].end());
      maxScores.push_back(maxScore);
    }

    dataio::writeToOutputFile(outputFile, scores, minScores, maxScores, pointSize, firstWrite);
    
    firstWrite = false;
    
    for (int i = 0; i < scores.size(); i++){
      scores[i].first = 0;
      scores[i].second = 0;
    }

    for (int i = 0; i < boundMinScores.size(); i++){
      boundMinScores[i].clear();
    }
    for (int i = 0; i < boundMaxScores.size(); i++){
      boundMaxScores[i].clear();
    }

  }

}

// executes preprocess for every input file of a the given dir 
void showCasedAlgos::initPreprocess(std::vector<Point> points, Preprocessor& processor, const std::vector<std::pair<int, std::string>>& filesNPoints){

  for (const std::pair<int, std::string>& f : filesNPoints){
    processor.preprocessInput(dataio::readPoints(f.second));
  }
  
}

