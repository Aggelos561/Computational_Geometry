#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <algorithm>
#include <chrono>
#include <iostream>
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




void  showCasedAlgos::runAlgorithm(const std::string& algorithm, const std::vector<Point>& points, std::vector<std::pair<double, double>>& scores, std::vector<std::vector<double>>& boundMinScores, std::vector<std::vector<double>>& boundMaxScores, const std::chrono::milliseconds cutOff, int& scoreIndex, const std::pair<int, std::string>& f, Preprocessor& processor){

  for (int mode = 2; mode <= 3; mode++){
    
    ft optimizedArea, optimizedRatio;
    std::chrono::milliseconds duration;

    if (algorithm == "INC+GLOBAL+LOCAL"){

      auto start = std::chrono::high_resolution_clock::now();

      Incremental incremental = Incremental(points, mode, "1a");
      incremental.start();

      std::vector<Segment_2> polygonLine = incremental.getPolygonLine();

      simulatedAnnealing simulatedGlobal = simulatedAnnealing(points, incremental.getPolygonLine(), incremental.getArea(), incremental.getRatio(), processor.getSimGlobal_L(f.first), mode - 1, 2);
      simulatedGlobal.startAnnealing();
      
      simulatedAnnealing simulatedLocal = simulatedAnnealing(points, simulatedGlobal.getPolygonLine(), simulatedGlobal.getOptimisedArea(), simulatedGlobal.getOptimisedRatio(), processor.getSimLocal_L(f.first), mode - 1, 1);
      simulatedLocal.startAnnealing();

      auto stop = std::chrono::high_resolution_clock::now();

      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

      optimizedArea = simulatedLocal.getOptimisedArea();
      optimizedRatio = simulatedLocal.getOptimisedRatio();

    }
    else if (algorithm == "SUBDIVISION"){

      auto start = std::chrono::high_resolution_clock::now();
      
      simulatedAnnealing simulatedSubdivision = simulatedAnnealing(points, processor.getSimSubDiv_L(f.first), mode, mode - 1, "1a", processor.getSimSubDiv_M(f.first), 1);
      simulatedSubdivision.startSubdivision();

      auto stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

      optimizedArea = simulatedSubdivision.getOptimisedArea();
      optimizedRatio = simulatedSubdivision.getOptimisedRatio();
    }
    else if (algorithm == "INC+LOCAL"){
      
      auto start = std::chrono::high_resolution_clock::now();
  
      Incremental incremental = Incremental(points, mode, "1a");
      incremental.start();

      std::vector<Segment_2> polygonLine = incremental.getPolygonLine();

      simulatedAnnealing simulatedLocal = simulatedAnnealing(points, incremental.getPolygonLine(), incremental.getArea(), incremental.getRatio(), processor.getSimLocal_L(f.first), mode - 1, 1);
      simulatedLocal.startAnnealing();

      auto stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

      optimizedArea = simulatedLocal.getOptimisedArea();
      optimizedRatio = simulatedLocal.getOptimisedRatio();
    
    }
    else if (algorithm == "CONVEX+LOCAL"){

      auto start = std::chrono::high_resolution_clock::now();
  
      convexHull convex = convexHull(points, mode);
      convex.start();

      std::vector<Segment_2> polygonLine = convex.getPolygonLine();

      simulatedAnnealing simulatedLocal = simulatedAnnealing(points, convex.getPolygonLine(), convex.getArea(), convex.getRatio(), processor.getSimLocal_L(f.first), mode - 1, 1);
      simulatedLocal.startAnnealing();

      auto stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

      optimizedArea = simulatedLocal.getOptimisedArea();
      optimizedRatio = simulatedLocal.getOptimisedRatio();

    }

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



void showCasedAlgos::partialWrite(const std::vector<std::pair<int, std::string>>& filesNPoints, std::vector<std::pair<double, double>>& scores, std::vector<std::vector<double>>& boundMinScores, std::vector<std::vector<double>>& boundMaxScores, const std::pair<int, std::string>& f,  bool& firstWrite, int pointSize, const std::string& outputFile, int fileIndex, bool lastTime){

  if ( (!lastTime && f.first != pointSize) || lastTime){
      
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
      boundMinScores[i].clear();
    }

  }

}


void showCasedAlgos::initPreprocess(std::vector<Point> points, Preprocessor& processor, const std::vector<std::pair<int, std::string>>& filesNPoints){

  processor.defaultInput(points);
  
  for (int i = 0; i < 5; i++){
    
    int size = 0;

    for (const std::pair<int, std::string>& f : filesNPoints){

      if (i == 0){
        if (f.first <= 10){
          size++;

          if (size > 3)
            break;

          points = dataio::readPoints(f.second);
              
          processor.preprocessInput(points);
        }
      }

      if (i == 1){
        if (f.first <= 100 && f.first > 10){
          size++;

          if (size > 3)
            break;

          points = dataio::readPoints(f.second);
          processor.preprocessInput(points);
        }
      }
      if (i == 2){
        if (f.first <= 1000 && f.first > 100){
          size++;

          if (size > 3)
            break;

          points = dataio::readPoints(f.second);
          processor.preprocessInput(points);
        }
      }
      if (i == 3){
        if (f.first <= 10000 && f.first > 1000){
          size++;

          if (size > 3)
            break;

          points = dataio::readPoints(f.second);
          processor.preprocessInput(points);
        }
      }

      if (i == 4){
        if (f.first <= 100000 && f.first > 10000){
          size++;

          if (size > 3)
            break;

          points = dataio::readPoints(f.second);
          processor.preprocessInput(points);
        }
      }

    }
  }

}

