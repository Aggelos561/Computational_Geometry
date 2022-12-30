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

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;




int main(int argc, char **argv){

  std::string nameOfDirectory;
  bool preprocessEnabled;
  std::string outputFile;


  bool parsed = dataio::getParameters(nameOfDirectory, outputFile, preprocessEnabled, argc, argv);

  if (!parsed){
    std::cout << "Input parameters invalid" << std::endl;
    return -1;
  }

  std::vector<std::pair<int, std::string>> filesNPoints = dataio::findFiles(nameOfDirectory);
  int previousPointsSize = filesNPoints[0].first;

  
  std::vector<Point> points;
  
  int fileIndex = 0;
  int scoreIndex = 0;
  bool firstWrite = true;
  std::vector<std::pair<double, double>> scores(3);
  std::vector<std::vector<double>> boundMinScores(3);
  std::vector<std::vector<double>> boundMaxScores(3);

  Preprocessor processor(points, 0.25);
      
  processor.preProcessInput();

  for (const std::pair<int, std::string>& f : filesNPoints){

    fileIndex++;
    std::cout << "FILE: " << f.second << std::endl;
    scoreIndex = 0;

    if (f.first != previousPointsSize){
      
      std::vector<double> minScores, maxScores;

      for (int i = 0; i < boundMinScores.size(); i++){
        auto minScore = *std::max_element(boundMinScores[i].begin(), boundMinScores[i].end());
        minScores.push_back(minScore);
      }

      for (int i = 0; i < boundMaxScores.size(); i++){
        auto maxScore = *std::min_element(boundMaxScores[i].begin(), boundMaxScores[i].end());
        maxScores.push_back(maxScore);
      }

      dataio::writeToOutputFile(outputFile, scores, minScores, maxScores, previousPointsSize, firstWrite);
      
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


    points = dataio::readPoints(f.second);

    for (int mode = 2; mode <= 3; mode++){

      auto start = std::chrono::high_resolution_clock::now();
      
      convexHull convex = convexHull(points, mode);
      convex.start();

      std::vector<Segment_2> polygonLine = convex.getPolygonLine();

      localSearch local = localSearch(points, polygonLine, convex.getArea(), convex.getRatio(), processor.getLocal_L(), 0, mode - 1);

      local.start();
      
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

      if (mode == 2){
        scores[scoreIndex].first += local.getOptimisedRatio();
        boundMinScores[scoreIndex].push_back(local.getOptimisedRatio());
      }
      else{
        scores[scoreIndex].second += local.getOptimisedRatio();
        boundMaxScores[scoreIndex].push_back(local.getOptimisedRatio());
      }
    }

    scoreIndex++;

    for (int mode = 2; mode <= 3; mode++){

      auto start = std::chrono::high_resolution_clock::now();
      
      Incremental incremental = Incremental(points, mode, "1a");
      incremental.start();

      std::vector<Segment_2> polygonLine = incremental.getPolygonLine();

      localSearch local = localSearch(points, polygonLine, incremental.getArea(), incremental.getRatio(), processor.getLocal_L(), 0, mode - 1);

      local.start();
      
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

      if (mode == 2){
        scores[scoreIndex].first += local.getOptimisedRatio();
        boundMinScores[scoreIndex].push_back(local.getOptimisedRatio());
      }
      else{
        scores[scoreIndex].second += local.getOptimisedRatio();
        boundMaxScores[scoreIndex].push_back(local.getOptimisedRatio());
      }

    }

    scoreIndex++;
    
    for (int mode = 2; mode <= 3; mode++){

      auto start = std::chrono::high_resolution_clock::now();

      Incremental incremental = Incremental(points, mode, "1a");
      incremental.start();

      std::vector<Segment_2> polygonLine = incremental.getPolygonLine();

      localSearch local = localSearch(points, polygonLine, incremental.getArea(), incremental.getRatio(), processor.getLocal_L(), 0, mode - 1);

      local.start();
      
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);


      std::cout << "Initial Area: " << incremental.getArea() << ", Optimized Area: " << local.getOptimisedArea() << std::endl;

      if (mode == 2){
        scores[scoreIndex].first += local.getOptimisedRatio();
        boundMinScores[scoreIndex].push_back(local.getOptimisedRatio());
      }
      else{
        scores[scoreIndex].second += local.getOptimisedRatio();
        boundMaxScores[scoreIndex].push_back(local.getOptimisedRatio());
      }

    }

    if (fileIndex == filesNPoints.size()){
      
      std::vector<double> minScores, maxScores;
      
      for (int i = 0; i < boundMinScores.size(); i++){
        auto minScore = *std::max_element(boundMinScores[i].begin(), boundMinScores[i].end());
        minScores.push_back(minScore);
      }

      for (int i = 0; i < boundMaxScores.size(); i++){
        auto maxScore = *std::min_element(boundMaxScores[i].begin(), boundMaxScores[i].end());
        maxScores.push_back(maxScore);
      }

      dataio::writeToOutputFile(outputFile, scores, minScores, maxScores, f.first, false);
    }

    previousPointsSize = f.first;


  }
  
  return 0;
}