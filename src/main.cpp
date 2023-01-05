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
#include "../include/preprocessor.hpp"
#include "../include/showCasedAlgos.hpp"


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
  
  std::chrono::milliseconds cutOff;

  std::vector<std::pair<double, double>> scores(4);
  std::vector<std::vector<double>> boundMinScores(4);
  std::vector<std::vector<double>> boundMaxScores(4);

  Preprocessor processor;
  
  if (preprocessEnabled){
    showCasedAlgos::initPreprocess(points, processor, filesNPoints);
  }
  else{
    processor.defaultInput(points);
  }

  for (const std::pair<int, std::string>& f : filesNPoints){

    cutOff = std::chrono::milliseconds(500 * f.first);

    fileIndex++;
    scoreIndex = 0;

    std::cout << "File: " << f.second << std::endl;
    
    showCasedAlgos::partialWrite(filesNPoints, scores, boundMinScores, boundMaxScores, f, firstWrite, previousPointsSize, outputFile, fileIndex, false);

    points = dataio::readPoints(f.second);

    showCasedAlgos::runAlgorithm("INC+GLOBAL+LOCAL", points, scores, boundMinScores, boundMaxScores, cutOff, scoreIndex, f, processor);

    showCasedAlgos::runAlgorithm("SUBDIVISION", points, scores, boundMinScores, boundMaxScores, cutOff, scoreIndex, f, processor);

    showCasedAlgos::runAlgorithm("INC+LOCAL", points, scores, boundMinScores, boundMaxScores, cutOff, scoreIndex, f, processor);

    showCasedAlgos::runAlgorithm("CONVEX+LOCAL", points, scores, boundMinScores, boundMaxScores, cutOff, scoreIndex, f, processor);

    if (fileIndex == filesNPoints.size())
      showCasedAlgos::partialWrite(filesNPoints, scores, boundMinScores, boundMaxScores, f, firstWrite, f.first, outputFile, fileIndex, true);

    previousPointsSize = f.first;

  }
  
  return 0;
}
