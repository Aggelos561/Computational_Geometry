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
typedef CGAL::Epick::FT ft;



int main(int argc, char **argv){

  std::string nameOfDirectory;
  bool preprocessEnabled;
  std::string outputFile;

  // Reading parameters from command line
  bool parsed = dataio::getParameters(nameOfDirectory, outputFile, preprocessEnabled, argc, argv);

  // Checking if parameters are parsed
  if (!parsed){
    std::cout << "Input parameters invalid" << std::endl;
    return -1;
  }

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
    std::cout << "CONVEX+LOCAL"  << std::endl;

    if (fileIndex == filesNPoints.size())
      showCasedAlgos::partialWrite(filesNPoints, scores, boundMinScores, boundMaxScores, f, firstWrite, f.first, outputFile, fileIndex, true);

    previousPointsSize = f.first;

  }
  
  return 0;
}
