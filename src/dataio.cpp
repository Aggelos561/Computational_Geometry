#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <dirent.h>
#include <iomanip>
#include "../include/dataio.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;


// Get parameters
bool dataio::getParameters(std::string& nameOfDirectory, std::string& outputFile, bool& preprocessEnabled, int argc, char** argv) {
  
  preprocessEnabled = false;
  for(int i = 0; i < argc- 1; i++){

    std::string param = (argv[i]);
    std::string mode = (argv[i+1]);

    if(param == "-i"){
      nameOfDirectory = mode;
    }
    else if(param == "-o"){
      outputFile = mode;
    }
    else if(mode == "-preprocess"){
      preprocessEnabled = true;
    }
  }

  return true;
}


std::vector<std::pair<int, std::string>> dataio::findFiles(const std::string& dirName){

  std::vector<std::pair<int, std::string>> files;

  DIR* directory = opendir(dirName.c_str());

  if (directory == NULL){
    std::cout << "Failed to open directory " << dirName << std::endl;
    std::exit(EXIT_FAILURE);
  }

  struct dirent* entryDir;

  while((entryDir = readdir(directory)) != NULL){
    if (entryDir->d_type != DT_DIR){
      
      std::string fileName = dirName + entryDir->d_name;

      
      std::string firstLine;

      std::ifstream pointsFile(fileName);
      std::getline(pointsFile, firstLine);

      std::stringstream streamString(firstLine);

      std::getline(streamString, firstLine, '(');
      std::getline(streamString, firstLine, ' ');

      int numberOfPoints = atoi(firstLine.c_str());

      std::cout << "Number OF Points " << numberOfPoints << std::endl;
      files.push_back(std::pair<int, std::string>(numberOfPoints, dirName + entryDir->d_name));

      pointsFile.close();

    }
  }

  closedir(directory);

  std::sort(files.begin(), files.end());

  return files;

}



// Read from data.instance file
std::vector<Point> dataio::readPoints(const std::string& name) {

  std::ifstream inFile(name.c_str());
  std::string strInput;

  std::vector<Point> points;
  
  int currentIndex= 0;

  while (std::getline(inFile, strInput)) {

    if (strInput[0] == '#' || strInput.size() < 1)
      continue;
    
    int fileIndex, x, y;
    
    std::vector<std::string> strings;
    std::string s;
    std::istringstream f(strInput);
    
    while (std::getline(f, s, '\t')) {
      strings.push_back(s);
    }
    
    fileIndex = stoi(strings[0]);
    
    if (fileIndex != currentIndex)
      break;

    x = stoi(strings[1]);
    y = stoi(strings[2]);
    currentIndex++;

    points.push_back(Point(x,y));
  }

  inFile.close();

  return points;
}



// Write resulta data into a spesific file
void dataio::writeToOutputFile(const std::string& output, const std::vector<std::pair<double, double>>& scores, const std::vector<double>& minScores, const std::vector<double>& maxScores, int filePointsSize, bool initOutput) {

  std::ofstream outdata;

  if (initOutput)
    outdata.open(output.c_str(), std::ofstream::out);
  else
    outdata.open(output.c_str(), std::ofstream::app);


  if (!outdata) {
    std::cout << "Error: file could not be opened" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string spaces = "";

  for (int i = 0; i < 18; i++){
    spaces += " ";
  }

  if (initOutput){
    outdata  << std::setw(10) << "||" << spaces << "Algorithm 1" <<  spaces << "||" << spaces << "Algorithm 2" << spaces << "||" << spaces << "Algorithm 3" << spaces << "||" << spaces << "Algorithm 4" << spaces << "||" << std::endl;
  
    outdata  << std::setw(8) << "Size " << "||" << std::setw(2) << " min score" << " | " << "max score" << " | " << "min bound" << " | " << "max bound";
    outdata << " ||" << std::setw(2) << " min score" << " | " << "max score" << " | " << "min bound" << " | " << "max bound";
    outdata << " ||" << std::setw(2) << " min score" << " | " << "max score" << " | " << "min bound" << " | " << "max bound";
    outdata << " ||" << std::setw(2) << " min score" << " | " << "max score" << " | " << "min bound" << " | " << "max bound" << " ||" << std::endl;
    
  }
  
  outdata << std::setw(6) << filePointsSize << std::setw(5);
  
  for (int i = 0; i < scores.size(); i++)
    outdata  << "  || " << std::setw(2) << std::fixed << std::setprecision(6) << scores[i].first << "  | " << scores[i].second << "  | " << minScores[i] << "  | " << maxScores[i];

  outdata << "  || " << std::endl;

  outdata.close();
  
}