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


// Check if we need to exevute polygomization algorithm or scoring tables
std::string dataio::getRunType(int argc, char** argv){

  for(int i = 1; i < argc; i++){

    std::string param = argv[i];

    if(param == "-scores")
      return "SCORES";

  }

  return "EXECUTION";

}


// Get execution type parameters
bool dataio::getExecutionParameters(std::string& nameOfFile, std::string& outputFile, std::string& algorithm, std::string& algorithm_initial,std::string& initial, int& edge_selection_int, int& polygon_edge_selection, double& threshold,std::string& annealing,int& L, int argc, char** argv, int& m) {
  
  polygon_edge_selection  = -1;
  m = -1;
  threshold = -1.0;
  int max_min = 0;
  std::string edge_selection;

  for(int i = 0; i < argc- 1; i++){

    std::string param = (argv[i]);
    std::string mode = (argv[i+1]);

    if(param == "-i"){
      nameOfFile = mode;
    }
    else if(param == "-o"){
      outputFile = mode;
    }
    else if(param == "-algorithm_initial"){
      algorithm_initial = mode;
    }
    else if(param == "-algorithm"){
      algorithm = mode;
    }
    else if(param == "-max"){
      edge_selection = "-max";
      if(max_min == 1)
        return false;
      max_min++;  
    }
    else if(param == "-min"){
      edge_selection = "-min";
      if(max_min == 1)
       return false;
      max_min++;
    }
    else if(param == "-L"){
      L = stod(mode);
    }
    else if(param == "-threshold"){
      threshold = stod(mode);
    }
    else if(param == "-annealing"){
      annealing = mode;
    }
    else if(param == "-initialization"){
      initial = mode;
    }
    else if(param == "-m"){
      m = stoi(mode);
    }
    else if(param == "-edge_selection"){
      polygon_edge_selection = stoi(mode);
    }
  }
  if(m < 0){
    m = 10;
  }
  if ((algorithm_initial != "incremental") && (algorithm_initial != "convex_hull")){
    return false;
  }

  if (max_min == 0){ 
    return false;
  }

  if(edge_selection == "-max"){
    edge_selection_int = 3;
  }

  else if(edge_selection == "-min"){
    edge_selection_int = 2;
  }

  if(polygon_edge_selection < 0 || polygon_edge_selection > 3)
    polygon_edge_selection = edge_selection_int;

  if ((algorithm == "local_search") && (threshold < 0.0))
    return false;

  if ((algorithm == "simulated_annealing") && (annealing != "local") && (annealing != "global") && (annealing != "subdivision"))
    return false;  

  if ((algorithm_initial == "incremental") && (initial != "1a") && (initial != "1b") && (initial != "2a") && (initial != "2b"))
    return false;

  if ((algorithm != "local_search") && (algorithm != "simulated_annealing"))
    return false;

  return true;
}



// Get scoring type parameters
bool dataio::getScoreParameters(std::string& nameOfDirectory, std::string& outputFile, bool& preprocessEnabled, int argc, char** argv) {
  
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
void dataio::createResultsFile(const std::vector<Segment_2> &polygLine, const ft& area, const ft& areaBefore, const ft& ratio,const ft& ratioBefore, const std::chrono::milliseconds& polygonizationDuration, const std::string& output, const std::string& algorithm, const int& edgeSelection, const std::string& initialization, const std::string& max_min) {

  std::ofstream outdata;

  outdata.open(output.c_str());

  if (!outdata) {
    std::cout << "Error: file could not be opened" << std::endl;
    exit(EXIT_FAILURE);
  }

  outdata << "Optimal Area Polygonization" << std::endl;

  for (const Segment_2& line : polygLine)
    outdata << line.source() << std::endl;

  outdata << std::endl;

  for (const Segment_2& line : polygLine)
    outdata << line << std::endl;

  outdata << std::endl;

  outdata << "Algorithm: "<< algorithm << "_" << max_min << std::endl;

  outdata << "Area_initial: " << (long int)areaBefore << std::endl;

  outdata << "Area: " << (long int)area << std::endl;

  outdata << "ratio_initial: " << ratioBefore << std::endl;

  outdata << "ratio: " << ratio << std::endl;

  outdata << "Construction time: " << polygonizationDuration.count() << std::endl;

  outdata.close();
  
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

  if (initOutput){
    outdata  << std::setw(10) << "||" << std::setw(29) << "Incr+Global+Local" <<  std::setw(20) << "||" << std::setw(29) << "Subdivision" << std::setw(20) << "||" << std::setw(29) << "Incr+Local" << std::setw(20) << "||" << std::setw(29) << "Convex+Local" << std::setw(20) << "||" << std::endl;
  
    outdata  << std::setw(8) << "Size " << "||" << std::setw(2) << " min score" << " | " << "max score" << " | " << "min bound" << " | " << "max bound";
   
    for (int i = 0; i < 3; i++){
      outdata << " ||" << std::setw(2) << " min score" << " | " << "max score" << " | " << "min bound" << " | " << "max bound";
    }

   outdata << " ||" << std::endl;

  }
  
  outdata << std::setw(6) << filePointsSize << std::setw(5);
  
  for (int i = 0; i < scores.size(); i++)
    outdata  << "  || " << std::setw(2) << std::fixed << std::setprecision(6) << scores[i].first << "  | " << scores[i].second << "  | " << minScores[i] << "  | " << maxScores[i];

  outdata << "  || " << std::endl;

  outdata.close();
  
}
