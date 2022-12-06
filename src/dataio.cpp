#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <string>
#include <vector>
#include "../include/dataio.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;


// Get parameters
bool dataio::getParameters(std::string& nameOfFile, std::string& outputFile, std::string& algorithm, std::string& algorithm_initial,std::string& initial, int& edge_selection_int,double& threshold,int& annealing_int,int& L, int argc, char** argv) {
  int max_min = 0;
  std::string edge_selection;
  std::string annealing;
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
      edge_selection = mode;
      if(max_min == 1)
        return false;
      max_min++;  
    }
    else if(param == "-min"){
      edge_selection = mode;
      if(max_min == 1)
       return false;
      max_min++;
    }
    else if(param == "-L"){
      L = stoi(mode);
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
  }

  if ((algorithm_initial != "incremental") && (algorithm_initial != "convex_hull")){
    printf("here\n");
    return false;
  }

  if (max_min == 0){ 
    printf("edw!\n");
    return false;
  }
  if(annealing == "local"){
    annealing_int = 1;
  }
  else if(annealing == "global"){
    annealing_int = 2;
  }
  else if(annealing == "simulated"){
    annealing_int = 0;
  }
  else{
    printf("hah got em!\n");
    return false;
  }
  if(edge_selection == "-max"){
    edge_selection_int = 3;
  }
  else if(edge_selection == "-min"){
    edge_selection_int = 2;
  }

  if ((algorithm_initial == "incremental") && (initial != "1a") && (initial != "1b") && (initial != "2a") && (initial != "2b"))
    return false;

  if ((algorithm != "local_search") && (algorithm != "simulated_annealing"))
    return false;

  return true;
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
void dataio::createResultsFile(const std::vector<Segment_2> &polygLine, const ft& area, const std::chrono::milliseconds& polygonizationDuration, const ft& ratio, const std::string& output, const std::string& algorithm, const int& edgeSelection, const std::string& initialization) {

  std::ofstream outdata;

  outdata.open(output.c_str());

  if (!outdata) {
    std::cout << "Error: file could not be opened" << std::endl;
    exit(1);
  }

  outdata << "Polygonization" << std::endl;

  for (const Segment_2& line : polygLine)
    outdata << line.source() << std::endl;

  for (const Segment_2& line : polygLine)
    outdata << line << std::endl;

  outdata << std::endl;

  outdata << "Algorithm: "<< algorithm << "_edge_selection_" << edgeSelection;

  if (algorithm == "incremental")
    outdata << "_initialization_" << initialization << std::endl;
  else
   outdata << std::endl;

  outdata << "Area: " << (long int)area << std::endl;

  outdata << "ratio: " << ratio << std::endl;

  outdata << "Construction time: " << polygonizationDuration.count() << std::endl;

  outdata.close();
  
}