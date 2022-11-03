#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <string>
#include <vector>
#include "../include/dataio.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;

// Read from data.instance file
std::vector<Point> dataio::readPoints(std::string name) {

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
    std::cout << points.back().x() << " " << points.back().y() << std::endl;
  }

  inFile.close();

  return points;
}


// Write resulta data into a spesific file
void dataio::createResultsFile(const std::vector<Segment_2> &polygLine,const ft& area, const std::chrono::milliseconds& polygonizationDuration, const ft& ratio,std::string output) {

  std::ofstream outdata;

  int i;

  outdata.open(output.c_str());

  if (!outdata) {
    std::cerr << "Error: file could not be opened" << std::endl;
    exit(1);
  }

  for (i = 0; i < polygLine.size(); ++i)
    outdata << polygLine[i] << std::endl;

  outdata.close();

  outdata.open("Polygon_Results.txt");

  outdata << "Polygonization" << std::endl;

  for (Segment_2 line : polygLine)
    outdata << line.source() << std::endl;

  for (Segment_2 line : polygLine)
    outdata << line << std::endl;

  outdata << std::endl;

  outdata << "Algorithm: "<< "incremental" << "_edge_selection_" << "1" "_initialization_" << "1a" << std::endl;

  outdata << "Area: " << (int)area << std::endl;

  outdata << "ratio: " << ratio << std::endl;

  outdata << "Construction time: " << polygonizationDuration.count() << std::endl;

  outdata.close();
  
}