#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <chrono>
#include <iostream>
#include <vector>
#include "../include/polygonization.hpp"
#include "../include/convexHull.hpp"
#include "../include/dataio.hpp"
#include "../include/incremental.hpp"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;



int main(int argc,char** argv) {
  std::string nameOfFile;
  std::string outputFile;
  std::string algorithm;
  int edge_selection;
  std::string initialization;
  for(int i = 0; i < argc; i++){
    std::string temp(argv[i+1]);
    if(!strcmp("-i",argv[i])){
      nameOfFile = temp;
    }
    else if(!strcmp("-o",argv[i])){
      outputFile = temp;
    }
    else if(!strcmp("-algorithm",argv[i])){
      algorithm = temp;
    }
    else if(!strcmp("-edge_selection",argv[i])){
      edge_selection = std::stoi(temp);
    }
    else if(!strcmp("initialization",argv[i])){
      initialization = temp;
    }
  }
  
  std::vector<Point> points = dataio::readPoints(nameOfFile);
  if(algorithm == "incremental"){
    Incremental pol = Incremental(points,edge_selection,initialization);

    auto start = std::chrono::high_resolution_clock::now();

    pol.start();

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    dataio::createResultsFile(pol.getPolygonLine(), pol.getArea(), duration, pol.getRatio(),outputFile);
  }
  else if(algorithm == "convex_hull"){
    convexHull pol = convexHull(points,edge_selection);

    auto start = std::chrono::high_resolution_clock::now();

    pol.start();
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    dataio::createResultsFile(pol.getPolygonLine(), pol.getArea(), duration, pol.getRatio(),outputFile);
  }
  return 0;
}