#include "../include/dataio.hpp"
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/enum.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>
#include "../include/polygonization.hpp"
#include "../include/simulatedAnnealing.hpp"
#include "../include/incremental.hpp"
#include "../include/convexHull.hpp"



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point>::type> Convex_hull_traits_2;
typedef CGAL::Epick::FT ft;



void simulatedAnnealing::mergePolygons(std::vector<std::vector<Point>>& subPolPoints, std::vector<std::vector<Segment_2>>& allPolygLines){

  for (int i = 0; i < allPolygLines.size(); i++){
    
    int indexStart;
    int indexEnd;

    for (int j = 0; j < allPolygLines[i].size(); j++){
      
      if (allPolygLines[i][j].source() == subPolPoints[i][0]){
        indexStart = j + 1 > allPolygLines[i].size() - 1 ? 0 : j + 1;
      }
      else if (allPolygLines[i][j].target() == subPolPoints[i][subPolPoints[i].size() - 2]){
         indexEnd = j;
      }
    }
    
     if (i != 0)
      polygLine.push_back(Segment_2(polygLine[polygLine.size() - 1].target(), allPolygLines[i][indexStart].source()));
    
    
    while (indexStart != indexEnd){
       polygLine.push_back(allPolygLines[i][indexStart]);
       std::cout << "Pushed back " << allPolygLines[i][indexStart] << std::endl;
       indexStart = indexStart + 1 > allPolygLines[i].size() - 1 ? 0 : indexStart + 1;
    }

    

    if (i < allPolygLines.size() - 1){
      polygLine.push_back(allPolygLines[i][indexStart]);
      std::cout << "Pushed back " << allPolygLines[i][indexStart] << std::endl;
    }

  }

  for (int i = allPolygLines.size() - 1; i >= 0; i--){
    
    int indexStart;
    int indexEnd;

    for (int j = 0; j < allPolygLines[i].size(); j++){

      if (allPolygLines[i][j].source() == subPolPoints[i][subPolPoints[i].size() - 1]){
        indexStart = j;
      }
      else if (allPolygLines[i][j].target() == subPolPoints[i][0]){
         indexEnd = j;
      }
    }

    
    while (indexStart != indexEnd){
       polygLine.push_back(allPolygLines[i][indexStart]);
       indexStart = indexStart + 1 > allPolygLines[i].size() - 1 ? 0 : indexStart + 1;
    }
    polygLine.push_back(allPolygLines[i][indexStart]);

     if (i == 0){
      polygLine.push_back(Segment_2(allPolygLines[i][indexEnd].target(), polygLine[0].source()));
    } 
  }

}


void simulatedAnnealing::createSubsetPoints(std::vector<std::vector<Point>>& subPolPoints){

  int m = 15, pIndex = 0, threshold = 5, n = points.size(), k = std::ceil((double)(n - 1)/(double)(m - 1));


  // Create m-teams
  for (int i = 0; i < k; i++){

    if (pIndex == n) break;

    std::vector<Point> vecInit;
    vecInit.clear();
    subPolPoints.push_back(vecInit);

    for (int j = 0; j < m; j++){
      subPolPoints[i].push_back(points[pIndex]);
      pIndex = j == m - 1 ? pIndex : pIndex + 1;
      if (pIndex == n) break;
    }

    int lastIndex = subPolPoints[i].size() - 1;
    
    while (pIndex < n - 1 && !(subPolPoints[i][lastIndex - 1].y() < subPolPoints[i][lastIndex].y() && subPolPoints[i][lastIndex].y() > points[pIndex + 1].y())){
      pIndex++;
      subPolPoints[i].push_back(points[pIndex]);
      lastIndex = subPolPoints[i].size() - 1;
    }

  }

  // Delete Last Team if < threashold
  int lastPointsSetIndex = subPolPoints.size() - 1;
  if (subPolPoints[lastPointsSetIndex].size() < threshold){
    for (int i = 1; i < subPolPoints[lastPointsSetIndex].size(); i++){
      subPolPoints[lastPointsSetIndex - 1].push_back(subPolPoints[lastPointsSetIndex][i]);
    }
    subPolPoints.erase(subPolPoints.begin() + lastPointsSetIndex);
    k--;
  }


  // Prints points subsets
  for (int i = 0; i < subPolPoints.size(); i++){
    for (int j = 0; j < subPolPoints[i].size(); j++){
      std::cout << "p --> " << subPolPoints[i][j] << std::endl;
    }
    std::cout << std::endl;
  }

}

void simulatedAnnealing::subPolygonization(std::vector<std::vector<Point>>& subPolPoints, std::vector<std::vector<Segment_2>>& allPolygLines, int edge_selection){

   for (int i = 0; i < subPolPoints.size(); i++){
    // convexHull pol = convexHull(subPolPoints[i], edge_selection,Segment_2(subPolPoints[i][0],subPolPoints[i][1]),Segment_2(subPolPoints[i][subPolPoints[i].size() - 2],subPolPoints[i][subPolPoints[i].size() - 1]), i, subPolPoints.size() - 1);
    // pol.start();
    // std::vector<Segment_2> vec = pol.getPolygonLine();
    // allPolygLines[i] = vec;

    Incremental pol = Incremental(subPolPoints[i], edge_selection, "1a", subPolPoints[i][1]);
    pol.start();
    std::vector<Segment_2> vec = pol.getPolygonLine();
    allPolygLines[i] = vec;

    for (const Segment_2& segment : pol.getPolygonLine()){
      std::cout << "Segment --> " << segment << std::endl;
    }

    std::cout << "-------------------------------------------" << std::endl;
  }

}


void simulatedAnnealing::subGlobalTransitions(std::vector<std::vector<Point>>& subPolPoints, std::vector<std::vector<Segment_2>>& allPolygLines){
  
  for (int i = 0; i < allPolygLines.size(); i++){
    
    std::vector<Segment_2> untouchable;

    if (i != 0)
      untouchable.push_back(Segment_2(subPolPoints[i][0],subPolPoints[i][1]));
    
    if (i != allPolygLines.size() - 1)
      untouchable.push_back(Segment_2(subPolPoints[i][subPolPoints[i].size() - 2],subPolPoints[i][subPolPoints[i].size() - 1]));

    simulatedAnnealing annealing = simulatedAnnealing(points , allPolygLines[i], calcArea(allPolygLines[i]), calcRatio(allPolygLines[i], calcArea(allPolygLines[i])), 1000, 2, 2, untouchable);
    annealing.startAnnealing();

    allPolygLines[i] = annealing.getPolygonLine();

  }

}

void simulatedAnnealing::startSubdivision(int edge_selection, const std::string& initialization){

  // Sorting lex order
  std::sort(points.begin(), points.end(), lexOrderPoints);

  // Create subset teams
  std::vector<std::vector<Point>> subPolPoints;
  createSubsetPoints(subPolPoints);

  // Polygonization for every subset
  std::vector<std::vector<Segment_2>> allPolygLines(subPolPoints.size());
  subPolygonization(subPolPoints, allPolygLines, edge_selection);

  // Global transitions for every sub polygon
  subGlobalTransitions(subPolPoints, allPolygLines);

  // Merge all subset polygons
  mergePolygons(subPolPoints, allPolygLines);

  // Local transitions for final polygon
  simulatedAnnealing annealing = simulatedAnnealing(points , polygLine, calcArea(polygLine), calcRatio(polygLine, calcArea(polygLine)), 1000, 2, 1);
  annealing.startAnnealing();

  polygLine = annealing.getPolygonLine();

  for (const Segment_2& segment : polygLine){
      std::cout << "CONSTRUCTED ==> " << segment << std::endl;
  }

  std::cout << "Polygon Simplicity: " << checkPolygonSimplicity(polygLine) << std::endl;
}


