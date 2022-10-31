#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <string>
#include <CGAL/squared_distance_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point>::type> Convex_hull_traits_2;
typedef CGAL::Epick::FT ft;

Segment_2 findVisibleSegment(std::vector<Segment_2> &polygLine, Segment_2 &convexSegment, Point &nextPoint) {
      
  // Convex hull segment == polygon line segment and not collinear check
  if ((std::find(polygLine.begin(), polygLine.end(), convexSegment) != polygLine.end())){ 
    
    if (!CGAL::collinear(nextPoint, convexSegment.point(0), convexSegment.point(1)))
      return convexSegment;
    else 
      throw std::runtime_error("No Visible Segment Available");
  
  }

  // Convex hull segment not on polygon line segment
  Segment_2 segment = convexSegment;

  // std::cout << "convex segment is " << segment << std::endl;
  int index = 0;
  
  for (int i = 0; i < polygLine.size(); i++){
    if (polygLine[i].point(0) == segment.point(0)){
      index = i;
      break;
    }
  }
  
  std::vector<Segment_2> visibleSegments;

  for (int i = index; i < polygLine.size(); i++){
    
    if (polygLine[i].point(1) != segment.point(1)){
      visibleSegments.push_back(polygLine[i]);
    }
    else{
      visibleSegments.push_back(polygLine[i]);
      break;
    }
  }

  auto it = visibleSegments.begin();

  while(it != visibleSegments.end()) {
    // std::cout << "Visible ==> " << vSeg << std::endl;
    auto indexing = it;
    indexing++;

    if (CGAL::collinear(nextPoint, (*it).point(0), (*it).point(1))) {
      int prevSize = visibleSegments.size();
      visibleSegments.erase(it);
      if (visibleSegments.size() != prevSize)
        indexing--;
    }

    it = indexing;
  }

  auto at = visibleSegments.begin();

  while (at != visibleSegments.end()) {

    auto indexing = at;
    indexing++;
    
    Segment_2 segmentsArray[] = {Segment_2(nextPoint, (*at).point(0)),
                                 Segment_2(nextPoint, (*at).point(1)),
                                 Segment_2(nextPoint, Point((((*at).point(0).x() + (*at).point(1).x()) / 2), (((*at).point(0).y() + (*at).point(1).y()) /2)))
                                };

    bool flag = false;
    for (Segment_2 polygSeg : polygLine){
      
      if ((*at) == polygSeg) continue;

      for (int k = 0; k < 3; k++) {
        const auto result = intersection(segmentsArray[k], polygSeg);

        if (result) {
          
          if (const Point *p = boost::get<Point>(&*result)) {
            if (*p == polygSeg.point(0) || *p == polygSeg.point(1)){
              continue;
            }
            else {
              visibleSegments.erase(at);

              indexing--;
              
              flag = true;
              break;
            }
          }
          else {
            visibleSegments.erase(at);
           
            indexing--;

            flag = true;
            break;
          }
        }

      }
      if (flag)
        break;
    }

    at = indexing;
  }

  srand(time(NULL));
  
  int randomIndex = visibleSegments.size() ? rand() % (visibleSegments.size()) : throw std::runtime_error("No Visible Segment Available");

  return visibleSegments[randomIndex];
}


void createResultsFile(std::vector<Segment_2> &polygLine){
  
  std::ofstream outdata;              
  
  int i;                         

  outdata.open("results.txt"); 
  
  if (!outdata) {               
    std::cerr << "Error: file could not be opened" << std::endl;
    exit(1);
  }

  for (i = 0; i < polygLine.size(); ++i)
    outdata << polygLine[i] << std::endl;
  
  outdata.close();

}

std::vector<Point> readPoints() {

  std::ifstream inFile("data.in");
  std::string strInput;

  std::vector<Point> points;
  
  int currentIndex= 0;

  while (inFile){ 

    getline(inFile, strInput);
    
    if (strInput[0] == '#' || strInput.size() < 2)
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

// Get all points from a vector of segments
std::vector<Point> getPolyLinePoints(std::vector<Segment_2> &polygLine) {

  std::vector<Point> linePoints;

  for (int i = 0; i < polygLine.size(); i++) {
    linePoints.push_back(polygLine[i].point(0));
  }

  return linePoints;
}

// std::vector<Segment_2> getPolyLineSegs(std::vector<Point> &polygLinePoints) {
//   Polygon_2 Polygon = Polygon_2();

//   for (auto it = polygLinePoints.begin(); it != polygLinePoints.end(); ++it) {
//     Polygon.push_back(*it);
//   }

//   std::vector<Segment_2> polSegments;

//   // Create A Vector That Stores Convex Hull segments
//   for (Segment_2 seg : Polygon.edges()) {
//     polSegments.push_back(seg);
//   }

//   return polSegments;
// }

std::vector<Segment_2> getConvexHull(std::vector<Point> &polygLinePoints,std::vector<Point> &remainingPoints) {

  std::vector<Point> convexHullPoints;

  // Find Complex Hull Points
  CGAL::convex_hull_2(polygLinePoints.begin(), polygLinePoints.end(), std::back_inserter(convexHullPoints));

  Polygon_2 convexHullPolygon = Polygon_2();

  for (auto it = convexHullPoints.begin(); it != convexHullPoints.end(); ++it) {
    convexHullPolygon.push_back(*it);
    remainingPoints.erase(std::remove(remainingPoints.begin(), remainingPoints.end(), *it), remainingPoints.end());
  }

  std::vector<Segment_2> convexHullSegments;

  // Create A Vector That Stores Convex Hull segments
  for (Segment_2 seg : convexHullPolygon.edges()) {
    convexHullSegments.push_back(seg);
  }

  return convexHullSegments;
}

void initializeConvexHull(std::vector<Segment_2> &polygLine, std::vector<Point> &points, std::vector<Point> &remainingPoints) {

  std::vector<Point> convexPoints;

  for (int i = 0; i < points.size(); i++)
    convexPoints.push_back(points[i]);

  polygLine = getConvexHull(convexPoints, remainingPoints);

}

void deleteSegment(std::vector<Segment_2> &polygLine, Segment_2 &visibleSegment){
   polygLine.erase(std::remove(polygLine.begin(), polygLine.end(), visibleSegment), polygLine.end());
}


void expandPolygonLine(std::vector<Segment_2> &polygLine, Segment_2 &visibleSegment, Point &nextPoint){
  // Insert the two new segments in the right place in polygon line
  int index = 0;
  for (int i = 0; i < polygLine.size(); i++) {
    if (polygLine[i].point(1) == visibleSegment.point(0)){
      
      polygLine.insert(polygLine.begin() + index + 1,Segment_2(visibleSegment.point(0), nextPoint));
      polygLine.insert(polygLine.begin() + index + 2, Segment_2(nextPoint, visibleSegment.point(1)));
      break;
      
    } 
    else if (polygLine[i].point(1) == visibleSegment.point(1)) {
      polygLine.insert(polygLine.begin() + index + 1,Segment_2(visibleSegment.point(1), nextPoint));
      polygLine.insert(polygLine.begin() + index + 2, Segment_2(nextPoint, visibleSegment.point(0)));
      break;
    }
    index++;
  }
}

typedef struct visible_point {
  Point cor;
  double distance;
} visPoint;

typedef struct pairPointSeg {
  Point cor;
  Segment_2 seg;
} pair;

void initialRun(std::vector<Segment_2>& currConvexHullSegments,std::vector<Point>& remainingPoints,std::vector<Segment_2>& polygLine){
    bool done = false;
    for(int i = 0; i < currConvexHullSegments.size(); i++){
      std::vector<visPoint> visPoints;
      if(done){
        printf("done\n");
        return;
      }
      for(int j = 0; j < remainingPoints.size(); j++){
        Segment_2 segmentsArray[] = {Segment_2(remainingPoints[j], currConvexHullSegments[i].point(0)),
                                 Segment_2(remainingPoints[j], currConvexHullSegments[i].point(1)),
                                 Segment_2(remainingPoints[j], Point((currConvexHullSegments[i].point(0).x() + currConvexHullSegments[i].point(1).x()) / 2, ((currConvexHullSegments[i].point(0).y() + currConvexHullSegments[i].point(1).y()) /2)))
                                };
        bool flag = false;
        for (Segment_2 polygSeg : polygLine){
          if ( currConvexHullSegments[i] == polygSeg) continue;

          for (int k = 0; k < 3; k++) {
            const auto result = intersection(segmentsArray[k], polygSeg);

            if (result) {
          
              if (const Point *p = boost::get<Point>(&*result)) {
                if (*p == currConvexHullSegments[i].point(0) || *p == currConvexHullSegments[i].point(1)){
                  continue;
                }
                else {
                  flag = true;
                  break;
                }
            }else {
            flag = true;
            break;
          }
        }

      }
      if (flag){
        break;
      }
      
    }
    if(!flag) {
      printf("inserted!\n");
      double distance =
          squared_distance(remainingPoints[j], currConvexHullSegments[i]);
      Point insPoint = remainingPoints[j];
      visPoint inserted;
      inserted.cor = insPoint;
      inserted.distance = distance;
      visPoints.push_back(inserted);
      done = true;
    }
  }
  if(!visPoints.size())
    continue;
  int index = 0;
  Point bestPoint = visPoints[0].cor;
  double bestDist = visPoints[0].distance;
  for(int k = 0;k < visPoints.size();k++){
        if(bestDist > visPoints[i].distance){
          bestPoint = visPoints[i].cor;
          bestDist = visPoints[i].distance;
        }
    }
    for(int m = 0; m < remainingPoints.size(); m++){
      if(remainingPoints[m] == bestPoint){
        remainingPoints.erase(remainingPoints.begin()+m);
        Segment_2 visSeg = currConvexHullSegments[i];
        deleteSegment(polygLine, visSeg);
        expandPolygonLine(polygLine, visSeg,bestPoint);
        // printf("poped!%d\n",remainingPoints.size());
        break;
      }
    }
  }
}

bool forceInsertPoint(std::vector<Segment_2> &polygLine, Point &nextPoint) {

  bool forceInserted = false;

  for (int i = 0; i < polygLine.size(); i++) {
    if (polygLine[i].has_on(nextPoint)) {
      polygLine.insert(polygLine.begin() + i + 1, Segment_2(polygLine[i].source(), nextPoint));
      polygLine.insert(polygLine.begin() + i + 2, Segment_2(nextPoint, polygLine[i].target()));
      polygLine.erase(polygLine.begin() + i);
      forceInserted = true;
      break;
    }
  }

  if (forceInserted)
    return true;
  
  return false;
}

int polCheck(std::vector<Segment_2> polygLine){
  Polygon_2 pol_result = Polygon_2();

  for (Segment_2 segment : polygLine) {
    pol_result.push_back(segment.point(0));
  }

  std::cout << "Polygon Is Simple: " << pol_result.is_simple() << std::endl;
  if(!pol_result.is_simple()){
    createResultsFile(polygLine);
    printf("ending...\n");
    exit(1);
  }
  printf("all ok!\n");
  return 1;
}

void findVisiblePoints(std::vector<visPoint>& visPoints,Point remainingPoint,Segment_2 seg,std::vector<Segment_2> polygLine){
    Segment_2 segmentsArray[] = {Segment_2(remainingPoint, seg.point(0)),
                                 Segment_2(remainingPoint, seg.point(1)),
                                 Segment_2(remainingPoint, Point((seg.point(0).x() + seg.point(1).x()) / 2, ((seg.point(0).y() + seg.point(1).y()) /2)))
                                };
        bool flag = false;
        for (Segment_2 polygSeg : polygLine){
          
          for (int k = 0; k < 3; k++) {
            const auto result = intersection(segmentsArray[k], polygSeg);

            if (result) {
              if (const Point *p = boost::get<Point>(&*result)) {
                if (*p ==  seg.point(0) || *p ==  seg.point(1) || *p == Point((seg.point(0).x() + seg.point(1).x()) / 2, ((seg.point(0).y() + seg.point(1).y()) /2))){
                  continue;
                }
                else {
                  flag = true;
                  break;
                }
              } 
              else {
                flag = true;
                break;
              }
            }
          }

          if (flag){
            break;
          }
          
        }

        if (!flag) {
          double distance = squared_distance(remainingPoint, seg);
          Point insPoint = remainingPoint;
          visPoint inserted;
          inserted.cor = insPoint;
          inserted.distance = distance;
          visPoints.push_back(inserted);
        }
}

Point findBestPoint(std::vector<visPoint> visPoints,std::vector<Point> remainingPoints,std::vector<Segment_2>& polygLine,Segment_2 visSeg){
    int index = 0;
    Point bestPoint = visPoints[0].cor;
    double bestDist = visPoints[0].distance;

    for(int k = 0; k < visPoints.size(); k++){
        if(bestDist > visPoints[k].distance){
          bestPoint = visPoints[k].cor;
          bestDist = visPoints[k].distance;
        }
    }
    return bestPoint;
    // for(int m = 0; m < remainingPoints.size(); m++){
    //     if(remainingPoints[m] == bestPoint){

    //       remainingPoints.erase(remainingPoints.begin()+m);
          
    //       deleteSegment(polygLine, visSeg);
    //       expandPolygonLine(polygLine, visSeg, bestPoint);
    //       polCheck(polygLine);
    //       break;
    //     }
    // }
}

void insertBestPoint(std::vector<pair> bestPoints,std::vector<Point>& remainingPoints,std::vector<Segment_2>& polygLine,int mode){
    if(mode != 1){
        std::vector<Segment_2> testPolyg = polygLine;
        deleteSegment(testPolyg, bestPoints[0].seg);
        expandPolygonLine(testPolyg, bestPoints[0].seg, bestPoints[0].cor);
        std::vector<Point> polygLinePoints = getPolyLinePoints(testPolyg);
        ft chosenArea = CGAL::polygon_area_2(polygLinePoints.begin(), polygLinePoints.end(), Convex_hull_traits_2(CGAL::make_property_map(polygLinePoints)));
        int index = 0;
        pair bestPair;
        bestPair.cor = bestPoints[0].cor;
        bestPair.seg = bestPoints[0].seg;
        for(int i = 0; i < bestPoints.size(); i++){
            std::vector<Segment_2> testPolyg = polygLine;
            deleteSegment(testPolyg, bestPoints[i].seg);
            expandPolygonLine(testPolyg, bestPoints[i].seg, bestPoints[i].cor);
            std::vector<Point> polygLinePoints = getPolyLinePoints(testPolyg);
            ft polArea = CGAL::polygon_area_2(polygLinePoints.begin(), polygLinePoints.end(), Convex_hull_traits_2(CGAL::make_property_map(polygLinePoints)));
            if(mode == 2 && polArea > chosenArea){
                index = i;
                chosenArea = polArea;
                bestPair.cor = bestPoints[i].cor;
                bestPair.seg = bestPoints[i].seg;
            }
            else if(mode == 3 && polArea < chosenArea){
                index = i;
                chosenArea = polArea;
                bestPair.cor = bestPoints[i].cor;
                bestPair.seg = bestPoints[i].seg;
            }
        }
        for(int m = 0; m < remainingPoints.size(); m++){
            if(remainingPoints[m] == bestPair.cor){

                remainingPoints.erase(remainingPoints.begin()+m);
                printf("erased!\n");
                deleteSegment(polygLine, bestPair.seg);
                expandPolygonLine(polygLine, bestPair.seg, bestPair.cor);
                polCheck(polygLine);
                break;
            }
        }
    }
    else{
        printf("at least here!\n");
        int rindex = ( rand() % (bestPoints.size() - 1));
        pair bestPair;
        printf("RAND:%d\n",bestPoints.size());
        bestPair.cor = bestPoints[rindex].cor;
        bestPair.seg = bestPoints[rindex].seg;
        for(int m = 0; m < remainingPoints.size(); m++){
            if(remainingPoints[m] == bestPair.cor){
                remainingPoints.erase(remainingPoints.begin()+m);
                deleteSegment(polygLine, bestPair.seg);
                expandPolygonLine(polygLine, bestPair.seg, bestPair.cor);
                polCheck(polygLine);
                printf("here!\n");
                break;
            }
        }
    }
}

int main(){
    srand(time(NULL));
  // Read Points
  std::vector<Point> points = readPoints();

  std::vector<Point> remainingPoints = points;

  std::vector<Segment_2> polygLine;

  initializeConvexHull(polygLine, points, remainingPoints);

  std::cout << "Size of remainig " << remainingPoints.size() << std::endl;

  for(int i = 0; i < remainingPoints.size(); i++)
    std::cout << "Point " << remainingPoints[i] << std::endl;

  std::vector<Point> polygLinePoints = getPolyLinePoints(polygLine);

  std::vector<Segment_2> currConvexHullSegments = getConvexHull(polygLinePoints,remainingPoints);
  initialRun(currConvexHullSegments, remainingPoints, polygLine);
  while(remainingPoints.size() > 0){
    std::vector<pair> bestPoints;
    for(int i = 0; i < polygLine.size(); i++){
      std::vector<visPoint> visPoints;
      for(int j = 0; j < remainingPoints.size(); j++){
        findVisiblePoints(visPoints,remainingPoints[j],polygLine[i],polygLine);
      }
      
      if(visPoints.size() == 0){
    
        // for (int i = 0; remainingPoints.size(); i++){
        //   if (forceInsertPoint(polygLine, remainingPoints[i])) {
        //     remainingPoints.erase(remainingPoints.begin() + i);
        //     i--;
        //   }

        // }
        continue;
      }

        Point bestPoint = findBestPoint(visPoints,remainingPoints,polygLine,polygLine[i]);
        Segment_2 bestSeg = polygLine[i];
        pair best = {bestPoint,bestSeg};
        bestPoints.push_back(best);
    }
    insertBestPoint(bestPoints,remainingPoints,polygLine,1);
  }

  Polygon_2 pol_result = Polygon_2();

  for (Segment_2 segment : polygLine) {
    pol_result.push_back(segment.point(0));
  }

  std::cout << "Polygon Is Simple: " << pol_result.is_simple() << std::endl;
  createResultsFile(polygLine);
  return 0;
}
