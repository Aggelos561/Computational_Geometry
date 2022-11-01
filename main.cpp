#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point>::type> Convex_hull_traits_2;
typedef CGAL::Epick::FT ft;


// Read from data.instance file
std::vector<Point> readPoints() {

  std::ifstream inFile("data.in");
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

void createResultsFile(std::vector<Segment_2> &polygLine, ft& area, std::chrono::milliseconds& polygonizationDuration, ft& ratio) {

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



// Sort Y Ascending
bool sortYAsc(Point& a, Point& b){
  if (a.y() < b.y())
    return true;
  
  return false;
}

//Sort Y descending
bool sortYDesc(Point &a, Point &b) {
  if (a.y() > b.y())
    return true;

  return false;
}

// Sort points in vector
void sortPoints(std::vector<Point> &points) {
  std::sort(points.begin(), points.end(), sortYDesc);
}



// Get all points from a vector of segments
std::vector<Point> getPolyLinePoints(std::vector<Segment_2> &polygLine) {

  std::vector<Point> linePoints;

  for (int i = 0; i < polygLine.size(); i++) {
    linePoints.push_back(polygLine[i].point(0));
  }

  return linePoints;
}


std::vector<Segment_2> getConvexHull(std::vector<Point> &polygLinePoints) {

  std::vector<Point> convexHullPoints;

  // Find Complex Hull Points and rerurn a vector of segments
  CGAL::convex_hull_points_2(polygLinePoints.begin(), polygLinePoints.end(), std::back_inserter(convexHullPoints));

  Polygon_2 convexHullPolygon = Polygon_2();

  for (Point p : convexHullPoints) {
    convexHullPolygon.push_back(p);
  }

  std::vector<Segment_2> convexHullSegments;

  // Create A Vector That Stores Convex Hull segments
  for (Segment_2 seg : convexHullPolygon.edges()) {
    convexHullSegments.push_back(seg);
  }

  return convexHullSegments;
}


// Initialize triangle 
void initializeTriangle(std::vector<Segment_2> &polygLine, std::vector<Point> &points, std::vector<Point> &remainingPoints) {

  std::vector<Point> trianglePoints;

  for (int i = 0; i < 3; i++){
    trianglePoints.push_back(points[i]);
  }

  int index = 3;
  polygLine = getConvexHull(trianglePoints);
  while (polygLine.size() != 3){
    trianglePoints.erase(trianglePoints.end() - 1);
    trianglePoints.push_back(points[index]);
    index++;
    
    polygLine = getConvexHull(trianglePoints);
  }

  for (int i = 0; i < polygLine.size(); i++){

    for (int j = 0; j < remainingPoints.size(); j++){
      if (polygLine[i].source() == remainingPoints[j]){
        remainingPoints.erase(remainingPoints.begin() + j);
        j--;
      }
    }

  }

}


//Deletes a spesific segment value from a vector of segments
void deleteSegment(std::vector<Segment_2> &polygLine, Segment_2 &visibleSegment){
  polygLine.erase(std::remove(polygLine.begin(), polygLine.end(), visibleSegment), polygLine.end());
}


// Get triangle area from 3 points: source, target and next point
ft getTriangleArea(Segment_2& segment, Point& nextPoint){

  std::vector<Point>trianglePoints{segment.source(), segment.target(), nextPoint};

  std::vector<Segment_2> traingleSegments = getConvexHull(trianglePoints);
  
  if (traingleSegments.size() <= 2){
    return 0;
  }
   
  std::vector<Point> triangleArea = getPolyLinePoints(traingleSegments);

  return CGAL::area(triangleArea[0], triangleArea[1], triangleArea[2]);

}


// Pick Random red segment and delete from the vector
Segment_2 chooseVisibleSegment(std::vector<Segment_2> &visibleSegments, Point &nextPoint, ft &area) {

  srand(time(NULL));
  // int randomRedIndex = rand() % visibleSegments.size();

  // area += getTriangleArea(visibleSegments[randomRedIndex], nextPoint);

  // return visibleSegments[randomRedIndex];

  std::vector<std::pair<ft, Segment_2>> areaToSegment;

  for (Segment_2 segment : visibleSegments){

    ft newArea = getTriangleArea(segment, nextPoint);

    areaToSegment.push_back(std::pair<ft, Segment_2>(newArea, segment));
  }

  // ft minArea = std::numeric_limits<ft>::max();
  // Segment_2 chosenSegment;

  // for (std::pair<ft, Segment_2> i : areaToSegment)
  //   if (i.first < minArea) {
  //     minArea = i.first;
  //     chosenSegment = i.second;
  //   }

  // area += minArea;
  // return chosenSegment;

  ft maxArea = -1;
  Segment_2 chosenSegment;

  for (std::pair<ft, Segment_2> i : areaToSegment)
    if (i.first > maxArea){
      maxArea = i.first;
      chosenSegment = i.second;
    }

  area += maxArea;

  return chosenSegment;
}

void forceInsertPoint(std::vector<Segment_2> &polygLine, Point &nextPoint) {

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

  if (!forceInserted) {
    std::cout << "Could not force insert the Point: " << nextPoint << std::endl;
    std::cout << "Algorithm terminated" << std::endl;
  }
}

// Insert the two new segments in the right place in polygon line
void expandPolygonLine(std::vector<Segment_2> &polygLine, Segment_2 &visibleSegment, Point &nextPoint, std::vector<Point> &remainingPoints){

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


//Compares current convex hull with a convex hull that has the new point inserted
// Returns a vector with all the red segments inserted
std::vector<Segment_2> getRedSegments(std::vector<Segment_2> &currConvexHullSegments, std::vector<Segment_2> &nextConvexHullSegments, Point& nextPoint) {

  std::vector<Segment_2> convexRedSegments;

  for (Segment_2 seg : currConvexHullSegments) {

    int segmentCounter = std::count(nextConvexHullSegments.begin(), nextConvexHullSegments.end(), seg);
    
    // If the new convex hull does not have this segment then this segment is RED
    if (!segmentCounter){
      convexRedSegments.push_back(seg);
    }
  }

  // if no new convex hull segments found then red segment is the on that has the point on it
  if (convexRedSegments.size() == 0){
    for (int i = 0; i < nextConvexHullSegments.size(); i++) {
      if (nextConvexHullSegments[i].has_on(nextPoint)) {
        convexRedSegments.push_back(nextConvexHullSegments[i]);
        break;
      }
    }
  }

  return convexRedSegments;
}


std::vector<Segment_2> findVisibleSegments(std::vector<Segment_2> &polygLine, std::vector<Segment_2> &redSegments, Point &nextPoint) {

  std::vector<Segment_2> totalVisibleSegments;

  for (Segment_2 convexSegment : redSegments){

    std::vector<Segment_2> visibleSegments;
    visibleSegments.clear();

    // Convex hull segment == polygon line segment and not collinear check
    if ((std::find(polygLine.begin(), polygLine.end(), convexSegment) != polygLine.end())){ 
      
      if (!CGAL::collinear(nextPoint, convexSegment.point(0), convexSegment.point(1))){
        totalVisibleSegments.push_back(convexSegment);
        continue;
      }
      else {
        continue;
      }
    }

    // Convex hull segment not on polygon line segment
    Segment_2 segment = convexSegment;

    int index = 0;
    
    for (int i = 0; i < polygLine.size(); i++){
      if (polygLine[i].point(0) == segment.point(0)){
        index = i;
        break;
      }
    }

    for (int i = index; i < polygLine.size(); i++){

      if (polygLine[i].point(1) != segment.point(1)){
        visibleSegments.push_back(polygLine[i]);
      }
      else{
        visibleSegments.push_back(polygLine[i]);
        break;
      }
    }

    auto iterator = visibleSegments.begin();

    //If New Point is collinear with a polygon segment then this segment is not visible
    while(iterator != visibleSegments.end()) {
      auto indexing = iterator;
      indexing++;

      if (CGAL::collinear(nextPoint, (*iterator).point(0), (*iterator).point(1))) {
        visibleSegments.erase(iterator);
        indexing--;
      }

      iterator = indexing;
    }

    iterator = visibleSegments.begin();

    while (iterator != visibleSegments.end()) {

      auto indexing = iterator;
      indexing++;
      
      Segment_2 segmentsArray[] = {Segment_2(nextPoint, (*iterator).point(0)),
                                  Segment_2(nextPoint, (*iterator).point(1)),
                                  Segment_2(nextPoint, Point((((*iterator).point(0).x() + (*iterator).point(1).x()) / 2), (((*iterator).point(0).y() + (*iterator).point(1).y()) /2)))
                                  };

      bool intersectionFound = false;
      for (Segment_2 polygSeg : polygLine){
        
        if ((*iterator) == polygSeg) continue;

        for (int k = 0; k < 3; k++) {
          const auto result = intersection(segmentsArray[k], polygSeg);

          if (result) {
            
            if (const Point *p = boost::get<Point>(&*result)) {
              if (*p == (*iterator).point(0) || *p == (*iterator).point(1)) {
                continue;
              } else {
                visibleSegments.erase(iterator);
                indexing--;

                intersectionFound = true;
                break;
              }
            }
            else {
              visibleSegments.erase(iterator);
              indexing--;

              intersectionFound = true;
              break;
            }
          }

        }
        if (intersectionFound)
          break;
      }
      

      iterator = indexing;
    }

    for (Segment_2 visible : visibleSegments)
      totalVisibleSegments.push_back(visible);

  }

  return totalVisibleSegments;
}



ft calcRatio(std::vector<Segment_2> &convexHull, ft &area){

  std::vector<Point> convexHullPoints = getPolyLinePoints(convexHull);

  ft convexHullArea = CGAL::polygon_area_2(convexHullPoints.begin(), convexHullPoints.end(), Convex_hull_traits_2(CGAL::make_property_map(convexHullPoints)));

  return area/convexHullArea;
}




int main() {

  // Read Points from data.instance file
  std::vector<Point> points = readPoints();

  auto start = std::chrono::high_resolution_clock::now();
  // Sort points in a spesific order
  sortPoints(points);

  std::vector<Point> remainingPoints = points;

  // Polygon Line Vector
  std::vector<Segment_2> polygLine;

  // Create A Starting Triangle
  initializeTriangle(polygLine, points, remainingPoints);

  // Get Polygon Line Points To Compute Convex Hull
  std::vector<Point> polygLinePoints = getPolyLinePoints(polygLine);
  
  //Calculate triangle area
  ft area = polygLinePoints.size() == 3 ? CGAL::area(polygLinePoints[0], polygLinePoints[1], polygLinePoints[2]) : 0;

  // Calculate current convex hull
  std::vector<Segment_2> currConvexHullSegments = getConvexHull(polygLinePoints);


  // Loop until there are no remaining points left
  while (remainingPoints.size() > 0) {
   
    // Get next point and delete it from remaining
    Point nextPoint = remainingPoints[0];
    std::cout << "Next Point --> " << nextPoint << std::endl;
    remainingPoints.erase(remainingPoints.begin());

    
    // Push next point
    polygLinePoints.push_back(nextPoint);

    // Calculate convex hull with the new point inside
    std::vector<Segment_2> nextConvexHullSegments = getConvexHull(polygLinePoints);
    
    // To get red segments compare current convex hull and convex hull with the
    // new point if the new novex hull DOES NOT have any segments from the
    // current convex hull then these segments are RED
    std::vector<Segment_2> redSegments = getRedSegments(currConvexHullSegments, nextConvexHullSegments, nextPoint);

    currConvexHullSegments = nextConvexHullSegments;

    while(true){
      
      std::vector<Segment_2>  visibleSegments = findVisibleSegments(polygLine, redSegments, nextPoint);
      
      // if NO visible segments where found then the next point is ON the polygon
      if (visibleSegments.size() == 0) {

        forceInsertPoint(polygLine, nextPoint);
        break;
      }

      Segment_2 visibleSegment = chooseVisibleSegment(visibleSegments, nextPoint, area);
  
      deleteSegment(polygLine, visibleSegment);
      expandPolygonLine(polygLine, visibleSegment, nextPoint, remainingPoints);

      break;

    }

  }

  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::milliseconds polygonizationDuration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  ft ratio = calcRatio(currConvexHullSegments, area);

  Polygon_2 pol_result = Polygon_2();

  for (Segment_2 segment : polygLine) {
    pol_result.push_back(segment.point(0));
  }

  std::cout << "Polygon Is Simple: " << pol_result.is_simple() << std::endl;
  createResultsFile(polygLine, area, polygonizationDuration, ratio);

  return 0;
}
