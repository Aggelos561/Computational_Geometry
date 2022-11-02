#include "../include/dataio.hpp"
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <chrono>
#include <iostream>
#include <vector>
#include "../include/polygonization.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point>::type> Convex_hull_traits_2;
typedef CGAL::Epick::FT ft;



Polygonization::Polygonization(std::vector<Point>& points){

	this->points = points;
	this->remainingPoints = points;
	this->duration = (std::chrono::milliseconds) 0;
}


const ft& Polygonization::getArea(){
	return totalArea;
}

const ft& Polygonization::getRatio(){
	return ratio;
}

const std::chrono::milliseconds& Polygonization::getDuration() {
  return duration;
}

const std::vector<Segment_2>& Polygonization::getPolygonLine(){
	return polygLine;
}




void Polygonization::incremental(){

	// Read Points from data.instance file
  auto start = std::chrono::high_resolution_clock::now();
  
	// Sort points in a spesific order
  this->sortPoints(points);

  remainingPoints = this->points;

  // Create A Starting Triangle
  initializeTriangle(polygLine, points, remainingPoints);

  // Get Polygon Line Points To Compute Convex Hull
  polygLinePoints = getPolyLinePoints(polygLine);
  
  //Calculate triangle area
  totalArea = polygLinePoints.size() == 3 ? CGAL::area(polygLinePoints[0], polygLinePoints[1], polygLinePoints[2]) : 0;

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

      Segment_2 visibleSegment = chooseVisibleSegment(visibleSegments, nextPoint, totalArea);
  
      deleteSegment(polygLine, visibleSegment);
      expandPolygonLine(polygLine, visibleSegment, nextPoint);

      break;

    }

  }

  auto stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  ratio = calcRatio(currConvexHullSegments, totalArea);

  Polygon_2 pol_result = Polygon_2();

  for (Segment_2 segment : polygLine) {
    pol_result.push_back(segment.point(0));
  }

  std::cout << "Polygon Is Simple: " << pol_result.is_simple() << std::endl;

}



void Polygonization::convexHullAlg(){

  srand(time(NULL));

 	remainingPoints = points;

  initializeConvexHull(polygLine, points, remainingPoints);

  std::cout << "Size of remainig " << remainingPoints.size() << std::endl;

  for(int i = 0; i < remainingPoints.size(); i++)
    std::cout << "Point " << remainingPoints[i] << std::endl;

 	polygLinePoints = getPolyLinePoints(polygLine);

  std::vector<Segment_2> currConvexHullSegments = getConvexHull(polygLinePoints,remainingPoints);

  initialRun(currConvexHullSegments, remainingPoints, polygLine);
  printf("size before:%d",remainingPoints.size());

  for (int i = 0; i < remainingPoints.size(); i++){
          printf("i=%d remeaing size:%lu\n",i,remainingPoints.size());
          if (forceInsertPoint(polygLine, remainingPoints[i])) {
            printf("seg!%d\n",i);
            std::cout << "Point " << remainingPoints[i] << std::endl;
            remainingPoints.erase(remainingPoints.begin() + i);
            i--;
            printf("insert reamainging,%lu!\n",remainingPoints.size());
          }
  }

  printf("going for second part!\n");
  while(remainingPoints.size() > 0){
    printf("second run!\n");
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

        Point bestPoint = findBestPoint(visPoints, remainingPoints, polygLine, polygLine[i]);
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

}






// Sort Y Ascending
bool Polygonization::sortYAsc(Point &a, Point &b) {
  if (a.y() < b.y())
    return true;
  
  return false;
}

//Sort Y descending
bool Polygonization::sortYDesc(Point &a, Point &b) {
  if (a.y() > b.y())
    return true;

  return false;
}

// Sort points in vector
void Polygonization::sortPoints(std::vector<Point> &points) {
  std::sort(points.begin(), points.end(), this->sortYDesc);
}


// Get all points from a vector of segments
std::vector<Point> Polygonization::getPolyLinePoints(std::vector<Segment_2> &polygLine) {

  std::vector<Point> linePoints;

  for (int i = 0; i < polygLine.size(); i++) {
    linePoints.push_back(polygLine[i].point(0));
  }

  return linePoints;
}

// Find Complex Hull Points and rerurn a vector of segments
std::vector<Segment_2> Polygonization::getConvexHull(std::vector<Point> &polygLinePoints) {

  std::vector<Point> convexHullPoints;

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
void Polygonization::initializeTriangle(std::vector<Segment_2> &polygLine, std::vector<Point> &points, std::vector<Point> &remainingPoints) {

  std::vector<Point> trianglePoints;

  for (int i = 0; i < 3; i++){
    trianglePoints.push_back(points[i]);
  }

  //Trying to create a convex hull triangle that contains 3 points  
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
void Polygonization::deleteSegment(std::vector<Segment_2> &polygLine, Segment_2 &visibleSegment){
  polygLine.erase(std::remove(polygLine.begin(), polygLine.end(), visibleSegment), polygLine.end());
}


// Get triangle area from 3 points: source, target and next point
ft Polygonization::getTriangleArea(Segment_2& segment, Point& nextPoint){

  std::vector<Point>trianglePoints{segment.source(), segment.target(), nextPoint};

  std::vector<Segment_2> traingleSegments = getConvexHull(trianglePoints);
  
  if (traingleSegments.size() <= 2){
    return 0;
  }
   
  std::vector<Point> triangleArea = getPolyLinePoints(traingleSegments);

  return CGAL::area(triangleArea[0], triangleArea[1], triangleArea[2]);

}


// Pick Random red segment and delete from the vector
Segment_2 Polygonization::chooseVisibleSegment(std::vector<Segment_2> &visibleSegments, Point &nextPoint, ft &area) {

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


// if new point already on polygon line the force insert into polygon line
bool Polygonization::forceInsertPoint(std::vector<Segment_2> &polygLine, Point &nextPoint) {

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

  if (forceInserted) {
    return true;
  }
	return false;
}


// Insert the two new segments in the right place in polygon line
void Polygonization::expandPolygonLine(std::vector<Segment_2> &polygLine, Segment_2 &visibleSegment, Point &nextPoint){

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
std::vector<Segment_2> Polygonization::getRedSegments(std::vector<Segment_2> &currConvexHullSegments, std::vector<Segment_2> &nextConvexHullSegments, Point& nextPoint) {

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

// Searching and finding all visible segments based on the new point
std::vector<Segment_2> Polygonization::findVisibleSegments(std::vector<Segment_2> &polygLine, std::vector<Segment_2> &redSegments, Point &nextPoint) {

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


// Calculating area / convex hull area ratio
ft Polygonization::calcRatio(std::vector<Segment_2> &convexHull, ft &area){

  std::vector<Point> convexHullPoints = getPolyLinePoints(convexHull);

  ft convexHullArea = CGAL::polygon_area_2(convexHullPoints.begin(), convexHullPoints.end(), Convex_hull_traits_2(CGAL::make_property_map(convexHullPoints)));

  return area/convexHullArea;
}

std::vector<Segment_2> Polygonization::getConvexHull(std::vector<Point> &polygLinePoints,std::vector<Point> &remainingPoints) {

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



void Polygonization::initializeConvexHull(std::vector<Segment_2> &polygLine, std::vector<Point> &points, std::vector<Point> &remainingPoints) {

  std::vector<Point> convexPoints;

  for (int i = 0; i < points.size(); i++)
    convexPoints.push_back(points[i]);

  polygLine = getConvexHull(convexPoints, remainingPoints);

}




void Polygonization::initialRun(std::vector<Segment_2>& currConvexHullSegments,std::vector<Point>& remainingPoints,std::vector<Segment_2>& polygLine){
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


void Polygonization::findVisiblePoints(std::vector<visPoint>& visPoints, Point& remainingPoint, Segment_2& seg,std::vector<Segment_2>& polygLine){
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



Point Polygonization::findBestPoint(std::vector<visPoint>& visPoints, std::vector<Point>& remainingPoints, std::vector<Segment_2>& polygLine, Segment_2& visSeg){
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
}


void Polygonization::insertBestPoint(std::vector<pair>& bestPoints, std::vector<Point>& remainingPoints,std::vector<Segment_2>& polygLine, int mode){
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
                break;
            }
        }
    }
    else{
        printf("size:%lu\n",remainingPoints.size());
        printf("at least here!\n");
        int rindex = ( rand() % (bestPoints.size() - 1));
        pair bestPair;
        printf("RAND:%lu\n",bestPoints.size());
        bestPair.cor = bestPoints[rindex].cor;
        bestPair.seg = bestPoints[rindex].seg;
        for(int m = 0; m < remainingPoints.size(); m++){
            if(remainingPoints[m] == bestPair.cor){
                remainingPoints.erase(remainingPoints.begin()+m);
                deleteSegment(polygLine, bestPair.seg);
                expandPolygonLine(polygLine, bestPair.seg, bestPair.cor);
                printf("here!\n");
                break;
            }
        }
    }
}
