#include "../include/dataio.hpp"
#include <CGAL/Bbox_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/enum.h>
#include <CGAL/intersections.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>
#include "../include/polygonization.hpp"

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Splitters.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Kd_tree_rectangle.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point>::type> Convex_hull_traits_2;
typedef CGAL::Epick::FT ft;
typedef CGAL::Epick Cartesian;
typedef CGAL::Search_traits_2<Cartesian> TreeTraits;
typedef CGAL::Kd_tree<TreeTraits> Kd_tree;
typedef Kd_tree::Tree Tree;
typedef CGAL::Fuzzy_iso_box<TreeTraits> FuzzyBox;


// Contructor for base polygonization class
Polygonization::Polygonization(const std::vector<Point>& points, int edgeSelection){
	this->points = points;
	this->remainingPoints = points;
  this->edgeSelection = edgeSelection;
  this->totalArea = 0;
  this->ratio = 0;
}

// Return total area of polygon
const ft& Polygonization::getArea(){
	return totalArea;
}

// Return ratio for results file
const ft& Polygonization::getRatio(){
	return ratio;
}


// Return polygon line vector
const std::vector<Segment_2>& Polygonization::getPolygonLine(){
	return polygLine;
}



// Get all points from a vector of segments
std::vector<Point> Polygonization::getPolyLinePoints(const std::vector<Segment_2> &polygLine) {

  std::vector<Point> linePoints;

  for (int i = 0; i < polygLine.size(); i++) {
    linePoints.push_back(polygLine[i].point(0));
  }

  return linePoints;
}



// Deletes a spesific segment value from a vector of segments
void Polygonization::deleteSegment(std::vector<Segment_2> &polygLine, const Segment_2 &visibleSegment){
  polygLine.erase(std::remove(polygLine.begin(), polygLine.end(), visibleSegment), polygLine.end());
}



// if new point already on polygon line the force insert into polygon line because no visible segment exists
bool Polygonization::forceInsertPoint(std::vector<Segment_2> &polygLine, const Point &nextPoint) {

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

  // Check if point was force inserted
  if (forceInserted) {
    return true;
  }
	return false;
}



// Insert the two new segments in the right place in polygon line for the new point
void Polygonization::expandPolygonLine(std::vector<Segment_2> &polygLine, const Segment_2 &visibleSegment, const Point &nextPoint){

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


ft Polygonization::calcArea(const std::vector<Segment_2>& polygLine) {
  polygLinePoints = this->getPolyLinePoints(polygLine);
  ft totalArea = CGAL::polygon_area_2(polygLinePoints.begin(), this->polygLinePoints.end(), Convex_hull_traits_2(CGAL::make_property_map(polygLinePoints)));

  return totalArea;
}

ft Polygonization::calcPointsArea(std::vector<Point> polygPoints) {

  ft totalArea = abs(CGAL::polygon_area_2(polygPoints.begin(), polygPoints.end(), Convex_hull_traits_2(CGAL::make_property_map(polygPoints))));

  return totalArea;
}


// Calculating area / convex hull area ratio for results
ft Polygonization::calcRatio(const std::vector<Segment_2> &convexHull, const ft &area){

  std::vector<Point> convexHullPoints = getPolyLinePoints(convexHull);

  // Calculating convex hull area and divide with total polygon area
  ft convexHullArea = CGAL::polygon_area_2(convexHullPoints.begin(), convexHullPoints.end(), Convex_hull_traits_2(CGAL::make_property_map(convexHullPoints)));

  return area/convexHullArea;
}


// Find Complex Hull Points and rerurn a vector of segments
std::vector<Segment_2> Polygonization::getConvexHull(const std::vector<Point> &polygLinePoints) {

  std::vector<Point> convexHullPoints;

  CGAL::convex_hull_points_2(polygLinePoints.begin(), polygLinePoints.end(), std::back_inserter(convexHullPoints));

  Polygon_2 convexHullPolygon = Polygon_2();

  for (const Point& p : convexHullPoints) {
    convexHullPolygon.push_back(p);
  }

  std::vector<Segment_2> convexHullSegments;

  // Create A Vector That Stores Convex Hull segments
  for (const Segment_2& seg : convexHullPolygon.edges()) {
    convexHullSegments.push_back(seg);
  }

  return convexHullSegments;
}





std::vector<Point> Polygonization::getPathK(std::vector<Segment_2>& polygLine, int segmentIndex, int k, std::pair<Point, Point>& kPairSequence){
  
  std::vector<Point> k_path;

  int segIndex;
  for (segIndex = segmentIndex; segIndex < polygLine.size(); segIndex++){
    k_path.push_back(polygLine[segIndex].source());

    if (segIndex == (polygLine.size() - 1)){
      segIndex = -1;
    }

    k--;

    if (k == 0)
      break;
  }

  if (segmentIndex - 1 < 0){
    kPairSequence.first = polygLine[polygLine.size() - 1].source();
  }
  else {
    kPairSequence.first = polygLine[segmentIndex - 1].source();
  }
  
  kPairSequence.second = polygLine[segIndex < 0 ? 0 : segIndex].target();
  

 return k_path;

}

bool Polygonization::isValidPath(const std::vector<Segment_2>& polygLine, const std::vector<Point>& kPoints, const Segment_2& e, const std::pair<Point, Point>& kPairSequence){

  Point source = e.source();
  Point target = e.target();

  for (const Point& point : kPoints){
    if ((point == source) || (point == target))
      return false;
  }

  return true;
}


ft Polygonization::calculateDeletedArea(std::vector<Point> points, const Segment_2& e){

  points.push_back(e.target());
  points.push_back(e.source());
 

  return calcPointsArea(points);

}


ft Polygonization::calculateAddedArea(std::vector<Point> points, const std::pair<Point, Point>& pairSequencePoints){

  points.push_back(pairSequencePoints.second);
  points.push_back(pairSequencePoints.first);

  return calcPointsArea(points);

}


void Polygonization::findChanges(std::vector<Changes>& possibleChanges, std::vector<Point>& points , const Segment_2& e, const std::pair<Point, Point>& pairSequencePoints){

  ft removedArea = calculateDeletedArea(points, e);
  
  ft addedArea = calculateAddedArea(points, pairSequencePoints);

  if (addedArea > removedArea){
    Changes change = Changes();

    change.pairPointsSeq = pairSequencePoints;
    change.path = points;
    change.segToRemove = e;
    change.areaDiff = addedArea - removedArea;

    possibleChanges.push_back(change);
  }

}

bool Polygonization::sortAreaChanges(Changes& a, Changes& b){
  
  if (a.areaDiff > b.areaDiff)
    return true;

  return false;
  
}

bool Polygonization::applyBlueRemoval(std::vector<Segment_2>& polygLine, Changes& change){

  bool segmentFound = false;
  for (int i = 0; i < polygLine.size(); i++){

    if (polygLine[i] == change.segToRemove){

      segmentFound = true;
      Segment_2 blueSeg = polygLine[i];

      if (change.path.size() == 1) {

        polygLine.insert(polygLine.begin() + (i + 1) % polygLine.size(), Segment_2(polygLine[i].source(), change.path[0]));
        polygLine.insert(polygLine.begin() + (i + 2) % polygLine.size(), Segment_2(change.path[0], polygLine[i].target()));

      }
      else if (change.path.size() == 2){
        
        for (int j = 0; j < change.path.size(); j++){
          
          if (j == 0){
            polygLine.insert(polygLine.begin() + (i + 1) % polygLine.size(), Segment_2(polygLine[i].source(), change.path[0]));
          }
          else if (j == change.path.size() - 1){
            polygLine.insert(polygLine.begin() + (i + j + 1) % polygLine.size(), Segment_2(change.path[j - 1], change.path[j]));
            polygLine.insert(polygLine.begin() + (i + j + 2) % polygLine.size(), Segment_2(change.path[change.path.size() - 1], polygLine[i].target()));
            
          }
        
        }

      }
      else{

        for (int j = 0; j < change.path.size(); j++){
          
          if (j == 0){
            polygLine.insert(polygLine.begin() + (i + 1) % polygLine.size(), Segment_2(polygLine[i].source(), change.path[0]));
          }
          else if (j == change.path.size() - 1){        
            polygLine.insert(polygLine.begin() + (i + j + 1) % polygLine.size(), Segment_2(change.path[j - 1], change.path[j]));
            polygLine.insert(polygLine.begin() + (i + j + 2) % polygLine.size(), Segment_2(change.path[change.path.size() - 1], polygLine[i].target()));
        
          }
          else{
            polygLine.insert(polygLine.begin() + (i + j + 1) % polygLine.size(), Segment_2(change.path[j - 1], change.path[j]));
          }
        
        }

      }

      polygLine.erase(polygLine.begin() + i);
      break;
    }
  }

  return segmentFound;

}


bool Polygonization::applyKPathRemoval(std::vector<Segment_2>& polygLine, Changes& change){
  
  int startIndex = 0;

  for (int i = 0; i < polygLine.size(); i++){
    if (polygLine[i].target() == change.path[0]){
      startIndex = i;
    }
  }

  int k = change.path.size() + 1;

  while(k != 0){
    
    polygLine.erase(polygLine.begin() + startIndex);
    k--;

    if (startIndex >= polygLine.size() - 1)
      startIndex = 0;
  }
  
  int insertionIndex = 0;
  for (insertionIndex = 0; insertionIndex < polygLine.size(); insertionIndex++){

    if (polygLine[insertionIndex].target() == change.pairPointsSeq.first){

      polygLine.insert(polygLine.begin() + (insertionIndex + 1) % polygLine.size(), Segment_2(change.pairPointsSeq.first, change.pairPointsSeq.second));      
      break;
      
    }
  
  }


  Polygon_2 polygon = Polygon_2();

  for (Segment_2& segment : polygLine){
    polygon.push_back(segment.source());
  }

  for (Point& p : change.path){

    K traits;

    CGAL::Bounded_side result = CGAL::bounded_side_2(polygon.begin(), polygon.end(), p, traits);

    if (result != CGAL::ON_BOUNDED_SIDE){
      return true;
    }

  }

  return false;
}


bool Polygonization::checkPolygonSimplicity(std::vector<Segment_2>& polygLine){

  Polygon_2 polygon = Polygon_2();

  for (Segment_2& segment : polygLine){
    polygon.push_back(segment.source());
  }


  for (int i = 0; i < polygLine.size(); i++){

    if (polygLine[i % polygLine.size()].target() != polygLine[(i+1) % polygLine.size()].source()){
      return false;
    }
  }
  return polygon.is_simple();

}

void Polygonization::applyChanges(std::vector<Segment_2>& polygLine, std::vector<Changes>& possibleChanges){

  std::sort(possibleChanges.begin(), possibleChanges.end(), this->sortAreaChanges);

  std::vector<Segment_2> changedPolygLine = polygLine, prevPolygon = polygLine;

  bool changeApplied = false;

  for (Changes& change : possibleChanges){

    bool reductionFound = applyKPathRemoval(changedPolygLine, change);

    bool blueRemoved = applyBlueRemoval(changedPolygLine, change);

    if (!checkPolygonSimplicity(changedPolygLine) || !blueRemoved || reductionFound){
      changedPolygLine = prevPolygon;
    }
    else{

      changeApplied = true;
      polygLine = changedPolygLine;

      ft newArea = calcArea(polygLine);

      this->areaDiff = abs(this->totalArea + change.areaDiff - this->totalArea);
      this->totalArea += change.areaDiff;

      std::cout << "New Area is " << (long int)newArea << std::endl;
      std::cout << "Blue is " << change.segToRemove << std::endl;
      std::cout << "Red Path is " << change.path[0] << ", " << change.path[1] << std::endl;

      for (Segment_2& segment : changedPolygLine){
        std::cout << "SIMPLE --> " << segment << std::endl;
      }

      std::cout << std::endl;

      break;
      
    }
  }

  if (!changeApplied){
    this->areaDiff = 0;
  }

}


void Polygonization::localSearch(std::vector<Segment_2>& polygLine){

  srand(time(NULL));

  ft Da;
  threshold = 10;
  this->L = 5;

  // polygLine = {Segment_2(Point(40, 20), Point(30, 30)), Segment_2(Point(30, 30), Point(22, 24)), Segment_2(Point(22, 24), Point(19, 22)),
  //               Segment_2(Point(19, 22), Point(15, 26)), Segment_2(Point(15, 26), Point(0, 50)), Segment_2(Point(0, 50), Point(-2, 25)),
  //               Segment_2(Point(-2, 25), Point(16, 4)), Segment_2(Point(16, 4), Point(25, 3)), Segment_2(Point(25, 3), Point(33, 13)),
  //               Segment_2(Point(33, 13), Point(40, 20))};

   for (Segment_2& v : polygLine){
    std::cout << v << std::endl;
   }
  
  std::cout << "Area Before ==> " << (long int)calcArea(polygLine)<< std::endl;
  
  do{

    std::vector<Changes> possibleChanges;
    
    for (Segment_2& e: polygLine){
    
      // Segment_2 e = polygLine[rand() % polygLine.size()];
      int segmentIndex = 0;

      for (Segment_2& v : polygLine){

        for(int k = 1; k <= this->L; k++){

          std::pair<Point, Point> kPairSequence;
          std::vector<Point> kPoints = getPathK(polygLine, segmentIndex, k, kPairSequence);

          if (!isValidPath(polygLine, kPoints, e, kPairSequence))
            break;

          findChanges(possibleChanges, kPoints, e, kPairSequence);

        }
        segmentIndex++;
      }

    }

    applyChanges(polygLine, possibleChanges);

    Da = this->areaDiff;

  } while(Da >= threshold);


  std::cout << std::endl << "Area After ==> " << (long int) calcArea(polygLine) << std::endl;


}



/*---------------------------------------------------------------------------------------------------------------------*/





void Polygonization::KdTreeInit(const std::vector<Segment_2>& polygLine, Tree& tree){ 

  for(const Point& point : getPolyLinePoints(polygLine)){
    tree.insert(point);
  }

  tree.build();
}



bool Polygonization::pointsYAscending(const Point& p1, const Point& p2){

  if (p1.y() < p2.y())
    return true;
  
  return false;

}

FuzzyBox Polygonization::getRectangeBox(std::vector<Point>& points){

  std::sort(points.begin(), points.end());

  Point left = points[0], right = points[3];

  std::sort(points.begin(), points.end(), pointsYAscending);

  Point bottom = points[0], top = points[3];

  return FuzzyBox(Point(left.x(),bottom.y()), Point(right.x(),top.y()));

}


bool Polygonization::validityCheck(const Tree& tree, const std::vector<Segment_2>& polygLine, const Segment_2& nextSegment, const Segment_2& middleSegment, const Segment_2& prevSegment){

  const auto result = intersection(nextSegment, prevSegment);

  if (result)
    return false;

  std::cout << "Middle Segment is " << middleSegment << std::endl;


  std::vector<Point> points = { prevSegment.source(), prevSegment.target(), nextSegment.source(), nextSegment.target() };

  FuzzyBox box = getRectangeBox(points);

  std::vector<Point> pointsResult;
  tree.search(std::back_inserter(pointsResult), box);

  std::vector<Point> intersectionPointsVector;
  for (std::vector<Point>::iterator it = pointsResult.begin(); it!= pointsResult.end(); ++it){
    
    if (*it != nextSegment.target() && *it != nextSegment.source() && *it != prevSegment.target() && *it != prevSegment.source() && *it != middleSegment.target() && *it != middleSegment.target()){
      intersectionPointsVector.push_back(*it);
    }
    std::cout << "==> " << *it << std::endl;
  }
  

  std::vector<Segment_2> intersectionSegments;

  for (int i = 0; i < polygLine.size(); i++){
    for (int j = 0; j < intersectionPointsVector.size(); j++){

      if (polygLine[i].source() == intersectionPointsVector[j] || polygLine[i].target() == intersectionPointsVector[j]){
        if (std::find(intersectionSegments.begin(), intersectionSegments.end(), polygLine[i]) == intersectionSegments.end()){
          intersectionSegments.push_back(polygLine[i]);
        } 
      }
    }

    if (polygLine[i].source() == nextSegment.target() || polygLine[i].target() == prevSegment.source())
      intersectionSegments.push_back(polygLine[i]);
  }


  for (const Segment_2& segment : intersectionSegments){
    
    if (segment.target() != prevSegment.source()){
      
      auto prevInterResult = intersection(segment, prevSegment);
      
      if (prevInterResult){
        return false;
      }

    }
    
    if (segment.source() != nextSegment.target()){
      auto nextInterResult = intersection(segment, nextSegment);

      if  (nextInterResult){
        return false;
      }
    }
  }

  return true;

}


ft Polygonization::calcAreaDiff(const std::vector<Segment_2>& polygLine, const Point& PrevSegSource, const Point& PrevSegTarget, const Point& nextSegSource, const Point& nextSegTarget){
  
  ft addedArea = abs(CGAL::area(PrevSegTarget, nextSegSource, nextSegTarget));
  ft deletedArea = abs(CGAL::area(PrevSegSource, nextSegSource, PrevSegTarget));

  Point triangleCenterP = CGAL::centroid (PrevSegSource, nextSegSource, PrevSegTarget);

  std::vector<Point> points = getPolyLinePoints(polygLine);

  K traits;
  CGAL::Bounded_side result = CGAL::bounded_side_2(points.begin(), points.end(), triangleCenterP, traits);


  if (result == CGAL::ON_BOUNDED_SIDE){
     return addedArea - deletedArea;
  }
  else{
    return deletedArea - addedArea;
  }
  
}


ft Polygonization::energyCalc(const ft& convexHullArea, const ft& totalArea){

  int n = polygLine.size();
  
  // return n * (1 - totalArea / convexHullArea);
  return n * (totalArea / convexHullArea);
  //Minimization
  // if(edgeSelection == 2){
  //   return n * (totalArea / convexHullArea);
  // }
  // //Maximization
  // else if(edgeSelection == 3){
  //   return n * (1 - totalArea / convexHullArea);
  // }

  return 0;
}


void Polygonization::replace(const Segment_2& prevPolygonPrevSeg, const Segment_2& prevPolygonNextSeg, const Segment_2& middleSegment, const Segment_2& newPolygonNextSeg, const Segment_2& newPolygonPrevSeg, std::vector<Segment_2>& polygLine, int middleSegIndex, int prevPolygonPrevIndex, int prevPolygonNextIndex){

  // Insert Segment After Polygon
  polygLine.insert(polygLine.begin() + (prevPolygonNextIndex + 1) % polygLine.size(), newPolygonNextSeg);
  polygLine.erase(polygLine.begin() + prevPolygonNextIndex % polygLine.size());
  
  Segment_2 middleSeg =  polygLine[middleSegIndex];

  Segment_2 newMiddleSeg(middleSeg.target(), middleSeg.source());

  polygLine.insert(polygLine.begin() + middleSegIndex + 1, newMiddleSeg);
  polygLine.erase(polygLine.begin() + middleSegIndex);

  // Insert Segment Before Polygon
  int insert_index = (prevPolygonPrevIndex ) > (polygLine.size() - 1) ? 0 : prevPolygonPrevIndex + 1;
  
  polygLine.insert(polygLine.begin() + insert_index, newPolygonPrevSeg);
  polygLine.erase(polygLine.begin() + prevPolygonPrevIndex);

}

transitionStep Polygonization::localTransition(std::vector<Segment_2>& polygLine, const Tree& tree){

  std::vector<Point> points = getPolyLinePoints(polygLine);

  int q_num = rand() % points.size();

  std::vector<Segment_2> polygLine2 = polygLine;

  for(int i = 0; i < polygLine2.size(); i++){

    if(polygLine2[i].target() == points[q_num]){

      Segment_2 middleSegment = polygLine2[i];
      int prevPolygonPrevIndex = (i-1) < 0 ? polygLine2.size() - 1 : i - 1;

      Segment_2 prevPolygonPrevSeg = polygLine2[prevPolygonPrevIndex];
      Segment_2 prevPolygonNextSeg = polygLine2[(i+1) % polygLine2.size()];

      Segment_2 newPolygonNextSeg(middleSegment.source(), prevPolygonNextSeg.target());
      Segment_2 newPolygonPrevSeg(prevPolygonPrevSeg.source(), middleSegment.target());

      replace(prevPolygonPrevSeg, prevPolygonNextSeg, middleSegment, newPolygonNextSeg, newPolygonPrevSeg, polygLine2, i, prevPolygonPrevIndex, (i+1) % polygLine2.size());
      
      if (validityCheck(tree, polygLine2, newPolygonNextSeg, middleSegment, newPolygonPrevSeg)){
          
          transitionStep local;
          local.newPolygLine = polygLine2;
          local.areaDiff = calcAreaDiff(polygLine, newPolygonPrevSeg.source(), newPolygonPrevSeg.target(), newPolygonNextSeg.source(), newPolygonNextSeg.target());
          std::cout << "Area before is " <<  totalArea << std::endl;
          std::cout << "Area diff is " << local.areaDiff << std::endl;
          std::cout << "Area step calc is " <<  totalArea + local.areaDiff << std::endl;
          std::cout << "Area step real is " <<  calcArea(local.newPolygLine) << std::endl;
          local.simple = true;

          return local;
        
      }
      else{
          transitionStep local;
          local.newPolygLine = polygLine2;
          local.areaDiff = 0;
          local.simple = false;

          return local;
      }
        
    }

  }   


  transitionStep local;
  local.newPolygLine = polygLine2;
  local.areaDiff = 0;
  local.simple = false;

  return local;
}



void Polygonization::findGlobalChanges(std::vector<Changes>& possibleChanges, std::vector<Point>& points , const Segment_2& e, const std::pair<Point, Point>& pairSequencePoints){

  ft removedArea = calculateDeletedArea(points, e);
  ft addedArea = calculateAddedArea(points, pairSequencePoints);

  Changes change = Changes();

  change.pairPointsSeq = pairSequencePoints;
  change.path = points;
  change.segToRemove = e;
  change.areaDiff = addedArea - removedArea;

  possibleChanges.push_back(change);
  
}



transitionStep Polygonization::globalTransition(std::vector<Segment_2>& polygLine){

  int k = 1;
  std::vector<Changes> possibleChanges;
  
  Segment_2 e;
  std::pair<Point, Point> kPairSequence;
  std::vector<Point> kPoints;


  do {
    
    kPoints.clear();

    int s_index = rand() % polygLine.size();

    e =  polygLine[s_index];

    int segmentIndex =  rand() % polygLine.size();

    kPoints = getPathK(polygLine, segmentIndex, k, kPairSequence);
  
  } while(!isValidPath(polygLine, kPoints, e, kPairSequence));

  findGlobalChanges(possibleChanges, kPoints, e, kPairSequence);
  
  return applyGlobalChanges(polygLine, possibleChanges);

}




transitionStep Polygonization::applyGlobalChanges(std::vector<Segment_2>& polygLine, std::vector<Changes>& possibleChanges){

  std::vector<Segment_2> changedPolygLine = polygLine;

  transitionStep localTransition;

  if (possibleChanges.size() == 0){
    localTransition.simple = false;
    return localTransition;
  }

  Changes change = possibleChanges[0];

  bool reductionFound = applyKPathRemoval(changedPolygLine, change);

  bool blueRemoved = applyBlueRemoval(changedPolygLine, change);

  if (!checkPolygonSimplicity(changedPolygLine) || !blueRemoved || reductionFound){
    localTransition.simple = false;
    return localTransition;
  }
  else{

    localTransition.areaDiff = change.areaDiff;
    localTransition.newPolygLine = changedPolygLine;
    localTransition.simple = true;

    return localTransition;
    
  }


}



ft Polygonization::metropolis(const ft& DEnergy, const ft& T){

  return pow(exp(1), (-DEnergy)/T);

}

bool Polygonization::lexOrderPoints(const Point& p1, const Point& p2){

  if (p2.x() > p1.x())
    return true;

  else if (p1.x() > p2.x())
    return false;
  
  else{

    if (p2.y() > p1.y())
      return true;
    
    else
     return false;
  }

}


void Polygonization::spatialSubdivision(std::vector<Point>& points, int edge_selection, const std::string& initialization){

  std::sort(points.begin(), points.end(), lexOrderPoints);

  for (const Point& p : points){
    std::cout << p << std::endl;
  }

  int m = 10;
  int n = points.size();

  int k = std::ceil((double)(n - 1)/(double)(m - 1));

  std::vector<std::vector<Point>> subPolPoints;

  std::cout << "K = " << k << std::endl;

  int threshold = 5;
  int pIndex = 0;

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


  // std::cout << "Sub Teams Size = " << subPolPoints.size() << std::endl;

  int lastPointsSetIndex = subPolPoints.size() - 1;
  std::cout << "Last vector points size is " << subPolPoints[lastPointsSetIndex].size() << std::endl;
  if (subPolPoints[lastPointsSetIndex].size() < threshold){
    for (int i = 1; i < subPolPoints[lastPointsSetIndex].size(); i++){
      subPolPoints[lastPointsSetIndex - 1].push_back(subPolPoints[lastPointsSetIndex][i]);
    }
    subPolPoints.erase(subPolPoints.begin() + lastPointsSetIndex);
    k--;
  }

  // Prints points subsets
  // for (int i = 0; i < subPolPoints.size(); i++){
  //   for (int j = 0; j < subPolPoints[i].size(); j++){
  //     std::cout << "p --> " << subPolPoints[i][j] << std::endl;
  //   }
  //   std::cout << std::endl;
  // }

}

void Polygonization::simulatedAnnealing(std::vector<Segment_2>& polygLine){

  // Print initial polygon line
  for (const Segment_2& segment : polygLine){
    std::cout << segment << std::endl;
  }

  std::cout << std::endl << std::endl;

  long int areaBefore = calcArea(polygLine);
  srand(time(NULL));

  Tree tree;
  KdTreeInit(polygLine, tree);

  double T = 1;

  double R = (double)rand()/(double)RAND_MAX;

  std::cout << "R = " << R << std::endl;

  this->L = 1000;

  totalArea = calcArea(polygLine); 
  std::vector<Segment_2> convexHullSegments = getConvexHull(getPolyLinePoints(polygLine));
  ft convexHullArea = calcArea(convexHullSegments);

  ft energyPrev = energyCalc(convexHullArea, totalArea);

  while(T >= 0){

    //Local transition 
    transitionStep transitionStep = localTransition(polygLine, tree);

    // transitionStep transitionStep = globalTransition(polygLine);

    if (transitionStep.simple){
      ft newEnergy = energyCalc(convexHullArea, totalArea + transitionStep.areaDiff);
      
      if (newEnergy - energyPrev < 0 || metropolis(newEnergy - energyPrev, T) >= R){
        
        std::cout << "Found Better Polygon!" << std::endl;

        for (const Segment_2& segment : transitionStep.newPolygLine){
          std::cout << segment << std::endl;
       }
        std::cout << std::endl;

        energyPrev = newEnergy;
        polygLine = transitionStep.newPolygLine;
        totalArea += transitionStep.areaDiff;
      }
    }
    else{
      continue;
    }

    
    T = T - 1.0/(double)L;
  }

  Polygon_2 polygon;

  for (const Segment_2& segment : polygLine){
    polygon.push_back(segment.source());
  }

  std::cout << "Final Check Simplicity: " << polygon.is_simple() << std::endl;
  long int areaAfter = calcArea(polygLine);  
  std::cout << "Area Before ==> " << areaBefore << std::endl;
  std::cout << " Real Area After ==> " << areaAfter << std::endl;
  std::cout << " Calculated Area After ==> " << (long int)totalArea << std::endl;
}


// Obtain greedy solution S as initial state; compute its "energy" E.
// 2: T ← 1
// 3: while T ≥ 0 do
// Perform transition step: local xor global according to parameter.
// 5:
// Check validity. If not valid goto 4.
// 6:
// If ∆E < 0 or Metropolis criterion holds, apply the transition.
// 7:
// T =T−1/L
// 8: end while