#include "../include/dataio.hpp"
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
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
#include "../include/simulatedAnnealing.hpp"
#include "../include/incremental.hpp"
#include "../include/convexHull.hpp"
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





// CLASSIC ANNEALLING LOCAL OR GLOBAL
simulatedAnnealing::simulatedAnnealing(const std::vector<Point>& points, const std::vector<Segment_2>& polygLine, const ft& area, const ft& ratio, int L, int mode, int transitionMode) : Polygonization(points, polygLine, area, ratio){

    this->subDAlgo = 0;
    this->L = L;
    this->optimisedArea = area;
    this->optimisedRatio = ratio;
    this->mode = mode; // 1 = min, 2 = max
    this->transitionMode = transitionMode;
    this->segmentImmunity = false;

}

// ONLY FOR SUBDIVISION - ANNEALING
simulatedAnnealing::simulatedAnnealing(const std::vector<Point>& points, int L, int edgeSelection, int mode, const std::string& initialization, int m,int subDAlgo = 0) : Polygonization(points, edgeSelection){

    this->subDAlgo = subDAlgo;
    this->initialization = initialization;
    this->L = L;
    this->optimisedArea = 0;
    this->optimisedRatio = 0;
    this->mode = mode; // 1 = min, 2 = max
    this->segmentImmunity = false;
    this->m = points.size() < m ? points.size() : m;

}

// ANNEALING WITH MARKED SEGMENTS FOR SUBDIVISION ONLY
simulatedAnnealing::simulatedAnnealing(const std::vector<Point>& points, const std::vector<Segment_2>& polygLine, const ft& area, const ft& ratio, int L, int mode, int transitionMode, const std::vector<Segment_2> untouchable) : Polygonization(points, polygLine, area, ratio){

    this->subDAlgo = 0;
    this->L = L;
    this->optimisedArea = area;
    this->optimisedRatio = ratio;
    this->mode = mode; // 1 = min, 2 = max
    this->transitionMode = transitionMode;
    this->untouchableVector = untouchable;
    this->segmentImmunity = true;

}


// Start Annealing
void simulatedAnnealing::startAnnealing(const std::chrono::_V2::system_clock::time_point startTime, const std::chrono::milliseconds cutOff, const bool measureTime){

  srand(time(NULL));

  Tree tree;
  if (transitionMode == 1)
    KdTreeInit(polygLine, tree);

  double T = 1;

  double R = (double)rand()/(double)RAND_MAX;

  std::vector<Segment_2> convexHullSegments = getConvexHull(getPolyLinePoints(polygLine));
  
  ft convexHullArea = calcArea(convexHullSegments);
  ft energyPrev = energyCalc(convexHullArea, optimisedArea);
  
  int pointsNumber = points.size();

  std::vector<Segment_2> bestPolygon = polygLine;
  ft bestArea = totalArea;

  while(T >= 0){

    if (measureTime && std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime) > cutOff){
      throw cutOffAbort("Cut off time exceeded");
    }

    // transition local or global based on parameter transition mode
    transitionStep transitionStep = transitionMode == 1 ? localTransition(polygLine, tree) : globalTransition(polygLine);

    if (transitionStep.simple){
      ft newEnergy = energyCalc(convexHullArea, optimisedArea + transitionStep.areaDiff);

      if (newEnergy - energyPrev < 0 || metropolis(newEnergy - energyPrev, T) >= R){

        energyPrev = newEnergy;
        polygLine = transitionStep.newPolygLine;
        optimisedArea += transitionStep.areaDiff;

        if ((mode == 1 && optimisedArea < bestArea) || (mode == 2 && optimisedArea > bestArea)){
          bestPolygon = transitionStep.newPolygLine;
          bestArea = optimisedArea;
        }

      }
    }
    else{
      if (pointsNumber <= 20){
        T = T - 1.0/((double)L * 0.7/pointsNumber);
      }

      if (transitionMode == 2 && pointsNumber > 10000){
          T = T - 1.0/((double)L / 2.0);
      }

      continue;
    }

    T = T - 1.0/(double)L;
  }

  polygLine = bestPolygon;
  optimisedArea = bestArea;
  this->optimisedRatio = this->calcRatio(this->getConvexHull(this->points) ,this->optimisedArea);

}




// Initialize and build KD TREE for queries
void simulatedAnnealing::KdTreeInit(const std::vector<Segment_2>& polygLine, Tree& tree){ 

  for(const Point& point : getPolyLinePoints(polygLine)){
    tree.insert(point);
  }

  tree.build();
}



bool simulatedAnnealing::pointsYAscending(const Point& p1, const Point& p2){

  if (p1.y() < p2.y())
    return true;
  
  return false;

}


// Get the rectangle box based on the points
FuzzyBox simulatedAnnealing::getRectangeBox(std::vector<Point>& points){

  std::sort(points.begin(), points.end());

  Point left = points[0], right = points[3];

  std::sort(points.begin(), points.end(), pointsYAscending);

  Point bottom = points[0], top = points[3];

  return FuzzyBox(Point(left.x(),bottom.y()), Point(right.x(),top.y()));

}


// Check if transition is simple based on KD TREE POINTS
bool simulatedAnnealing::validityCheck(const Tree& tree, const std::vector<Segment_2>& polygLine, const Segment_2& nextSegment, const Segment_2& middleSegment, const Segment_2& prevSegment){

  const auto result = intersection(nextSegment, prevSegment);

  if (result)
    return false;


  std::vector<Point> points = { prevSegment.source(), prevSegment.target(), nextSegment.source(), nextSegment.target() };

  FuzzyBox box = getRectangeBox(points);

  std::vector<Point> pointsResult;
  tree.search(std::back_inserter(pointsResult), box);

  std::vector<Point> intersectionPointsVector;
  for (std::vector<Point>::iterator it = pointsResult.begin(); it!= pointsResult.end(); ++it){
    
    if (*it != nextSegment.target() && *it != nextSegment.source() && *it != prevSegment.target() && *it != prevSegment.source() && *it != middleSegment.target() && *it != middleSegment.target()){
      intersectionPointsVector.push_back(*it);
    }
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


// Get area diffrence based on if centroid point of first triangle (CCW segments) is on bounded or not
ft simulatedAnnealing::calcAreaDiff(const std::vector<Segment_2>& polygLine, const Point& PrevSegSource, const Point& PrevSegTarget, const Point& nextSegSource, const Point& nextSegTarget){
  
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


// Enerfy calculation based on mode (minimize maximize)
ft simulatedAnnealing::energyCalc(const ft& convexHullArea, const ft& totalArea){

  int n = polygLine.size();

  //Minimization
  if(mode == 1){
    return n * (totalArea / convexHullArea);
  }
  //Maximization
  else if(mode == 2){
    return n * (1 - totalArea / convexHullArea);
  }

  return 0;
}



// Switch (replace) segments for local transition simulated annealing
void simulatedAnnealing::replace(const Segment_2& prevPolygonPrevSeg, const Segment_2& prevPolygonNextSeg, const Segment_2& middleSegment, const Segment_2& newPolygonNextSeg, const Segment_2& newPolygonPrevSeg, std::vector<Segment_2>& polygLine, int middleSegIndex, int prevPolygonPrevIndex, int prevPolygonNextIndex){

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




transitionStep simulatedAnnealing::localTransition(std::vector<Segment_2>& polygLine, const Tree& tree){

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


// vanilla (k = 1) only for global (only one change)
void simulatedAnnealing::findGlobalChanges(std::vector<Changes>& possibleChanges, std::vector<Point>& points , const Segment_2& e, const std::pair<Point, Point>& pairSequencePoints){

  ft removedArea = calculateDeletedArea(points, e);
  ft addedArea = calculateAddedArea(points, pairSequencePoints);

  Changes change = Changes();

  change.pairPointsSeq = pairSequencePoints;
  change.path = points;
  change.segToRemove = e;
  change.areaDiff = addedArea - removedArea;

  possibleChanges.push_back(change);
  
}


// Find a transition. only for global because at least global has to find a valid k path
transitionStep simulatedAnnealing::globalTransition(std::vector<Segment_2>& polygLine){

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

  } while(!isGlobalValidPath(polygLine, kPoints, e, kPairSequence));

  findGlobalChanges(possibleChanges, kPoints, e, kPairSequence);
  
  return applyGlobalChanges(polygLine, possibleChanges);

}



// Apply the transition found (change found)
transitionStep simulatedAnnealing::applyGlobalChanges(std::vector<Segment_2>& polygLine, std::vector<Changes>& possibleChanges){

  std::vector<Segment_2> changedPolygLine = polygLine;

  transitionStep localTransition;

  if (possibleChanges.size() == 0){
    localTransition.simple = false;
    return localTransition;
  }

  Changes change = possibleChanges[0];

  bool reductionFound = applyKPathRemoval(changedPolygLine, change);

  bool blueRemoved = applyBlueRemoval(changedPolygLine, change);

  if (!isSimple(change, changedPolygLine) || !blueRemoved || reductionFound){
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

// Check if the kpath (k = 1 in global) is valid for transition checking
bool simulatedAnnealing::isGlobalValidPath(const std::vector<Segment_2>& polygLine, const std::vector<Point>& kPoints, const Segment_2& e, const std::pair<Point, Point>& kPairSequence){

  Point source = e.source();
  Point target = e.target();

  if (segmentImmunity){
    for (const Segment_2& segmentUnt : untouchableVector){
      if (looseSegCompare(e, segmentUnt)){
        return false;
      }
    }
  }
  
  for (const Point& point : kPoints){
    if ((point == source) || (point == target))
      return false;
  }

  return true;
}


// Caclulate metropolis
ft simulatedAnnealing::metropolis(const ft& DEnergy, const ft& T){

  return pow(exp(1), (-DEnergy)/T);

}

// Lex order points (if points are not lex ordered in input file). This is for subdivision
bool simulatedAnnealing::lexOrderPoints(const Point& p1, const Point& p2){

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

ft simulatedAnnealing::getOptimisedArea(){
  return this->optimisedArea;
}

ft simulatedAnnealing::getOptimisedRatio(){
  return this->optimisedRatio;
}