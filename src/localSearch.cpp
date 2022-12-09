#include "../include/dataio.hpp"
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/enum.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include "../include/polygonization.hpp"
#include "../include/localSearch.hpp"


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



localSearch::localSearch(const std::vector<Point>& points, const std::vector<Segment_2>& polygLine, const ft& area, const ft& ratio, int L, double threshold, int mode) : Polygonization(points, polygLine, area, ratio) {
    this->optimisedArea = area;
    this->optimisedRatio = ratio;
    this->L = L;
    this->threshold = threshold;
    this->mode = mode;//Minimazation or Maximazation
}




void localSearch::start(){

  srand(time(NULL));

  ft Da;

  do{

    std::vector<Changes> possibleChanges;
    
    for (Segment_2& e: polygLine){
    
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

    Da = abs(this->areaDiff);


  } while(Da >= threshold);

  this->optimisedRatio = this->calcRatio(this->getConvexHull(this->points),this->optimisedArea);
}



void localSearch::applyChanges(std::vector<Segment_2>& polygLine, std::vector<Changes>& possibleChanges){

  std::sort(possibleChanges.begin(), possibleChanges.end(), this->sortAreaChanges);

  std::vector<Segment_2> changedPolygLine = polygLine, prevPolygon = polygLine;

  bool changeApplied = false;

  for (Changes& change : possibleChanges){

    bool reductionFound = applyKPathRemoval(changedPolygLine, change);

    bool blueRemoved = applyBlueRemoval(changedPolygLine, change);

    if (!isSimple(change, changedPolygLine) || !blueRemoved || reductionFound){
      changedPolygLine = prevPolygon;
    }
    else{

      changeApplied = true;
      polygLine = changedPolygLine;

      this->areaDiff = change.areaDiff;
      this->optimisedArea += this->areaDiff;

      break;
      
    }
  }

  if (!changeApplied){
    this->areaDiff = 0;
  }

}




void localSearch::findChanges(std::vector<Changes>& possibleChanges, std::vector<Point>& points , const Segment_2& e, const std::pair<Point, Point>& pairSequencePoints){

  ft removedArea = calculateDeletedArea(points, e);
  
  ft addedArea = calculateAddedArea(points, pairSequencePoints);

  bool optimisationCond = mode == 1 ? removedArea > addedArea : addedArea > removedArea;

  if (optimisationCond){
    Changes change = Changes();

    change.pairPointsSeq = pairSequencePoints;
    change.path = points;
    change.segToRemove = e;
    change.areaDiff = addedArea - removedArea;

    possibleChanges.push_back(change);
  }

}



bool localSearch::sortAreaChanges(Changes& a, Changes& b){
  
  if (a.areaDiff > b.areaDiff)
    return true;

  return false;
  
}

ft localSearch::getOptimisedArea(){
  return this->optimisedArea;
}

ft localSearch::getOptimisedRatio(){
  return this->optimisedRatio;
}