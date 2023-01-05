#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cstdlib>
#include <map>
#include <iostream>
#include <vector>
#include "../include/preprocessor.hpp"
#include "../include/polygonization.hpp"
#include "../include/convexHull.hpp"
#include "../include/dataio.hpp"
#include "../include/incremental.hpp"
#include "../include/localSearch.hpp"
#include "../include/simulatedAnnealing.hpp"
#include "../include/preprocessor.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;


Preprocessor::Preprocessor(){
}


void Preprocessor::generateSubPoints(){

    // std::vector<std::pair<int, Point>> indexToPoints;

    // int totalPoints = this->percentage * this->points.size();

    // int randomIndex;

    // std::set<int> uniqueIndexes;

    // while(uniqueIndexes.size() < totalPoints){
    //     randomIndex = rand() % this->points.size();
    //     uniqueIndexes.insert(randomIndex);
    // }

    // for (int index : uniqueIndexes){
    //     indexToPoints.push_back(std::pair<int, Point>(index, this->points[index]));
    // }

    // std::sort(indexToPoints.begin(), indexToPoints.end());

    // for (const std::pair<int, Point>& pr : indexToPoints){
    //     this->subPoints.push_back(pr.second);
    // }

}

void Preprocessor::defaultInput(const std::vector<Point>& pointsVec){

    points = pointsVec;

    int pointsSize = points.size();

    // Simulated Annealing Local Step and Global Step
    if (pointsSize <= 1000){
        sizeToLocalL[pointsSize] = 20000;
        sizeToGlobalL[pointsSize] = 3000;
        sizeToSubdivisionL[pointsSize] = 9000;
        sizeToLocalL[10] = 20000;
        sizeToGlobalL[10] = 3000;
        sizeToSubdivisionL[10] = 9000;
        sizeToLocalL[100] = 20000;
        sizeToGlobalL[100] = 3000;
        sizeToSubdivisionL[100] = 9000;
        sizeToLocalL[1000] = 20000;
        sizeToGlobalL[1000] = 3000;
        sizeToSubdivisionL[1000] = 9000;
        sizeToLocalL[1000] = 20000;
        sizeToGlobalL[1000] = 3000;
        sizeToSubdivisionL[1000] = 9000;
    }
    else if (pointsSize <= 10000){
        sizeToLocalL[pointsSize] = 15000;
        sizeToGlobalL[pointsSize] = 1500;
        sizeToSubdivisionL[pointsSize] = 7000;
        sizeToLocalL[10000] = 20000;
        sizeToGlobalL[10000] = 3000;
        sizeToSubdivisionL[10000] = 9000;
    }
    else {
        sizeToLocalL[pointsSize] = 10000;
        sizeToGlobalL[pointsSize] = 1000;
        sizeToSubdivisionL[pointsSize] = 6000;
        sizeToLocalL[100000] = 20000;
        sizeToGlobalL[100000] = 3000;
        sizeToSubdivisionL[100000] = 9000;
    }

}



void Preprocessor::preprocessInput(const std::vector<Point>& pointsVec){
    points =  pointsVec;

    int pointsSize = points.size();
    int starting_L, ending_L;

    Incremental incremental = Incremental(points, 3, "1a");
    incremental.start();

    std::vector<Segment_2> polygon = incremental.getPolygonLine();

    ft prevArea = incremental.getArea();

    if (pointsSize <= 1000){
        starting_L = 10000;
        ending_L = 25000;
    }
    else if (pointsSize <= 10000){
        starting_L = 7000;
        ending_L = 21000;
    }
    else {
        starting_L = 6000;
        ending_L = 21000;
    }

    int optimalLocal_L = 18000;

    for (int L = starting_L; L <= ending_L; L+=100){
      
      simulatedAnnealing simulatedLocal = simulatedAnnealing(points, polygon, incremental.getArea(), incremental.getRatio(), L, 2, 1);
      simulatedLocal.startAnnealing();

      ft optimizedArea = simulatedLocal.getOptimisedArea();

      if (optimizedArea > prevArea){
        prevArea = optimizedArea;
        optimalLocal_L = L;
      }

    }
    std::cout << "OMPTIMAL ==> " << optimalLocal_L << std::endl;
    sizeToLocalL[getSizeBucket(pointsSize)] = (optimalLocal_L + sizeToLocalL[getSizeBucket(pointsSize)]) / 2;
    
    
    if (pointsSize <= 1000){
        starting_L = 8000;
        ending_L = 12000;
    }
    else if (pointsSize <= 10000){
        starting_L = 7000;
        ending_L = 10000;
    }
    else {
        starting_L = 6000;
        ending_L = 8000;
    }

    int optimalsubDiv= starting_L;
    prevArea = 0;

    for (int L = starting_L; L <= ending_L; L+=100){
        simulatedAnnealing simulatedSubdivision = simulatedAnnealing(points, L, 3, 2, "1a", getSimSubDiv_M(pointsSize), 1);
        simulatedSubdivision.startSubdivision();

        ft optimizedArea = simulatedSubdivision.getOptimisedArea();

        if (optimizedArea > prevArea){
            prevArea = optimizedArea;
            optimalLocal_L = L;
        }

    }

    sizeToSubdivisionL[getSizeBucket(pointsSize)] = (optimalsubDiv + sizeToSubdivisionL[getSizeBucket(pointsSize)]) / 2;


    prevArea = incremental.getArea();

    if (pointsSize <= 1000){
        starting_L = 2800;
        ending_L = 3500;
    }
    else if (pointsSize <= 10000){
        starting_L = 1200;
        ending_L = 1800;
    }
    else {
        starting_L = 800;
        ending_L = 1300;
    }

    int optimalGlobal_L = starting_L;
    prevArea = incremental.getArea();

    for (int L = starting_L; L <= ending_L; L+=100){
      
      simulatedAnnealing simulatedGlobal = simulatedAnnealing(points, polygon, incremental.getArea(), incremental.getRatio(), L, 2, 1);
      simulatedGlobal.startAnnealing();

      ft optimizedArea = simulatedGlobal.getOptimisedArea();

      if (optimizedArea > prevArea){
        prevArea = optimizedArea;
        optimalGlobal_L = L;
      }

    }

    sizeToGlobalL[getSizeBucket(pointsSize)] = (optimalGlobal_L + sizeToGlobalL[getSizeBucket(pointsSize)]) / 2;

}


int Preprocessor::getSizeBucket(int size){

    if (sizeToLocalL.find(size) == sizeToLocalL.end()) {
        if (size <= 10){
            size = 10;
        }
        else if (size <= 100){
            size = 100;
        }
        else if (size <= 1000){
            size = 1000;
        }
        else if (size <= 10000){
            size = 10000;
        }
        else if (size <= 100000){
            size = 100000;
        }
    } 
    return size;
}

int Preprocessor::getSimLocal_L(int size){
    return sizeToLocalL[getSizeBucket(size)];
}

int Preprocessor::getSimGlobal_L(int size){
    return sizeToGlobalL[getSizeBucket(size)];
}

int Preprocessor::getSimSubDiv_L(int size){
    return sizeToSubdivisionL[getSizeBucket(size)];
}

int Preprocessor::getSimSubDiv_M(int size){
    // Spatial Subdivision
    if (size <= 10){
        optimal_M = size;
    }
    else{
        optimal_M = size + 0.2 * size;
        if (optimal_M > 100)
            optimal_M = 100;
    }

    return optimal_M;
}