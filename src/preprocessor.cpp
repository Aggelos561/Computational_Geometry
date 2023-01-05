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


std::vector<Point> Preprocessor::generateSubPoints(){

    std::vector<Point> subPoints;

    std::vector<std::pair<int, Point>> indexToPoints;

    int totalPoints = 1000;

    int randomIndex;

    std::set<int> uniqueIndexes;

    while(uniqueIndexes.size() < totalPoints){
        randomIndex = rand() % this->points.size();
        uniqueIndexes.insert(randomIndex);
    }

    for (int index : uniqueIndexes){
        indexToPoints.push_back(std::pair<int, Point>(index, this->points[index]));
    }

    std::sort(indexToPoints.begin(), indexToPoints.end());

    for (const std::pair<int, Point>& pr : indexToPoints){
        subPoints.push_back(pr.second);
    }

    return subPoints;

}

void Preprocessor::defaultInput(){

    // Simulated Annealing Local Step and Global Step
   
    sizeToLocalL[10] = 25000;
    sizeToGlobalL[10] = 3000;
    sizeToSubdivisionL[10] = 10000;
    sizeToLocalL[100] = 21000;
    sizeToGlobalL[100] = 3000;
    sizeToSubdivisionL[100] = 10000;
    sizeToLocalL[1000] = 21000;
    sizeToGlobalL[1000] = 1100;
    sizeToSubdivisionL[1000] = 10000;
    sizeToLocalL[1000] = 21000;
    sizeToGlobalL[1000] = 1100;
    sizeToSubdivisionL[1000] = 10000;


    sizeToLocalL[10000] = 21000;
    sizeToGlobalL[10000] = 1000;
    sizeToSubdivisionL[10000] = 10000;


    sizeToLocalL[100000] = 20000;
    sizeToGlobalL[100000] = 1000;
    sizeToSubdivisionL[100000] = 10000;


}



void Preprocessor::preprocessInput(const std::vector<Point>& pointsVec){
    
    points =  pointsVec;

    int pointsSize = points.size();

    if (pointsSize > 1000)
        points = generateSubPoints();

    std::cout << "Points Size = " << pointsSize << std::endl;
    int starting_L, ending_L;

    Incremental incremental = Incremental(points, 3, "1a");
    incremental.start();

    std::vector<Segment_2> polygon = incremental.getPolygonLine();

    ft prevArea = incremental.getArea();

    if (pointsSize <= 100){
        starting_L = 18000;
        ending_L = 30000;
    } 
    else if (pointsSize <= 1000){
        starting_L = 16000;
        ending_L = 25000;
    }
    else if (pointsSize <= 10000){
        starting_L = 15000;
        ending_L = 20000;
    }
    else {
        starting_L = 15000;
        ending_L = 20000;
    }

    int optimalLocal_L = 18000;

    for (int L = starting_L; L <= ending_L; L+=1000){
      
      simulatedAnnealing simulatedLocal = simulatedAnnealing(points, polygon, incremental.getArea(), incremental.getRatio(), L, 2, 1);
      simulatedLocal.startAnnealing();

      ft optimizedArea = simulatedLocal.getOptimisedArea();

      if (optimizedArea > prevArea){
        prevArea = optimizedArea;
        optimalLocal_L = L;
      }

    }

    
    sizeToLocalL[pointsSize] = (optimalLocal_L + sizeToLocalL[getSizeBucket(pointsSize)]) / 2;    
    sizeToLocalL[getSizeBucket(pointsSize)] = (optimalLocal_L + sizeToLocalL[getSizeBucket(pointsSize)]) / 2;    

    if (pointsSize <= 100){
        starting_L = 10000;
        ending_L = 15000;
    } 
    if (pointsSize <= 1000){
        starting_L = 10000;
        ending_L = 12000;
    }
    else if (pointsSize <= 10000){
        starting_L = 9000;
        ending_L = 12000;
    }
    else {
        starting_L = 10000;
        ending_L = 12000;
    }

    int optimalsubDiv= starting_L;
    prevArea = 0;

    for (int L = starting_L; L <= ending_L; L+=1000){
        simulatedAnnealing simulatedSubdivision = simulatedAnnealing(points, L, 3, 2, "1a", getSimSubDiv_M(pointsSize), 1);
        simulatedSubdivision.startSubdivision();

        ft optimizedArea = simulatedSubdivision.getOptimisedArea();

        if (optimizedArea > prevArea){
            prevArea = optimizedArea;
            optimalLocal_L = L;
        }

    }

    sizeToSubdivisionL[pointsSize] = (optimalsubDiv + sizeToSubdivisionL[getSizeBucket(pointsSize)]) / 2;
    sizeToSubdivisionL[getSizeBucket(pointsSize)] = (optimalsubDiv + sizeToSubdivisionL[getSizeBucket(pointsSize)]) / 2;

    prevArea = incremental.getArea();

    if (pointsSize <= 100){
        starting_L = 4000;
        ending_L = 6000;
    }
    if (pointsSize <= 1000){
        starting_L = 1000;
        ending_L = 1600;
    }
    else if (pointsSize <= 10000){
        starting_L = 800;
        ending_L = 1400;
    }
    else {
        starting_L = 700;
        ending_L = 1300;
    }

    int optimalGlobal_L = starting_L;
    prevArea = incremental.getArea();

    for (int L = starting_L; L <= ending_L; L+=200){
      
      simulatedAnnealing simulatedGlobal = simulatedAnnealing(points, polygon, incremental.getArea(), incremental.getRatio(), L, 2, 1);
      simulatedGlobal.startAnnealing();

      ft optimizedArea = simulatedGlobal.getOptimisedArea();

      if (optimizedArea > prevArea){
        prevArea = optimizedArea;
        optimalGlobal_L = L;
      }

    }

    sizeToGlobalL[pointsSize] = (optimalGlobal_L + sizeToGlobalL[getSizeBucket(pointsSize)]) / 2;
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