#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "../include/preprocessor.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;


Preprocessor::Preprocessor(const std::vector<Point>& points, float percentage){
    this->points = points;
    this->percentage = percentage;
}


void Preprocessor::generateSubPoints(){

    std::vector<std::pair<int, Point>> indexToPoints;

    int totalPoints = this->percentage * this->points.size();

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
        this->subPoints.push_back(pr.second);
    }

}

void Preprocessor::preProcessInput(){
    generateSubPoints();

    // for (const Point& p : this->subPoints){
    //     std::cout << "==> " << p << std::endl;
    // }

    

}

int Preprocessor::getLocal_L(){
    return optimalLocal_L;
}


int Preprocessor::getSimLocal_L(){
    return optimalSimLocal_L;
}


int Preprocessor::getSimGlobal_L(){
    return optimalSimGlobal_L;
}
