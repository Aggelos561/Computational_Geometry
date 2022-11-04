#pragma once

#include "polygonization.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;


// Class for incremental algorithm
// Inherits from polygonization base class

class Incremental: public Polygonization{

	private:
		std::string initialization;

		static bool sortYAsc(Point &, Point &);
		static bool sortYDesc(Point &, Point &);
		void sortPoints(std::vector<Point> &);
		

		std::vector<Segment_2> getConvexHull(const std::vector<Point>&);

		void initializeTriangle(std::vector<Segment_2>&, const std::vector<Point>&, std::vector<Point>&);

		ft getTriangleArea(const Segment_2&, const Point&);

		Segment_2 chooseVisibleSegment(const std::vector<Segment_2>&, const Point&, ft&);

		std::vector<Segment_2> getRedSegments(const std::vector<Segment_2>&, const std::vector<Segment_2>&, const Point&);

		std::vector<Segment_2> findVisibleSegments(const std::vector<Segment_2>&, const std::vector<Segment_2>&, const Point&);


	public:
		Incremental(const std::vector<Point>&, int, const std::string&);
		void start();
};