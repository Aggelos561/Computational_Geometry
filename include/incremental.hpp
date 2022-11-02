#pragma once

#include "polygonization.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Segment_2 Segment_2;
typedef CGAL::Epick::FT ft;

class Incremental: public Polygonization{

	private:
		static bool sortYAsc(Point &, Point &);
		static bool sortYDesc(Point &, Point &);
		void sortPoints(std::vector<Point> &);

		std::vector<Segment_2> getConvexHull(std::vector<Point>&);

		void initializeTriangle(std::vector<Segment_2>&, std::vector<Point>&, std::vector<Point>&);

		ft getTriangleArea(Segment_2&, Point&);

		Segment_2 chooseVisibleSegment(std::vector<Segment_2>&, Point&, ft&);

		std::vector<Segment_2> getRedSegments(std::vector<Segment_2>&, std::vector<Segment_2>&, Point&);

		std::vector<Segment_2> findVisibleSegments(std::vector<Segment_2>&, std::vector<Segment_2>&, Point&);


	public:
		Incremental(std::vector<Point>&);
		void start();
};