#include <vector>
#include <Rcpp.h>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

//typedef boost::geometry::model::d2::point_xy<double> point;
//typedef boost::geometry::model::linestring<point> line;
typedef std::vector<boost::geometry::model::linestring<boost::geometry::model::d2::point_xy<double> > > line_collection;
