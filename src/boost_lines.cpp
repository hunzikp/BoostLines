// [[Rcpp::depends(BH)]]


// MNEMONICS
#include <vector>
#include <Rcpp.h>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/assign/list_of.hpp>
using namespace Rcpp;
namespace bg = boost::geometry;


// TYPEDEFS & CLASSES
typedef bg::model::d2::point_xy<double> point;
typedef bg::model::linestring<point> line;
class line_collection : public std::vector<line>{
};
struct dpoint {

  dpoint(point inp, double indist) {
    p = inp;
    dist = indist;
  }

  point p;
  double dist;
};
struct closest {
  bool operator()(dpoint const &a, dpoint const &b) {
    return a.dist < b.dist;
  }
};



// R CONSTRUCTORS
line make_line(NumericMatrix coord) {

  line new_line;
  for (int i=0; i < coord.nrow(); i++) {
    new_line.push_back(bg::make<point>(coord(i,0), coord(i,1)));
  }

  return(new_line);
}

line_collection make_lc(List coord_ls) {

  line_collection lc;
  for (int i = 0; i < coord_ls.size(); i++) {
    line new_line = make_line(coord_ls(i));
    lc.push_back(new_line);
  }
  return(lc);
}



// NODING
bool point_within (const point &p, const line & l, const point &startpoint, const point &endpoint) {

  bool ret = bg::intersects(p, l) & !bg::intersects(startpoint, p) & !bg::intersects(endpoint,p);
  return(ret);

}

line_collection break_line_on_points (const line &l, const std::vector<point> &pvec) {

  std::vector<point> pvec_unique;
  std::unique_copy(pvec.begin(), pvec.end(), std::back_inserter(pvec_unique), boost::geometry::equal_to<point>());

  int line_len = l.size();
  int num_points = pvec_unique.size();

  line_collection out_lc;
  line ongoing_line;
  ongoing_line.push_back(l[0]);

  for (int i = 0; i < line_len - 1; i++) {
    line this_segment;
    point segment_start_point = l[i];
    point segment_end_point = l[i+1];
    this_segment.push_back(segment_start_point);
    this_segment.push_back(segment_end_point);

    // Get points on segment and their distance from this point
    // Also check whether any of the break points lie on segment end point
    bool any_intersects_endpoint = false;
    std::vector<dpoint> points_on_segment;
    for (int k = 0; k < num_points; k++) {
      point this_break_point = pvec_unique[k];

      if (point_within(this_break_point, this_segment, segment_start_point, segment_end_point)) {
        double dist = bg::distance(segment_start_point, this_break_point);
        dpoint this_dpoint(this_break_point, dist);
        points_on_segment.push_back(this_dpoint);
      } else {
        if (!any_intersects_endpoint) {
          any_intersects_endpoint = bg::intersects(segment_end_point, this_break_point);
        }
      }
    }

    // Make new segments in order
    if (points_on_segment.size() > 0) {

      // Sort points on segment
      std::sort(points_on_segment.begin(), points_on_segment.end(), closest());

      for (int k = 0; k < points_on_segment.size(); k++) {
        // Checkout current line
        point this_break_point = points_on_segment[k].p;
        ongoing_line.push_back(this_break_point);
        line finished_line = ongoing_line;
        out_lc.push_back(finished_line);
        // Start new line
        ongoing_line.clear();
        ongoing_line.push_back(this_break_point);
      }
    }

    // Add segment end point to ongoing line
    ongoing_line.push_back(segment_end_point);

    if (i == line_len - 2) {
      // If this is the last segment: checkout current line
      line finished_line = ongoing_line;
      out_lc.push_back(finished_line);
    } else {
      // If a break point intersects with segment end point:
      // Checkount ongoing line and start new one
      if (any_intersects_endpoint) {
        line finished_line = ongoing_line;
        out_lc.push_back(finished_line);
        ongoing_line.clear();
        ongoing_line.push_back(segment_end_point);
      }
    }
  }

  return(out_lc);
}

// line_collection break_line_on_line (const line &l1, const line &l2) {
//   // Does not break lines with overlapping segments!
//
//   line_collection broken_lc;
//   if (bg::intersects(l1, l2) & !bg::overlaps(l1, l2)) {
//
//     std::vector<point> intersection_points;
//     bg::intersection(l1, l2, intersection_points);
//     broken_lc = break_line_on_points(l1, intersection_points);
//
//   } else {
//
//     line copy_line = l1;
//     broken_lc.push_back(copy_line);
//
//   }
//
//   return(broken_lc);
// }
//
//
// line_collection break_line_on_line_collection (const line &l, const line_collection &breaker_lc) {
//   // Ignores collinear segments
//
//   line_collection input_lc;
//   line_collection output_lc;
//   input_lc.push_back(l);
//
//   for (int i = 0; i < breaker_lc.size(); i++) {
//
//     line this_breaker = breaker_lc[i];
//
//     for (int j = 0; j < input_lc.size(); j++) {
//
//       line this_line = input_lc[j];
//
//       if (bg::intersects(this_line, this_breaker)) {
//
//         line_collection broken_blc = break_line_on_line(this_line, this_breaker);
//         for (int k = 0; k < broken_blc.size(); k++) {
//           output_lc.push_back(broken_blc[k]);
//         }
//
//       } else {
//         output_lc.push_back(this_line);
//       }
//     }
//
//     input_lc = output_lc;
//     output_lc.clear();
//   }
//
//   return(input_lc);
// }

line_collection break_line_on_line_collection (const line &l, const line_collection &lc) {
  // Does not break lines with overlapping segments!

  std::vector<point> intersection_points;
  for (int i = 0; i < lc.size(); i++) {
    line this_breaker = lc[i];
    if (bg::intersects(l, this_breaker) & !bg::overlaps(l, this_breaker)) {

      std::vector<point> this_intersection_points;
      bg::intersection(l, this_breaker, this_intersection_points);
      intersection_points.insert(intersection_points.end(), this_intersection_points.begin(), this_intersection_points.end());
    }
  }

  line_collection broken_lc;
  if (intersection_points.size() > 0) {
    broken_lc = break_line_on_points(l, intersection_points);
  } else {
    line copy_line = l;
    broken_lc.push_back(copy_line);
  }

  return(broken_lc);
}

List node_lc (const line_collection &in_lc) {
  // Ignores collinear segments

  int n_lines = in_lc.size();
  line_collection out_lc;
  std::vector<int> indices;

  for (int i = 0; i < n_lines; i++) {

    line this_line = in_lc[i];
    line_collection breaker_lc(in_lc);
    breaker_lc.erase(breaker_lc.begin() + i);

    line_collection broken_lc = break_line_on_line_collection(this_line, breaker_lc);

    for (int j = 0; j < broken_lc.size(); j++) {
      out_lc.push_back(broken_lc[j]);
      indices.push_back(i+1); // R indexing
    }
  }

  IntegerVector indices_nv = wrap(indices);
  List out_ls = List::create(out_lc, indices_nv);

  return(out_ls);
}


// OTHER GEOMETRIC OPERATIONS
LogicalMatrix intersects_lc(const line_collection &lc) {

  int n_lines = lc.size();
  LogicalMatrix lm(n_lines, n_lines);
  lm.fill(false);
  for (int j = 0; j < n_lines; j++) {
    for (int i = j; i < n_lines; i++) {
      bool cr;
      if (i != j) {
        line l_i = lc[i];
        line l_j = lc[j];
        cr = bg::intersects(l_i, l_j);
      } else {
        cr = true;
      }
      lm(i,j) = cr;
      lm(j,i) = cr;
    }
  }

  return(lm);
}

NumericMatrix distance_lc(const line_collection &lc) {

  int n_lines = lc.size();
  NumericMatrix dm(n_lines, n_lines);
  dm.fill(0);
  for (int j = 0; j < n_lines; j++) {
    for (int i = j; i < n_lines; i++) {
      if (i != j) {
        line l_i = lc[i];
        line l_j = lc[j];
        double d = bg::distance(l_i, l_j);
        dm(i,j) = d;
        dm(j,i) = d;
      }
    }
  }

  return(dm);
}



// R UNPACKERS
List unpack_lc(const line_collection &lc) {

  int n_lines = lc.size();
  List coord_ls(n_lines);

  for (int i = 0; i < n_lines; i++) {
    line this_line = lc[i];
    int len = this_line.size();
    NumericMatrix this_coord(len, 2);
    for (int j = 0; j < len; j++) {
      this_coord(j,0) = this_line[j].x();
      this_coord(j,1) = this_line[j].y();
    }
    coord_ls(i) = this_coord;
  }

  return coord_ls;
}


// RCPP STUFF
RCPP_EXPOSED_CLASS(line_collection)
RCPP_MODULE(mod) {
  class_<line_collection>( "line_collection" );
  function( "make_lc", &make_lc );
  function( "unpack_lc", &unpack_lc );
  function( "node_lc", &node_lc );
  function( "intersects_lc", &intersects_lc );
  function( "distance_lc", &distance_lc );
}
