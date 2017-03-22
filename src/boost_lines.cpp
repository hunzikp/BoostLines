// [[Rcpp::depends(BH)]]

// [[Rcpp::plugins("cpp11")]]

// MNEMONICS
#include <vector>
#include <Rcpp.h>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/geometry/index/rtree.hpp>
using namespace Rcpp;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


// TYPEDEFS
typedef bg::model::d2::point_xy<double> point;
typedef bg::model::linestring<point> line;
typedef bg::model::segment<point> segment;
typedef std::pair<segment, unsigned int> isegment;
typedef std::vector<line> line_collection;

// STRUCTS
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


// BL MAIN CLASS
class boost_line_collection {
public:

  boost_line_collection(line_collection in_lc) {
    line_collection new_lc = in_lc;
    orig_lc = new_lc;
    has_rtree = false;
  }

  line_collection get_line_collection() {
    line_collection new_lc = orig_lc;
    return(orig_lc);
  }

  void make_rtree() {

    if (!has_rtree) {

      std::vector<isegment> line_index_pairs;

      // Create vector of segment-index pairs
      for(int i = 0; i < orig_lc.size(); i++) {
        line this_line = orig_lc[i];
        for (int j = 0; j < (this_line.size()-1); j++) {
          segment new_segment;
          new_segment.first = this_line[j];
          new_segment.second = this_line[j+1];
          line_index_pairs.push_back(std::make_pair(new_segment, static_cast<unsigned int>(i)));
        }
      }

      // Create rtree using packing algorithm
      bgi::rtree<isegment, bgi::quadratic<16> > rtree_new(line_index_pairs);
      rtree_ = rtree_new;
      has_rtree = true;
    }
  }

  std::vector<int> get_intersecting_ids(line in_line) {

    std::vector<isegment> result_n;
    rtree_.query(bgi::intersects(in_line), std::back_inserter(result_n));
    std::vector<int> indexes;
    std::vector<isegment>::iterator itr;
    for (itr = result_n.begin(); itr != result_n.end(); ++itr) {
      isegment value = *itr;
      indexes.push_back(value.second);
    }

    std::sort(indexes.begin(), indexes.end());
    std::vector<int>::iterator last = std::unique(indexes.begin(), indexes.end());
    indexes.erase(last, indexes.end());

    return indexes;
  }

private:
  line_collection orig_lc;
  bool has_rtree;
  bgi::rtree<isegment, bgi::quadratic<16> > rtree_;
};


// R CONSTRUCTORS
line make_line(NumericMatrix coord) {

  line new_line;
  for (int i=0; i < coord.nrow(); i++) {
    new_line.push_back(bg::make<point>(coord(i,0), coord(i,1)));
  }

  return(new_line);
}

boost_line_collection make_blc(List coord_ls) {
  line_collection new_lc;
  for (int i = 0; i < coord_ls.size(); i++) {
    line new_line = make_line(coord_ls(i));
    new_lc.push_back(new_line);
  }
  boost_line_collection new_blc(new_lc);
  return(new_blc);
}


// R UNPACKERS
List unpack_blc(boost_line_collection blc) {

  line_collection lc = blc.get_line_collection();
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


// NODING FUNCTIONS
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

  line::const_iterator line_iter;
  line::const_iterator second_to_last_iter = --l.end();
  line::const_iterator third_to_last_iter = --(--l.end());
  for (line_iter = l.begin(); line_iter != second_to_last_iter; ++line_iter) {

    line this_segment;

    point segment_start_point = *line_iter;
    line::const_iterator next_iter = boost::next(line_iter, 1);
    point segment_end_point = *next_iter;

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

    if (line_iter == third_to_last_iter) {
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

List node_blc (boost_line_collection blc) {
  // Ignores collinear segments

  line_collection in_lc = blc.get_line_collection();

  // Make Rtree
  blc.make_rtree();

  // Get all intersection points
  std::vector<std::vector<point> > intersection_matrix(in_lc.size(), std::vector<point>());
  for (int i = 0; i < in_lc.size(); i++) {
    line i_line = in_lc[i];
    std::vector<int> j_idx = blc.get_intersecting_ids(i_line);
    if (j_idx.size() > 0) {
      for (int k = 0; k < j_idx.size(); k++) {
        int j = j_idx[k];
        if (j > i) {
          line j_line = in_lc[j];
          if (bg::intersects(i_line, j_line) & !bg::overlaps(i_line, j_line)) {

            std::vector<point> intersection_points;
            bg::intersection(i_line, j_line, intersection_points);

            intersection_matrix[j].insert(intersection_matrix[j].end(), intersection_points.begin(), intersection_points.end());
            intersection_matrix[i].insert(intersection_matrix[i].end(), intersection_points.begin(), intersection_points.end());
          }
        }
      }
    }
  }

  // line_collection::const_iterator i, j;
  // std::vector<std::vector<point> >::iterator mi, mj;
  // for (i = in_lc.begin(), mi = intersection_matrix.begin(); i != in_lc.end(); ++i, ++mi) {
  //   for (j = i, mj = mi; j != in_lc.end(); ++j, ++mj) {
  //
  //     line i_line = *i;
  //     line j_line = *j;
  //
  //     if (bg::intersects(i_line, j_line) & !bg::overlaps(i_line, j_line)) {
  //
  //       std::vector<point> intersection_points;
  //       bg::intersection(i_line, j_line, intersection_points);
  //
  //       (*mj).insert((*mj).end(), intersection_points.begin(), intersection_points.end());
  //       (*mi).insert((*mi).end(), intersection_points.begin(), intersection_points.end());
  //     }
  //   }
  // }

  // Break lines
  line_collection out_lc;
  std::vector<int> indices;

  int counter = 1;
  line_collection::const_iterator lc_it;
  std::vector<std::vector<point> >::iterator ipoint_it;
  for(lc_it = in_lc.begin(), ipoint_it = intersection_matrix.begin(); lc_it != in_lc.end(); ++lc_it, ++ipoint_it) {

    line this_line = *lc_it;
    std::vector<point> intersection_points = *ipoint_it;

    // Break this line on break points
    line_collection broken_lc;
    if (intersection_points.size() > 0) {
      broken_lc = break_line_on_points(this_line, intersection_points);
      for (int j = 0; j < broken_lc.size(); j++) {
        out_lc.push_back(broken_lc[j]);
        indices.push_back(counter);
      }
    } else {
      line copy_line = this_line;
      out_lc.push_back(copy_line);
      indices.push_back(counter);
    }

    // Increment line counter
    counter++;
  }

  IntegerVector indices_nv = wrap(indices);
  boost_line_collection out_blc(out_lc);
  List out_ls = List::create(out_blc, indices_nv);

  return(out_ls);
}



// GEOMETRIC OPERATIONS
LogicalMatrix intersects_blc(boost_line_collection blc) {

  blc.make_rtree();
  line_collection lc = blc.get_line_collection();
  int n_lines = lc.size();
  LogicalMatrix lm(n_lines, n_lines);
  lm.fill(false);
  for (int i = 0; i < n_lines; i++) {
    line i_line = lc[i];
    std::vector<int> j_idx = blc.get_intersecting_ids(i_line);
    if (j_idx.size() > 0) {
      for (int k = 0; k < j_idx.size(); k++) {
        int j = j_idx[k];
        if (j >= i) {
          lm(i,j) = true;
          lm(j,i) = true;
        }
      }
    }
  }
  return(lm);
}

NumericMatrix distance_blc(boost_line_collection blc) {

  line_collection lc = blc.get_line_collection();
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


// RCPP STUFF
RCPP_EXPOSED_CLASS(boost_line_collection)
  RCPP_MODULE(mod) {
    class_<boost_line_collection>( "boost_line_collection" );
    function( "make_blc", &make_blc );
    function( "unpack_blc", &unpack_blc );
    function( "node_blc", &node_blc );
    function( "intersects_blc", &intersects_blc );
    function( "distance_blc", &distance_blc );
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
//
// line_collection break_line_on_line_collection (const line &l, const line_collection &lc) {
//   // Does not break lines with overlapping segments!
//
//   std::vector<point> intersection_points;
//   for (int i = 0; i < lc.size(); i++) {
//     line this_breaker = lc[i];
//     if (bg::intersects(l, this_breaker) & !bg::overlaps(l, this_breaker)) {
//
//       std::vector<point> this_intersection_points;
//       bg::intersection(l, this_breaker, this_intersection_points);
//       intersection_points.insert(intersection_points.end(), this_intersection_points.begin(), this_intersection_points.end());
//     }
//   }
//
//   line_collection broken_lc;
//   if (intersection_points.size() > 0) {
//     broken_lc = break_line_on_points(l, intersection_points);
//   } else {
//     line copy_line = l;
//     broken_lc.push_back(copy_line);
//   }
//
//   return(broken_lc);
// }
//
// List node_lc (line_collection in_lc) {
//   // Ignores collinear segments
//
//   line_collection out_lc;
//   std::vector<int> indices;
//
//   int counter = 1;
//   line_collection::const_iterator lc_it;
//   for(lc_it = in_lc.begin(); lc_it != in_lc.end(); ++lc_it) {
//
//     line this_line = *lc_it;
//     line_collection breaker_lc(in_lc);
//
//     // Get break points for this line
//     std::vector<point> intersection_points;
//     for (int i = 0; i < breaker_lc.size(); i++) {
//       line this_breaker = breaker_lc[i];
//       if (bg::intersects(this_line, this_breaker) & !bg::overlaps(this_line, this_breaker)) {
//
//         std::vector<point> this_intersection_points;
//         bg::intersection(this_line, this_breaker, this_intersection_points);
//         intersection_points.insert(intersection_points.end(), this_intersection_points.begin(), this_intersection_points.end());
//       }
//     }
//
//     // Break this line on break points
//     line_collection broken_lc;
//     if (intersection_points.size() > 0) {
//       broken_lc = break_line_on_points(this_line, intersection_points);
//       for (int j = 0; j < broken_lc.size(); j++) {
//         out_lc.push_back(broken_lc[j]);
//         indices.push_back(counter);
//       }
//     } else {
//       line copy_line = this_line;
//       out_lc.push_back(copy_line);
//       indices.push_back(counter);
//     }
//
//     // Increment line counter
//     counter++;
//   }
//
//   IntegerVector indices_nv = wrap(indices);
//   List out_ls = List::create(out_lc, indices_nv);
//
//   return(out_ls);
// }


