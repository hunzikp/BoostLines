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
typedef std::vector<segment> segment_collection;

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
struct dpoint_equal {
  bool operator()(dpoint const &a, dpoint const &b) {
    return a.p.x() == b.p.x() & a.p.y() == b.p.y();
  }
};


// BL MAIN CLASS
class boost_line_collection {
public:

  boost_line_collection(line_collection in_lc) {
    line_collection new_lc = in_lc;
    main_lc = new_lc;
    has_rtree = false;
  }

  line_collection get_line_collection() {
    // Return a copy of the line collection
    line_collection new_lc = main_lc;
    return(main_lc);
  }

  void make_rtree() {
    // Make the rtree (if not already made)

    if (!has_rtree) {

      std::vector<isegment> segment_index_pairs;

      // Segmentize lines
      int segment_counter = 0;
      for(int i = 0; i < main_lc.size(); i++) {
        line this_line = main_lc[i];
        for (int j = 0; j < (this_line.size()-1); j++) {
          segment new_segment;
          new_segment.first = this_line[j];
          new_segment.second = this_line[j+1];

          // Add segment to segment collection
          main_sc.push_back(new_segment);

          // Make segment index pair
          segment_index_pairs.push_back(std::make_pair(new_segment, static_cast<unsigned int>(segment_counter)));

          // Update segment-line-ID map
          segment_map.insert(std::pair<int,int>(segment_counter,i));
          segment_counter++;
        }
      }

      // Create rtree using packing algorithm
      bgi::rtree<isegment, bgi::quadratic<16> > rtree_new(segment_index_pairs);
      main_rtree = rtree_new;
      has_rtree = true;
    }
  }

  segment_collection get_intersecting_segments(line in_line) {
    // Does not return overlapping segments!

    std::vector<isegment> result_n;
    main_rtree.query(bgi::intersects(in_line), std::back_inserter(result_n));
    segment_collection intersecting_segments;
    std::vector<isegment>::iterator itr;
    for (itr = result_n.begin(); itr != result_n.end(); ++itr) {
      isegment value = *itr;
      segment this_segment = main_sc[value.second];
      line this_line;
      bg::convert(this_segment, this_line);
      if (!bg::overlaps(in_line, this_line)) {
        intersecting_segments.push_back(this_segment);
      }
    }

    return intersecting_segments;
  }

  std::vector<int> get_intersecting_line_ids(line in_line) {

    std::vector<isegment> result_n;
    main_rtree.query(bgi::intersects(in_line), std::back_inserter(result_n));
    std::vector<int> line_ids;
    std::vector<isegment>::iterator itr;
    for (itr = result_n.begin(); itr != result_n.end(); ++itr) {
      isegment value = *itr;
      int segment_id = value.second;
      line_ids.push_back(segment_map.at(segment_id) );
    }

    std::sort(line_ids.begin(), line_ids.end());
    std::vector<int>::iterator last = std::unique(line_ids.begin(), line_ids.end());
    line_ids.erase(last, line_ids.end());

    return line_ids;
  }

private:
  line_collection main_lc;

  bool has_rtree;
  bgi::rtree<isegment, bgi::quadratic<16> > main_rtree;
  std::map<int,int> segment_map;  // maps segment IDs (key) onto line IDs (value)
  segment_collection main_sc;
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
bool proper_intersects (const line &l, const line &break_line, const point &startpoint, const point &endpoint) {
  // TRUE if break_line intersects or touches l in between start and end points.

  bool ret = bg::intersects(l, break_line) & !bg::intersects(startpoint, break_line) & !bg::intersects(endpoint, break_line);
  return(ret);

}

line_collection break_line_on_segments (const line &l, const segment_collection &sc) {

  int line_len = l.size();
  int num_segments = sc.size();

  line_collection out_lc;
  line ongoing_line;
  ongoing_line.push_back(l[0]);

  line::const_iterator line_iter;
  line::const_iterator second_to_last_iter = --l.end();
  line::const_iterator third_to_last_iter = --(--l.end());
  for (line_iter = l.begin(); line_iter != second_to_last_iter; ++line_iter) {

    // Make this segment
    line this_segment;

    point segment_start_point = *line_iter;
    line::const_iterator next_iter = boost::next(line_iter, 1);
    point segment_end_point = *next_iter;

    this_segment.push_back(segment_start_point);
    this_segment.push_back(segment_end_point);

    // Get intersection points on segment and their distance from segment start point
    // Also check whether any of the break points lie on segment end point
    bool any_intersects_endpoint = false;
    std::vector<dpoint> points_on_segment;
    for (int k = 0; k < num_segments; k++) {
      segment this_break_segment = sc[k];
      line this_break_line;
      bg::convert(this_break_segment, this_break_line);

      if (proper_intersects(this_segment, this_break_line, segment_start_point, segment_end_point)) {

        std::vector<point> intersection_points;
        bg::intersection(this_segment, this_break_line, intersection_points);
        point this_break_point = intersection_points[0];

        double dist = bg::distance(segment_start_point, this_break_point);
        dpoint this_dpoint(this_break_point, dist);
        points_on_segment.push_back(this_dpoint);

      } else {
        if (!any_intersects_endpoint) {
          any_intersects_endpoint = bg::intersects(segment_end_point, this_break_line);
        }
      }
    }

    // Make new segments in order
    if (points_on_segment.size() > 0) {

      // Sort points on segment
      std::sort(points_on_segment.begin(), points_on_segment.end(), closest());

      // Only keep unique points
      std::vector<dpoint> points_on_segment_copy;
      std::unique_copy(points_on_segment.begin(), points_on_segment.end(),
                       std::back_inserter(points_on_segment_copy), dpoint_equal());
      points_on_segment = points_on_segment_copy;

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
  line_collection out_lc;
  std::vector<int> indices;
  int counter = 1;

  // Make Rtree
  blc.make_rtree();

  // Get all intersecting segments and break
  std::vector<segment_collection > intersecting_segment_vec(in_lc.size(), segment_collection());
  for (int i = 0; i < in_lc.size(); i++) {
    line i_line = in_lc[i];
    segment_collection intersecting_segments = blc.get_intersecting_segments(i_line);
    if (intersecting_segments.size() > 0) {
      line_collection broken_lc = break_line_on_segments(i_line, intersecting_segments);

      for (int j = 0; j < broken_lc.size(); j++) {
        out_lc.push_back(broken_lc[j]);
        indices.push_back(counter);
      }
    } else {
      line copy_line = i_line;
      out_lc.push_back(copy_line);
      indices.push_back(counter);
    }
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
    std::vector<int> j_idx = blc.get_intersecting_line_ids(i_line);
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









