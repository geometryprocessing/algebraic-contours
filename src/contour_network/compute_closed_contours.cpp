// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "compute_closed_contours.h"

// Compute the squared distance between two contours
double
contour_segment_squared_distance(
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  int contour_segment_1,
  int contour_segment_2)
{
  // Get distance between the contour start and end points
  SpatialVector first_segment_end_point =
    contour_segments[contour_segment_1].end_point();
  SpatialVector second_segment_start_point =
    contour_segments[contour_segment_2].start_point();
  SpatialVector displacement =
    second_segment_start_point - first_segment_end_point;
  double squared_distance = displacement.dot(displacement);

  return squared_distance;
}

// Return true iff contour segment 2 follows contour segment 1
bool
is_adjacent_contour_segment(
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  int contour_segment_1,
  int contour_segment_2)
{
  // Get distance between the contour start and end points
  double squared_distance = contour_segment_squared_distance(
    contour_segments, contour_segment_1, contour_segment_2);

  return float_equal_zero(squared_distance);
}

// Return true iff the two points overlap
bool
are_overlapping_points(const SpatialVector& point_1,
                       const SpatialVector& point_2)
{
  SpatialVector displacement = point_1 - point_2;
  double squared_distance = displacement.dot(displacement);
  return float_equal_zero(squared_distance);
}

// Return true iff the contour chain is closed
bool
is_closed_contour(const std::vector<int>& contour,
                  const std::vector<SpatialVector>& contour_start_points,
                  const std::vector<SpatialVector>& contour_end_points)
{
  // Treat empty contour as closed
  if (contour.empty())
    return true;

  // Check if the start point of the first segment overlaps the end point of the
  // last segment
  return are_overlapping_points(contour_end_points[contour.back()],
                                contour_start_points[contour.front()]);
}

// Return true iff contours defines valid contour chains
bool
is_valid_contours(const std::vector<std::vector<int>>& contours,
                  const std::vector<SpatialVector>& contour_start_points,
                  const std::vector<SpatialVector>& contour_end_points)
{
  // Check contour segments are contiguous
  for (size_t i = 0; i < contours.size(); ++i) {
    for (size_t j = 1; j < contours[i].size(); ++j) {
      if (!are_overlapping_points(contour_end_points[contours[i][j - 1]],
                                  contour_start_points[contours[i][j]])) {
        spdlog::error(
          "Segment {} in contour {} not adjacent to segment {}", j - 1, i, j);
        return false;
      }
    }
  }

  return true;
}

// Get unused contour segment
int
get_starting_contour_segment(const std::vector<bool>& used_segments)
{
  // Find free segment
  for (size_t i = 0; i < used_segments.size(); ++i) {
    if (!used_segments[i])
      return i;
  }

  // Return -1 if no free segment
  return -1;
}

// Add next contour segment if possible and return false if not
bool
add_next_contour_segment(
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<int>& current_contour,
  std::vector<bool>& used_segments,
  const std::vector<SpatialVector>& contour_start_points,
  const std::vector<SpatialVector>& contour_end_points)
{
  for (size_t i = 0; i < contour_segments.size(); ++i) {
    // Skip used segments
    if (used_segments[i])
      continue;

    // Determine if the candidate contour segment is adjacent to the growing
    // chain
    if (are_overlapping_points(contour_end_points[current_contour.back()],
                               contour_start_points[i])) {
      current_contour.push_back(i);
      used_segments[i] = true;
      return true;
    }
  }

  // Return false if no segment found
  return false;
}

// Add previous contour segment to reverse list if possible and return false if
// not
bool
add_next_reverse_contour_segment(
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<int>& current_contour_reverse,
  std::vector<bool>& used_segments,
  const std::vector<SpatialVector>& contour_start_points,
  const std::vector<SpatialVector>& contour_end_points)
{
  for (size_t i = 0; i < contour_segments.size(); ++i) {
    // Skip used segments
    if (used_segments[i])
      continue;

    // Determine if the candidate contour segment is adjacent to the growing
    // chain
    if (are_overlapping_points(
          contour_end_points[i],
          contour_start_points[current_contour_reverse.back()])) {
      current_contour_reverse.push_back(i);
      used_segments[i] = true;
      return true;
    }
  }

  // Return false if no segment found
  return false;
}

// Combine a forward and reverse contour into one contour
void
combine_forward_and_reverse_contour(const std::vector<int>& contour,
                                    const std::vector<int>& contour_reverse,
                                    std::vector<int>& contour_full)
{
  assert(contour.size() > 0);
  assert(contour_reverse.size() > 0);
  assert(contour[0] == contour_reverse[0]);

  // Reserve space for the full contour
  size_t num_segments = contour.size() + contour_reverse.size() - 1;
  contour_full.reserve(num_segments);

  // Add reverse contour segments in reverse (excluding the last one as it is
  // redundant)
  for (size_t i = contour_reverse.size() - 1; i > 0; --i) {
    contour_full.push_back(contour_reverse[i]);
  }

  // Add contour segments
  for (size_t i = 0; i < contour.size(); ++i) {
    contour_full.push_back(contour[i]);
  }
}

// Add contour to the list of all contours and assign it a new label
void
add_contour(std::vector<int>& contour,
            std::vector<std::vector<int>>& contours,
            std::vector<int>& contour_labels)
{
  contours.push_back(contour);
  for (size_t i = 0; i < contour.size(); ++i) {
    contour_labels[contour[i]] = contours.size();
  }
}

void
compute_closed_contours(
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<std::vector<int>>& contours,
  std::vector<int>& contour_labels)
{
  contours.clear();
  contour_labels.clear();
  std::vector<bool> used_segments(contour_segments.size(), false);
  contour_labels.resize(contour_segments.size());

  std::vector<SpatialVector> contour_start_points(contour_segments.size());
  std::vector<SpatialVector> contour_end_points(contour_segments.size());

  for (size_t i = 0; i < contour_segments.size(); i++) {
    contour_start_points[i] = contour_segments[i].start_point();
    contour_end_points[i] = contour_segments[i].end_point();
  }

  while (true) {
    // Get next starting contour segment to process or return if none left
    int starting_segment_index = get_starting_contour_segment(used_segments);
    if (starting_segment_index == -1) {
      assert(
        is_valid_contours(contours, contour_start_points, contour_end_points));
      return;
    }

    // Initialize the next contour chain
    std::vector<int> current_contour = { starting_segment_index };
    used_segments[starting_segment_index] = true;

    // Traverse forward until the contour is closed or no new contour is found
    bool closed_contour = true; // Assume closed until proven otherwise
    while (!is_closed_contour(
      current_contour, contour_start_points, contour_end_points)) {
      bool adjacent_segment_found =
        add_next_contour_segment(contour_segments,
                                 current_contour,
                                 used_segments,
                                 contour_start_points,
                                 contour_end_points);

      // If no segment found, the contour is open
      if (!adjacent_segment_found) {
        closed_contour = false;
        break;
      }
    }

    if (closed_contour) {
      spdlog::debug("Closed contour of size {} found", current_contour.size());
      add_contour(current_contour, contours, contour_labels);
    }
    // Traverse an open contour in reverse to get the full chain
    else {
      std::vector<int> current_contour_reverse = { starting_segment_index };
      std::vector<int> current_contour_full(0);
      while (true) {
        bool adjacent_segment_found =
          add_next_reverse_contour_segment(contour_segments,
                                           current_contour_reverse,
                                           used_segments,
                                           contour_start_points,
                                           contour_end_points);

        // If no further segment found, combine the two chains and add it to the
        // list
        if (!adjacent_segment_found) {
          combine_forward_and_reverse_contour(
            current_contour, current_contour_reverse, current_contour_full);
          spdlog::debug("Open contour of size {} found",
                        current_contour_full.size());
          add_contour(current_contour_full, contours, contour_labels);
          break;
        }
      }
    }

    assert(
      is_valid_contours(contours, contour_start_points, contour_end_points));
  }
}

// Compute the distances between contour endpoints
void
compute_contour_closure_distances(
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<std::vector<int>>& contours,
  std::vector<double>& contour_closure_distances)
{
  size_t num_contour_segments = contour_segments.size();
  size_t num_contours = contours.size();
  contour_closure_distances.clear();
  contour_closure_distances.reserve(num_contour_segments);
  spdlog::info("Computing distances for {} contours", num_contours);

  // Iterate over contours and compute distances between the adjacent segments
  for (size_t contour_index = 0; contour_index < num_contours;
       ++contour_index) {
    for (size_t i = 1; i < contours[contour_index].size(); ++i) {
      size_t contour_segment_1 = contours[contour_index][i - 1];
      size_t contour_segment_2 = contours[contour_index][i];
      double squared_distance = contour_segment_squared_distance(
        contour_segments, contour_segment_1, contour_segment_2);
      contour_closure_distances.push_back(std::sqrt(squared_distance));
    }
  }

  spdlog::info("{} contour closure distances found",
               contour_closure_distances.size());
}

void
write_contour_closure_distances(
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<std::vector<int>>& contours,
  const std::string& filename)
{
  std::vector<double> contour_closure_distances;
  compute_contour_closure_distances(
    contour_segments, contours, contour_closure_distances);

  write_vector(contour_closure_distances, filename, ", ");
}
