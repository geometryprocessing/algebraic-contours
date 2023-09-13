// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

// #define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_OFF

#include "autodiff.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "spdlog/fmt/ostr.h"
#include "spdlog/spdlog.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <unsupported/Eigen/AutoDiff>

#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>

#include <igl/combine.h>
#include <igl/facet_components.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/remove_unreferenced.h>
#include <igl/vertex_components.h>

// Global epsilons
// const double MAX_PRECISION = 1e-8;
extern double FLOAT_EQUAL_PRECISION;      // Epsilon for default float
extern double ADJACENT_CONTOUR_PRECISION; // Epsilon for chaining contours
extern double
  PLANAR_BOUNDING_BOX_PRECISION; // Epsilon for curve-curve bounding box padding
extern double
  FIND_INTERSECTIONS_BEZIER_CLIPPING_PRECISION; // Epsilon for Bezier clipping
                                                // intersections

extern int DISCRETIZATION_LEVEL; // Spline surface discretization level
const int HASH_TABLE_SIZE = 70; // Size of spline surface hash table

// Real number representations
typedef Eigen::VectorXd VectorXr;
typedef Eigen::RowVectorXd OneFormXr;
typedef Eigen::Matrix<double, 1, 2> PlanarPoint;
typedef Eigen::Matrix<double, 1, 3> SpatialVector;
typedef Eigen::MatrixXd MatrixXr;
typedef Eigen::Matrix<double, 2, 3> Matrix2x3r;
typedef Eigen::Matrix<double, 2, 2> Matrix2x2r;
typedef Eigen::Matrix<double, 3, 2> Matrix3x2r;
typedef Eigen::Matrix<double, 3, 3> Matrix3x3r;
typedef std::array<int, 2> Edge;

// Bivariate quadratic matrices
typedef Eigen::Matrix<double, 6, 1> Vector6r;
typedef Eigen::Matrix<double, 6, 3> Matrix6x3r;
typedef Eigen::Matrix<double, 6, Eigen::Dynamic> Matrix6xXr;

// Typedefs from ConTess
typedef double real_t;
using Color = Eigen::Vector4f;

// Variable representation
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> DynamicGradient;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DynamicHessian;

// Six Split
typedef Eigen::Matrix<double, 27, 1> SixSplitGradient;
typedef Eigen::Matrix<double, 27, 27> SixSplitHessian;

// Twelve Split
typedef Eigen::Matrix<double, 36, 1> TwelveSplitGradient;
typedef Eigen::Matrix<double, 36, 36> TwelveSplitHessian;

// Colors
const Eigen::Matrix<double, 3, 1> MINT_GREEN(0.170, 0.673, 0.292);
const Eigen::Matrix<double, 3, 1> SKY_BLUE(0.297, 0.586, 0.758);
const Eigen::Matrix<double, 3, 1> OFF_WHITE(0.896, 0.932, 0.997);
const Eigen::Matrix<double, 3, 1> GOLD_YELLOW(0.670, 0.673, 0.292);

// Algebraic constrained values
const int MAX_PATCH_RAY_INTERSECTIONS = 4;

// ***********************
// Floating point equality
// ***********************

/// @brief  Check if some floating point value is numerically zero.
///
/// @param[in] x: value to compare with zero
/// @param[in] eps: threshold for equality
/// @return true iff x is below 1e-10
inline bool
float_equal_zero(double x, double eps = FLOAT_EQUAL_PRECISION)
{
  return abs(x) < eps;
}

/// @brief Check if two floating point values are numerically equal
///
/// @param[in] x: first value to compare
/// @param[in] y: second value to compare
/// @param[in] eps: threshold for equality
/// @return true iff x - y is numerically zero
inline bool
float_equal(double x, double y, double eps = FLOAT_EQUAL_PRECISION)
{
  return float_equal_zero(x - y, eps);
}

/// @brief Check if two row vectors of floating point values are numerically
/// equal
///
/// @param[in] v: first vector of values to compare
/// @param[in] w: second vector of values to compare
/// @param[in] eps: threshold for equality
/// @return true iff v - w is numerically the zero vector
template<int dimension>
bool
vector_equal(const Eigen::Matrix<double, 1, dimension>& v,
             const Eigen::Matrix<double, 1, dimension>& w,
             double eps = FLOAT_EQUAL_PRECISION)
{
  // Two vectors of different sizes are not equal
  if (v.size() != w.size())
    return false;

  // Check each component
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    if (!float_equal(v[i], w[i], eps))
      return false;
  }

  return true;
}

/// @brief Check if two column vectors of floating point values are numerically
/// equal
///
/// @param[in] v: first vector of values to compare
/// @param[in] w: second vector of values to compare
/// @param[in] eps: threshold for equality
/// @return true iff v - w is numerically the zero vector
template<int num_coeffs>
bool
column_vector_equal(const Eigen::Matrix<double, num_coeffs, 1>& v,
                    const Eigen::Matrix<double, num_coeffs, 1>& w,
                    double eps = FLOAT_EQUAL_PRECISION)
{
  // Two vectors of different sizes are not equal
  if (v.size() != w.size())
    return false;

  // Check each component
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    if (!float_equal(v[i], w[i], eps))
      return false;
  }

  return true;
}

/// @brief Check if two matrices of floating point values are numerically equal
///
/// @param[in] A: first matrix of values to compare
/// @param[in] B: second matrix of values to compare
/// @param[in] eps: threshold for equality
/// @return true iff A - B is numerically the zero matrix
inline bool
matrix_equal(const MatrixXr& A,
             const MatrixXr& B,
             double eps = FLOAT_EQUAL_PRECISION)
{
  // Two matrices of different sizes are not equal
  if (A.rows() != B.rows())
    return false;
  if (A.cols() != B.cols())
    return false;

  // Check each component
  for (Eigen::Index i = 0; i < A.rows(); ++i) {
    for (Eigen::Index j = 0; j < A.cols(); ++j) {
      if (!float_equal(A(i, j), B(i, j), eps))
        return false;
    }
  }

  return true;
}

// ******
// Viewer
// ******

/// @brief View a mesh in polyscope
///
/// @param[in] V: mesh vertices
/// @param[in] F: mesh faces
inline void
view_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
  polyscope::init();
  polyscope::registerSurfaceMesh("surface", V, F);
  polyscope::show();
  polyscope::removeAllStructures();
}

/// @brief View a mesh with uv coordinates in polyscope
///
/// @param[in] V: mesh vertices
/// @param[in] F: mesh faces
/// @param[in] V: mesh uv coordinates
inline void
view_parametrized_mesh(const Eigen::MatrixXd& V,
                       const Eigen::MatrixXi& F,
                       const Eigen::MatrixXd& uv)
{
  polyscope::init();
  polyscope::registerSurfaceMesh("surface", V, F);
  polyscope::getSurfaceMesh("surface")->addVertexParameterizationQuantity(
    "parameterization", uv);
  polyscope::show();
}

inline void
screenshot_mesh(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const std::string& filename,
  SpatialVector camera_position = SpatialVector(0, 0, 0),
  SpatialVector camera_target = SpatialVector(0, 0, 2),
  bool use_orthographic = false)
{
  polyscope::init();
  polyscope::registerSurfaceMesh("surface", V, F)
    ->setEdgeWidth(1)
    ->setSurfaceColor(glm::vec3(0.670, 0.673, 0.292));
  //polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
  glm::vec3 glm_camera_position = { camera_position[0],
                                    camera_position[1],
                                    camera_position[2] };
  glm::vec3 glm_camera_target = { camera_target[0],
                                  camera_target[1],
                                  camera_target[2] };
  polyscope::view::lookAt(glm_camera_position, glm_camera_target);
  if (use_orthographic) {
    polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
  }
  else {
    polyscope::view::projectionMode = polyscope::ProjectionMode::Perspective;
  }
  polyscope::screenshot(filename);
  polyscope::removeAllStructures();
}

// ****************
// Basic arithmetic
// ****************

/// @brief Check the sign of a value x
///
/// Returns -1 for negative, 0 for zero, and 1 for positive
///
/// @param[in] x: value to check sign of
/// @return sign of x
inline double
sgn(double x)
{
  return (x > 0) - (x < 0);
}

/// @brief Compute x^p using iteration.
///
/// @param[in] x: base
/// @param[in] p: exponent
/// @return result of x^p
inline double
power(double x, int p)
{
  double xp = 1.0;
  for (int i = 0; i < p; ++i) {
    xp *= x;
  }

  return xp;
}

/// @brief Compute the discriminant of the quadratic ax^2 + bx + c.
///
/// @param[in] a: x^2 coefficient
/// @param[in] b: x coefficient
/// @param[in] c: constant coefficient
/// @return discriminate b^2 - 4ac
inline double
compute_discriminant(double a, double b, double c)
{
  return (b * b) - (4.0 * a * c);
}

// ********************
// Basic linear algebra
// ********************

/// @brief  Compute the dot product of two vectors of arbitrary scalars.
///
/// @tparam Scalar: scalar field (must support addition and multiplication)
/// @param[in] v: first vector to dot product
/// @param[in] w: second vector to dot product
/// @return dot product of v and w
template<typename Scalar, int dimension>
Scalar
dot_product(Eigen::Matrix<Scalar, dimension, 1> v,
            Eigen::Matrix<Scalar, dimension, 1> w)
{
  return v.dot(w);
}

/// @brief  Compute the cross product of two vectors of arbitrary scalars.
///
/// @tparam Scalar: scalar field (must support addition and multiplication)
/// @param[in] v: first vector to cross product
/// @param[in] w: second vector to cross product
/// @return cross product v x w
template<typename Scalar>
Eigen::Matrix<Scalar, 3, 1>
cross_product(Eigen::Matrix<Scalar, 3, 1> v, Eigen::Matrix<Scalar, 3, 1> w)
{
  Eigen::Matrix<Scalar, 3, 1> n(3);
  n(0) = v(1) * w(2) - v(2) * w(1);
  n(1) = -(v(0) * w(2) - v(2) * w(0));
  n(2) = v(0) * w(1) - v(1) * w(0);
  return n;
}

/// @brief  Compute the triple product of three vectors.
///
/// @tparam Scalar: scalar field (must support addition and multiplication)
/// @param[in] u: first vector to triple product
/// @param[in] v: second vector to triple product
/// @param[in] w: third vector to triple product
/// @return triple product u * (v x w)
template<typename Scalar>
Scalar
triple_product(Eigen::Matrix<Scalar, 3, 1> u,
               Eigen::Matrix<Scalar, 3, 1> v,
               Eigen::Matrix<Scalar, 3, 1> w)
{
  Eigen::Matrix<Scalar, 3, 1> vxw = cross_product<Scalar>(v, w);
  return dot_product<Scalar, 3>(u, vxw);
}

/// @brief Normalize a vector to unit length
///
/// For the zero vector, return the zero vector.
///
/// @tparam Scalar: scalar field for the vector (must support all field
/// operations)
/// @param[in] v: vector to normalize
/// @return normalized vector v / ||v||
template<typename Scalar, int dimension>
Eigen::Matrix<Scalar, dimension, 1>
normalize(Eigen::Matrix<Scalar, dimension, 1> v)
{
  // Do nothing for empty vectors
  if (v.size() == 0)
    return v;

  // Get the norma and check if it's zero
  Scalar v_norm = v.norm();
  if (float_equal_zero(v.norm()))
    return v;

  return v / v_norm;
}

/// @brief  Generate the elementary basis vector for Euclidean space of given
/// dimension with all entries zero except a 1 in position i
///
/// @tparam dimension: dimension of the vector space
/// @param[in] i: index for the nonzero entry
/// @return basis vector ei
template<int dimension>
Eigen::Matrix<double, 1, dimension>
elementary_basis_vector(int i)
{
  // Initialize zero vector
  Eigen::Matrix<double, 1, dimension> basis_vector;
  basis_vector.setZero();

  // Return zero if i exceeds what is possible for the dimension
  if ((i >= dimension) || (i < 0))
    return basis_vector;

  // Set ith coordinate to 1
  basis_vector(i) = 1.0;
  return basis_vector;
}

/// @brief  Reflect a vector in the plane across the x-axis.
///
/// @param[in] vector: vector to reflect
/// @return reflected vector
inline PlanarPoint
reflect_across_x_axis(const PlanarPoint& vector)
{
  PlanarPoint reflected_vector;
  reflected_vector[0] = vector[0];
  reflected_vector[1] = -vector[1];
  return reflected_vector;
}

/// @brief  Rotate a vector in three-dimensional space a given angle around an
/// axis.
///
/// @param[in] vector: vector to rotate
/// @param[in] axis_of_rotation: vector giving the rotation axis
/// @param[in] angle: angle to rotate the vector by
/// @param[out] rotated_vector: rotated vector
inline void
rotate_vector(const SpatialVector& vector,
              const SpatialVector& axis_of_rotation,
              double angle,
              SpatialVector& rotated_vector)
{
  // Do nothing for the zero axis of rotation
  if (float_equal_zero(axis_of_rotation.norm())) {
    rotated_vector = vector;
    return;
  }

  // Rewrite variables in more compact notation and ensure the axis is normal
  SpatialVector const& v = vector;
  SpatialVector const& k = axis_of_rotation / axis_of_rotation.norm();
  double theta = angle;

  // Build rotated vector w from three components along v, k,
  SpatialVector n = cross_product<double>(k, v);
  SpatialVector w_v = std::cos(theta) * v;
  SpatialVector w_k = dot_product<double, 3>(k, v) * (1 - std::cos(theta)) * k;
  SpatialVector w_n = std::sin(theta) * n;
  rotated_vector = w_v + w_k + w_n;
}

/// @brief  Rotate a vector in three-dimensional space a given angle around an
/// axis.
///
/// @param[in] vector: vector to rotate
/// @param[in] axis_of_rotation: vector giving the rotation axis
/// @param[in] angle: angle to rotate the vector by
/// @return rotated vector
inline SpatialVector
rotate_vector(const SpatialVector& vector,
              const SpatialVector& axis_of_rotation,
              double angle)
{
  SpatialVector rotated_vector;
  rotate_vector(vector, axis_of_rotation, angle, rotated_vector);
  return rotated_vector;
}

/// @brief Project a vector to the plane defined by some normal.
///
/// @tparam Scalar: scalar for the vector field
/// @param[in] vector: vector to project
/// @param[in] plane_normal: normal to the projection plane
/// @param[out] projected_vector: projected vector
template<typename Scalar, int dimension>
Eigen::Matrix<Scalar, dimension, 1>
project_vector_to_plane(const Eigen::Matrix<Scalar, dimension, 1>& vector,
                        const Eigen::Matrix<Scalar, dimension, 1>& plane_normal)
{
  Scalar vector_normal_component =
    dot_product<double, dimension>(vector, plane_normal);
  Scalar normal_length_sq =
    dot_product<double, dimension>(plane_normal, plane_normal);

  // Do nothing for the zero plane normal
  if (normal_length_sq == Scalar(0.0))
    return vector;

  // Project vector
  return vector - (vector_normal_component / normal_length_sq) * plane_normal;
}

// **************
// Vector methods
// **************

template<int dimension>
double
vector_min(const Eigen::Matrix<double, 1, dimension>& v)
{
  if (v.size() == 0)
    return 0;

  double v_min = v[0];
  for (Eigen::Index i = 1; i < v.size(); ++i) {
    if (v[i] < v_min)
      v_min = v[i];
  }

  return v_min;
}

template<int dimension>
double
vector_max(const Eigen::Matrix<double, 1, dimension>& v)
{
  if (v.size() == 0)
    return 0;

  double v_max = v[0];
  for (Eigen::Index i = 1; i < v.size(); ++i) {
    if (v[i] > v_max)
      v_max = v[i];
  }

  return v_max;
}

template<int dimension>
double
column_vector_min(const Eigen::Matrix<double, dimension, 1>& v)
{
  if (v.size() == 0)
    return 0;

  double v_min = v[0];
  for (Eigen::Index i = 1; i < v.size(); ++i) {
    if (v[i] < v_min)
      v_min = v[i];
  }

  return v_min;
}

template<int dimension>
double
column_vector_max(const Eigen::Matrix<double, dimension, 1>& v)
{
  if (v.size() == 0)
    return 0;

  double v_max = v[0];
  for (Eigen::Index i = 1; i < v.size(); ++i) {
    if (v[i] > v_max)
      v_max = v[i];
  }

  return v_max;
}

template<typename T>
bool
vector_contains(const std::vector<T>& vec, const T item)
{
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i] == item)
      return true;
  }

  return false;
}

template<typename Index>
void
convert_index_vector_to_boolean_array(const std::vector<Index>& index_vector,
                                      Index num_indices,
                                      std::vector<bool>& boolean_array)
{
  boolean_array.resize(num_indices, false);
  for (size_t i = 0; i < index_vector.size(); ++i) {
    boolean_array[index_vector[i]] = true;
  }
}

/// @brief From a boolean array, build a vector of the indices that are true.
///
/// @param[in] boolean_array: array of boolean values
/// @param[out] index_vector: indices where the array is true
inline void
convert_boolean_array_to_index_vector(const std::vector<bool>& boolean_array,
                                      std::vector<size_t>& index_vector)
{
  size_t num_indices = boolean_array.size();
  index_vector.clear();
  index_vector.reserve(num_indices);
  for (size_t i = 0; i < num_indices; ++i) {
    if (boolean_array[i]) {
      index_vector.push_back(i);
    }
  }
}

template<typename Index>
void
index_vector_complement(const std::vector<Index>& index_vector,
                        Index num_indices,
                        std::vector<Index>& complement_vector)
{
  // Build index boolean array
  std::vector<bool> boolean_array;
  convert_index_vector_to_boolean_array(
    index_vector, num_indices, boolean_array);

  // Build complement
  complement_vector.clear();
  complement_vector.reserve(num_indices - index_vector.size());
  for (Index i = 0; i < num_indices; ++i) {
    if (!boolean_array[i]) {
      complement_vector.push_back(i);
    }
  }
}

inline void
convert_signed_vector_to_unsigned(const std::vector<int>& signed_vector,
                                  std::vector<size_t>& unsigned_vector)
{
  size_t vector_size = signed_vector.size();
  unsigned_vector.resize(vector_size);
  for (size_t i = 0; i < vector_size; ++i) {
    unsigned_vector[i] = signed_vector[i];
  }
}

inline void
convert_unsigned_vector_to_signed(const std::vector<size_t>& unsigned_vector,
                                  std::vector<int>& signed_vector)
{
  size_t vector_size = unsigned_vector.size();
  signed_vector.resize(vector_size);
  for (size_t i = 0; i < vector_size; ++i) {
    signed_vector[i] = unsigned_vector[i];
  }
}

template<typename T, typename Index>
void
remove_vector_values(const std::vector<Index>& indices_to_remove,
                     const std::vector<T>& vec,
                     std::vector<T>& subvec)
{
  // Remove faces adjacent to cones
  std::vector<Index> indices_to_keep;
  index_vector_complement<Index>(
    indices_to_remove, vec.size(), indices_to_keep);
  subvec.resize(indices_to_keep.size());
  for (size_t i = 0; i < indices_to_keep.size(); ++i) {
    subvec[i] = vec[indices_to_keep[i]];
  }
}

template<typename T>
void
copy_to_planar_point(const std::vector<T>& input_vector,
                     std::vector<PlanarPoint>& output_vector)
{
  output_vector.resize(input_vector.size());
  for (size_t i = 0; i < input_vector.size(); ++i) {
    output_vector[i] << input_vector[i][0], input_vector[i][1];
  }
}

template<typename T>
void
copy_to_spatial_vector(const std::vector<T>& input_vector,
                       std::vector<SpatialVector>& output_vector)
{
  output_vector.resize(input_vector.size());
  for (size_t i = 0; i < input_vector.size(); ++i) {
    output_vector[i] << input_vector[i][0], input_vector[i][1],
      input_vector[i][2];
  }
}

template<typename T>
inline std::string
formatted_vector(const std::vector<T>& vec, std::string delim = "\n")
{
  std::stringstream vector_string;
  for (size_t i = 0; i < vec.size(); ++i) {
    vector_string << vec[i] << delim;
  }

  return vector_string.str();
}

/// @brief Write vector to file in csv file.
///
/// @tparam T: vector data type
/// @param[in] vec: vector to write
/// @param[in] filename: filename to write the vector to
/// @param[in] delim: deliminator between vector entries
template<typename T>
void
write_vector(const std::vector<T>& vec,
             const std::string& filename,
             std::string delim = "\n")
{
  std::ofstream output_file;
  output_file.open(filename, std::ios::out | std::ios::trunc);
  for (size_t i = 0; i < vec.size(); ++i) {
    output_file << vec[i] << delim;
  }
  output_file.close();
}

/// @brief Write floating point vector to file in csv file with given precision.
///
/// @tparam T: vector data type
/// @param[in] vec: vector to write
/// @param[in] filename: filename to write the vector to
/// @param[in] delim: deliminator between vector entries
/// @param[in] prec: precision for the output
inline void
write_float_vector(const std::vector<double>& vec,
                   const std::string& filename,
                   std::string delim = "\n",
                   int prec = 17)
{
  std::ofstream output_file;
  output_file.open(filename, std::ios::out | std::ios::trunc);
  for (size_t i = 0; i < vec.size(); ++i) {
    output_file << std::setprecision(prec) << vec[i] << delim;
  }
  output_file.close();
}

template<typename T>
inline void
append(std::vector<T>& vec, std::vector<T>& vec_to_add)
{
  vec.insert(vec.end(), vec_to_add.begin(), vec_to_add.end());
}

template<typename T>
inline size_t
nested_vector_size(std::vector<std::vector<T>> v)
{
  size_t count = 0;

  for (size_t i = 0; i < v.size(); ++i) {
    count += v[i].size();
  }

  return count;
}

template<typename T>
inline MatrixXr
convert_nested_vector_to_matrix(const std::vector<T>& vec)
{
  size_t n = vec.size();
  if (n <= 0)
    return MatrixXr();
  MatrixXr matrix(vec.size(), vec[0].size());
  for (size_t i = 0; i < n; ++i) {
    size_t inner_vec_size = vec[i].size();
    for (size_t j = 0; j < inner_vec_size; ++j) {
      matrix(i, j) = vec[i][j];
    }
  }

  return matrix;
}

inline Eigen::MatrixXd
convert_nested_vector_to_matrix(const std::vector<std::vector<double>>& vec)
{
  size_t n = vec.size();
  if (n <= 0)
    return Eigen::MatrixXd();
  size_t m = vec[0].size();
  Eigen::MatrixXd matrix(vec.size(), vec[0].size());
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      matrix(i, j) = vec[i][j];
    }
  }

  return matrix;
}

inline Eigen::MatrixXi
convert_nested_vector_to_matrix(const std::vector<std::vector<int>>& vec)
{
  size_t n = vec.size();
  if (n <= 0)
    return Eigen::MatrixXi();
  size_t m = vec[0].size();
  Eigen::MatrixXi matrix(vec.size(), vec[0].size());
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      matrix(i, j) = vec[i][j];
    }
  }

  return matrix;
}

// **************
// Matrix methods
// **************

template<typename T>
void
append_matrix(T& matrix, T& matrix_to_add)
{
  int i = matrix.rows();
  int j = matrix.cols();
  int p = matrix_to_add.rows();
  int q = matrix_to_add.cols();

  if (j == 0) {
    matrix = matrix_to_add;
  } else {
    assert(j == q);
    matrix.conservativeResize(i + p, j);
    matrix.block(i, 0, p, q) = matrix_to_add;
  }
}

template<typename Derived>
void
flatten_matrix_by_row(const Eigen::EigenBase<Derived>& mat,
                      std::vector<double>& vec)
{
  size_t row_length = mat.cols();
  size_t num_entries = mat.rows() * mat.cols();
  vec.resize(num_entries);

  for (Eigen::Index i = 0; i < mat.rows(); ++i) {
    for (Eigen::Index j = 0; j < mat.cols(); ++j) {
      size_t flat_index = row_length * i + j;
      vec[flat_index] = mat(i, j);
    }
  }
}

// Read a 4x4 matrix from file
///
/// @param[in] filename: file with matrix to read
/// @param[out] vec: vector from file
inline void read_camera_matrix(
  const std::string &filename,
  Eigen::Matrix<double, 4, 4>& mat
) {
  // Open file
  std::ifstream input_file(filename);
  if (!input_file) return;

  // Read file
  std::string line;
  int row = 0;
  while (std::getline(input_file, line))
  {
    std::istringstream iss(line);
    std::string cell;
    int col = 0;
    while (std::getline(iss, cell, ','))
    {
      mat(row, col) = std::stod(cell);
      ++col;
    }
    ++row;
  }

  // Close file
  input_file.close();
}

// ****************
// Pythonic methods
// ****************

inline void
generate_linspace(double t_0, double t_1, size_t num_points, VectorXr& linspace)
{
  linspace.resize(num_points);
  if (num_points < 2)
    return;

  for (size_t i = 0; i < num_points; ++i) {
    double delta_t = (t_1 - t_0) / static_cast<double>(num_points - 1);
    linspace[i] = t_0 + delta_t * i;
  }
}

template<typename T>
void
arange(size_t size, std::vector<T>& arange_vec)
{
  arange_vec.resize(size);
  for (size_t i = 0; i < size; ++i) {
    arange_vec[i] = T(i);
  }
}

// *******************
// Basic mesh topology
// *******************

// Return true iff the face contains the given vertex
inline bool
contains_vertex(Eigen::VectorXi face, int vertex_index)
{
  for (Eigen::Index i = 0; i < face.size(); ++i) {
    if (face[i] == vertex_index)
      return true;
  }

  return false;
}

inline int
find_face_vertex_index(Eigen::VectorXi face, int vertex_index)
{
  for (Eigen::Index i = 0; i < face.size(); ++i) {
    if (face[i] == vertex_index)
      return i;
  }

  return -1;
}

/// @brief Check if F describes a manifold mesh with a single component
///
/// @param[in] F: mesh faces
/// @return true iff the mesh is manifold
inline bool
is_manifold(const Eigen::MatrixXi& F)
{
  // Check edge manifold condition
  if (!igl::is_edge_manifold(F)) {
    spdlog::error("Mesh is not edge manifold");
    return false;
  }

  // Check vertex manifold condition
  Eigen::VectorXi invalid_vertices;
  if (!igl::is_vertex_manifold(F, invalid_vertices)) {
    spdlog::error("Mesh is not edge manifold");
    return false;
  }

  // Manifold otherwise
  return true;
}

// *******************
// Basic mesh geometry
// *******************

/// @brief Compute the area of a triangle from the edge lengths.
///
/// @param[in] l0: first edge length
/// @param[in] l1: second edge length
/// @param[in] l2: third edge length
/// @return area of the triangle
inline double
area_from_length(double l0, double l1, double l2)
{
  // Return the area (or zero if there is a triangle inequality violation)
  double s = 0.5 * (l0 + l1 + l2); // semi-perimeter
  double area = std::sqrt(std::max(s * (s - l0) * (s - l1) * (s - l2), 0.0));
  assert(!std::isnan(area));
  return area;
}

/// @brief Compute the area of a triangle from the vertex positions
///
/// @param[in] p0: first vertex position
/// @param[in] p1: second vertex position
/// @param[in] p2: third vertex position
/// @return triangle area
template<int dimension>
double
area_from_positions(const Eigen::Matrix<double, 1, dimension>& p0,
                    const Eigen::Matrix<double, 1, dimension>& p1,
                    const Eigen::Matrix<double, 1, dimension>& p2)
{
  double l0 = (p2 - p1).norm();
  ;
  double l1 = (p0 - p2).norm();
  ;
  double l2 = (p1 - p0).norm();
  ;
  return area_from_length(l0, l1, l2);
}

/// @brief Compute the angle of a triangle corner with given edge lengths
///
/// @param[in] edge_length_opposite_corner: length of the edge opposite the
/// corner
/// @param[in] first_adjacent_edge_length: length of one of the edges adjacent
/// to the corner
/// @param[in] first_adjacent_edge_length: length of the other edge adjacent to
/// the corner
/// @return angle of the corner
inline double
angle_from_length(double edge_length_opposite_corner,
                  double first_adjacent_edge_length,
                  double second_adjacent_edge_length)
{
  // Rename variables for readability
  double l0 = edge_length_opposite_corner;
  double l1 = first_adjacent_edge_length;
  double l2 = second_adjacent_edge_length;

  // Compute the angle
  // FIXME Avoid potential division by 0
  double Ijk = (-l0 * l0 + l1 * l1 + l2 * l2);
  return acos(std::min<double>(std::max(Ijk / (2.0 * l1 * l2), -1.0), 1.0));
}

/// @brief Compute the angle of a triangle corner with given positions
///
/// @param[in] angle_corner_position: position of the corner to compute the
/// angle for
/// @param[in] second_corner_position: position of one of the other two corners
/// of the triangle
/// @param[in] third_corner_position: position of the final corner of the
/// triangle
/// @return angle of the corner
template<int dimension>
double
angle_from_positions(
  const Eigen::Matrix<double, 1, dimension>& angle_corner_position,
  const Eigen::Matrix<double, 1, dimension>& second_corner_position,
  const Eigen::Matrix<double, 1, dimension>& third_corner_position)
{
  double l0 = (third_corner_position - second_corner_position).norm();
  double l1 = (second_corner_position - angle_corner_position).norm();
  double l2 = (third_corner_position - angle_corner_position).norm();
  return angle_from_length(l0, l1, l2);
}

/// @brief Map [t_min_0, t_max_0] -> [t_min_1, t_max_1] with the unique linear
/// isomorphism
///
/// @param[in] t_min_0: minimum of the domain interval
/// @param[in] t_max_0: maximum of the domain interval
/// @param[in] t_min_1: minimum of the image interval
/// @param[in] t_max_1: maximum of the domain interval
/// @param[in] t_0: point in the domain interval
/// @return mapped point in the iamge
inline double
interval_lerp(double t_min_0,
              double t_max_0,
              double t_min_1,
              double t_max_1,
              double t_0)
{
  // Return the midpoint of the image if the input domain is trivial
  if (float_equal(t_min_0, t_max_0)) {
    return 0.5 * (t_min_1 + t_max_1);
  }

  // Perform the interpolation
  double r_0 = t_max_0 - t_min_0;
  double r_1 = t_max_1 - t_min_1;
  double t_1 = t_min_1 + (r_1 / r_0) * (t_0 - t_min_0);
  return t_1;
}

/// @brief Compute the bounding box for a matrix of points in R^n.
///
/// The points are assumed to be the rows of the points matrix.
///
/// @param[in] points: points to compute the bounding box for
/// @param[out] min_point: point with minimum coordinates for the bounding box
/// @param[out] max_point: point with maximum coordinates for the bounding box
template<typename Vector, typename Matrix>
void
compute_point_cloud_bounding_box(const Matrix& points,
                                 Vector& min_point,
                                 Vector& max_point)
{
  Eigen::Index num_points = points.rows();
  Eigen::Index dimension = points.cols();
  if (num_points == 0)
    return;
  if (dimension == 0)
    return;

  // Get minimum and maximum coordinates for the points
  min_point = points.row(0);
  max_point = points.row(0);
  for (Eigen::Index pi = 0; pi < num_points; ++pi) {
    for (Eigen::Index j = 0; j < dimension; ++j) {
      min_point[j] = std::min<double>(min_point[j], points(pi, j));
      max_point[j] = std::max<double>(max_point[j], points(pi, j));
    }
  }
}

template<typename Index>
void
remove_mesh_faces(const Eigen::MatrixXd& V,
                  const Eigen::MatrixXi& F,
                  const std::vector<Index>& faces_to_remove,
                  Eigen::MatrixXd& V_submesh,
                  Eigen::MatrixXi& F_submesh)
{
  std::vector<Index> faces_to_keep;
  index_vector_complement<Index>(faces_to_remove, F.rows(), faces_to_keep);
  Eigen::MatrixXi F_unsimplified_submesh(faces_to_keep.size(), F.cols());
  F_unsimplified_submesh.resize(faces_to_keep.size(), 3);
  for (size_t i = 0; i < faces_to_keep.size(); ++i) {
    F_unsimplified_submesh.row(i) = F.row(faces_to_keep[i]);
  }

  // Remove unreferenced vertices and update face indices
  Eigen::MatrixXi I, J;
  igl::remove_unreferenced(V, F_unsimplified_submesh, V_submesh, F_submesh, I);
  SPDLOG_TRACE("Final mesh has {} faces and {} vertices",
               F_submesh.rows(),
               V_submesh.rows());
}

template<typename VertexIndex, typename FaceIndex>
void
remove_mesh_vertices(const Eigen::MatrixXd& V,
                     const Eigen::MatrixXi& F,
                     const std::vector<VertexIndex>& vertices_to_remove,
                     Eigen::MatrixXd& V_submesh,
                     Eigen::MatrixXi& F_submesh,
                     std::vector<FaceIndex>& faces_to_remove)
{
  SPDLOG_TRACE("Removing {} vertices from mesh with {} faces and {} vertices",
               vertices_to_remove.size(),
               F.rows(),
               V.rows());

  // Tag faces adjacent to the vertices to remove
  faces_to_remove.clear();
  for (Eigen::Index face_index = 0; face_index < F.rows(); ++face_index) {
    for (size_t i = 0; i < vertices_to_remove.size(); ++i) {
      if (contains_vertex(F.row(face_index), vertices_to_remove[i])) {
        faces_to_remove.push_back(face_index);
        break;
      }
    }
  }
  spdlog::trace("Removing {} faces", faces_to_remove.size());

  // Remove faces adjacent to cones
  remove_mesh_faces(V, F, faces_to_remove, V_submesh, F_submesh);
}

/// Given a mesh, separate it into lists of separate connected components.
///
/// @param[in] V: initial mesh vertices
/// @param[in] F: initial mesh faces
/// @param[out] V_components: list of component mesh vertices
/// @param[out] F_components: list of component mesh faces
/// @param[out] I_components: list of component maps from old to new vertices
/// @param[out] J_components: list of component maps from new to old vertices
inline void
separate_mesh_components(const Eigen::MatrixXd& V,
                        const Eigen::MatrixXi& F,
                        std::vector<Eigen::MatrixXd>& V_components,
                        std::vector<Eigen::MatrixXi>& F_components,
                        std::vector<Eigen::VectorXi>& I_components,
                        std::vector<Eigen::VectorXi>& J_components)
{
  // Get labels for mesh components
  Eigen::VectorXi component_ids;
  igl::facet_components(F, component_ids);
  int num_components = component_ids.maxCoeff() + 1;
  std::vector<int> counts(num_components, 0);
  int num_faces = F.rows();
  for (int i = 0; i < num_faces; ++i)
  {
    ++counts[component_ids[i]];
  }

  // Split mesh faces into components
  std::vector<Eigen::MatrixXi> F_split(num_components);
  std::vector<int> F_counts(num_components, 0);
  for (int i = 0; i < num_components; ++i)
  {
    F_split[i].resize(counts[i], 3);
  }
  for (int fi = 0; fi < num_faces; ++fi)
  {
    // Get component and current row for the component
    int ci = component_ids[fi];
    int ri = F_counts[ci];

    // Set face component
    F_split[ci].row(ri) = F.row(fi);

    // Increment face count
    ++F_counts[ci];
  }

  // Split mesh vertices into components and reindex face lists
  V_components.resize(num_components);
  F_components.resize(num_components);
  I_components.resize(num_components);
  J_components.resize(num_components);
  for (int i = 0; i < num_components; ++i)
  {
    igl::remove_unreferenced(
      V,
      F_split[i],
      V_components[i],
      F_components[i],
      I_components[i],
      J_components[i]
    );
  }
}

/// Given a lists of separate mesh components, combine them into a single mesh
///
/// @param[in] V_components: list of component mesh vertices
/// @param[in] F_components: list of component mesh faces
/// @param[out] V: combined mesh vertices
/// @param[out] F: combined mesh faces
inline void
combine_mesh_components(const std::vector<Eigen::MatrixXd>& V_components,
                        const std::vector<Eigen::MatrixXi>& F_components,
                        Eigen::MatrixXd& V,
                        Eigen::MatrixXi& F)
{
  // Just use the libigl implementation
  igl::combine(V_components, F_components, V, F);
}



// *********************
// Filepath manipulation
// *********************

inline std::string
join_path(const std::string& first_path, const std::string& second_path)
{
  if (first_path.back() == '/') {
    return first_path + second_path;
  } else {
    std::string path = first_path + "/";
    return path + second_path;
  }
}

// **********
// Nan checks
// **********

inline bool
matrix_contains_nan(const MatrixXr& mat)
{
  for (Eigen::Index i = 0; i < mat.rows(); ++i) {
    for (Eigen::Index j = 0; j < mat.cols(); ++j) {
      if (std::isnan(mat(i, j)))
        return true;
    }
  }

  return false;
}

template<int dimension>
bool
vector_contains_nan(const Eigen::Matrix<double, 1, dimension>& vec)
{
  for (Eigen::Index i = 0; i < vec.size(); ++i) {
    if (std::isnan(vec[i]))
      return true;
  }

  return false;
}

// TODO Move elsewhere

inline std::vector<Edge>
convert_polylines_to_edges(const std::vector<std::vector<int>>& polylines)
{
  std::vector<Edge> edges;
  edges.reserve(10 * polylines.size());

  for (size_t i = 0; i < polylines.size(); ++i) {
    for (size_t j = 1; j < polylines[i].size(); ++j) {
      Edge edge = { polylines[i][j - 1], polylines[i][j] };
      edges.push_back(edge);
    }
  }

  return edges;
}
