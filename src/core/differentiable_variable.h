// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

inline size_t
generate_local_variable_matrix_index(size_t row,
                                     size_t col,
                                     size_t dimension = 3)
{
  return dimension * row + col;
}

/// Build an independent variable with a given value and variable index.
///
/// @param[in] value: initial value of the variable
/// @param[in] variable_index: global index of the variable in the full system
/// @param[in] total_independent_variables: total number of variables in the
/// system
/// @return constructed differentiable variable
template<typename Var>
Var
generate_independent_variable(double value,
                              int variable_index,
                              int total_independent_variables)
{
  DiffScalarBase::setVariableCount(total_independent_variables);
  return Var(variable_index, value);
}

/// Build a vector of independent variables with a given initial value and
/// contiguous variable indices from some starting index.
///
/// @param[in] value_vector: initial values of the variables
/// @param[out] variable_vector: constructed differentiable variables
/// @param[in] start_variable_index: global index of the first variable in the
/// full system
/// @param[in] total_independent_variables: total number of variables in the
/// system
template<typename Var, typename VectorXv>
void
build_independent_variable_vector(const VectorXr& value_vector,
                                  VectorXv& variable_vector,
                                  int start_variable_index,
                                  int total_independent_variables)
{
  DiffScalarBase::setVariableCount(total_independent_variables);
  size_t vector_size = value_vector.size();
  variable_vector.resize(vector_size);
  for (size_t i = 0; i < vector_size; ++i) {
    variable_vector[i] = generate_independent_variable<Var>(
      value_vector[i], start_variable_index + i, total_independent_variables);
  }
}

/// Build a matrix of independent variables with a given initial value and
/// contiguous row-major variable indices from some starting index.
///
/// @param[in] value_matrix: initial values of the variables
/// @param[out] value_matrix: constructed differentiable variables
/// @param[in] start_variable_index: global index of the first variable in the
/// full system
/// @param[in] total_independent_variables: total number of variables in the
/// system
template<typename Var, typename MatrixXv>
void
build_independent_variable_matrix(const MatrixXr& value_matrix,
                                  MatrixXv& variable_matrix,
                                  int start_variable_index,
                                  int total_independent_variables)
{
  DiffScalarBase::setVariableCount(total_independent_variables);
  size_t rows = value_matrix.rows();
  size_t cols = value_matrix.cols();
  variable_matrix.resize(rows, cols);
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      size_t local_index = generate_local_variable_matrix_index(i, j, cols);
      variable_matrix(i, j) =
        generate_independent_variable<Var>(value_matrix(i, j),
                                           start_variable_index + local_index,
                                           total_independent_variables);
    }
  }
}

/// Build a differentiable constant with a given value.
///
/// @param[in] value: value of the constant
/// @return constructed differentiable variable
template<typename Var>
Var
generate_constant_variable(double value)
{
  return Var(value);
}

/// Build a vector of differentiable constants with given values.
///
/// @param[in] value_vector: values of the constants
/// @param[out] constant_variable_vector: vector of constant variables
template<typename Var, typename VectorXv>
void
build_constant_variable_vector(const VectorXr& value_vector,
                               VectorXv& constant_variable_vector)
{
  size_t vector_size = value_vector.size();
  constant_variable_vector.resize(vector_size);
  for (size_t i = 0; i < vector_size; ++i) {
    constant_variable_vector[i] =
      generate_constant_variable<Var>(value_vector[i]);
  }
}

/// Build a matrix of differentiable constants with given values.
///
/// @param[in] value_matrix: values of the constants
/// @param[out] constant_variable_matrix: matrix of constant variables
template<typename Var, typename MatrixXv>
void
build_constant_variable_matrix(const MatrixXr& value_matrix,
                               MatrixXv& constant_variable_matrix)
{
  size_t rows = value_matrix.rows();
  size_t cols = value_matrix.cols();
  constant_variable_matrix.resize(rows, cols);
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      constant_variable_matrix(i, j) =
        generate_constant_variable<Var>(value_matrix(i, j));
    }
  }
}

/// Extract the value of a differentiable variable.
///
/// @param[in] variable: differentiable variable
/// @return value of the variable
template<typename Var>
double
compute_variable_value(const Var& variable)
{
  return variable.getValue();
}

/// Extract the gradient of a differentiable variable with respect to the
/// independent variables.
///
/// @param[in] variable: differentiable variable
/// @return gradient of the variable
template<typename Gradient, typename Var>
void
compute_variable_gradient(const Var& variable, Gradient& gradient)
{
  gradient = variable.getGradient();
}

/// Extract the hessian of a differentiable variable with respect to the
/// independent variables.
///
/// @param[in] variable: differentiable variable
/// @return hessian of the variable
template<typename Hessian, typename Var>
void
compute_variable_hessian(const Var& variable, Hessian& hessian)
{
  hessian = variable.getHessian();
}

/// Extract the values of a vector of differentiable variables.
///
/// @param[in] variable_vector: vector of differentiable variables
/// @param[out] values_vector: vector of the values of the variables
template<typename Var, typename VectorXv>
void
extract_variable_vector_values(const VectorXv& variable_vector,
                               VectorXr& values_vector)
{
  values_vector.resize(variable_vector.size());
  for (size_t i = 0; i < variable_vector.size(); ++i) {
    values_vector[i] = compute_variable_value<Var>(variable_vector[i]);
  }
}

/// Extract the values of a vector of differentiable variables.
///
/// @param[in] variable_vector: vector of differentiable variables
/// @return vector of the values of the variables
template<typename Var, typename VectorXv>
VectorXr
extract_variable_vector_values(const VectorXv& variable_vector)
{
  VectorXr values_vector;
  extract_variable_vector_values<Var, VectorXv>(variable_vector, values_vector);
  return values_vector;
}

/// Extract the values of a matrix of differentiable variables.
///
/// @param[in] variable_matrix: matrix of differentiable variables
/// @param[out] values_matrix: matrix of the values of the variables
template<typename Var, typename MatrixXv>
void
extract_variable_matrix_values(const MatrixXv& variable_matrix,
                               MatrixXr& values_matrix)
{
  values_matrix.resize(variable_matrix.rows(), variable_matrix.cols());
  for (size_t i = 0; i < variable_matrix.rows(); ++i) {
    for (size_t j = 0; j < variable_matrix.cols(); ++j) {
      values_matrix(i, j) = compute_variable_value<Var>(variable_matrix(i, j));
    }
  }
}

/// Extract the values of a matrix of differentiable variables.
///
/// @param[in] variable_matrix: matrix of differentiable variables
/// @return matrix of the values of the variables
template<typename Var, typename MatrixXv>
MatrixXr
extract_variable_matrix_values(const MatrixXv& variable_matrix)
{
  MatrixXr values_matrix;
  extract_variable_matrix_values<Var, MatrixXv>(variable_matrix, values_matrix);
  return values_matrix;
}

/// Determine if a vector of differentiable variables contains NaN.
///
/// @param[in] variable_vector: vector of differentiable variables
/// @return true iff the vector contains NaN
template<typename Var, typename VectorXv>
bool
vector_contains_nan(const VectorXv& variable_vector)
{
  for (Eigen::Index i = 0; i < variable_vector.size(); ++i) {
    if (std::isnan(compute_variable_value(variable_vector[i])))
      return true;
  }

  return false;
}

/// Determine if a matrix of differentiable variables contains NaN.
///
/// @param[in] variable_matrix: matrix of differentiable variables
/// @return true iff the matrix contains NaN
template<typename Var, typename MatrixXv>
bool
matrix_contains_nan(const MatrixXv& variable_matrix)
{
  for (Eigen::Index i = 0; i < variable_matrix.rows(); ++i) {
    for (Eigen::Index j = 0; j < variable_matrix.cols(); ++j) {
      if (std::isnan(compute_variable_value<Var>(variable_matrix(i, j))))
        return true;
    }
  }

  return false;
}

/// Compute the square of a variable
///
/// @param[in] x: differentiable variables
/// @return square of the variable
template<typename Var>
Var
variable_square(const Var& x)
{
  return x * x;
}

/// Compute the square norm of a variable vector.
///
/// @param[in] variable_vector: vector of differentiable variables
/// @return square norm of the variable
template<typename Var, typename VectorXv>
Var
variable_square_norm(const VectorXv& variable_vector)
{
  Var norm(0);
  size_t num_entries = variable_vector.size();
  for (size_t i = 0; i < num_entries; ++i) {
    norm += variable_square<Var>(variable_vector[i]);
  }

  return norm;
}
