// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "quadratic_spline_surface_patch.h"

// ******************************
// Quadratic Spline Surface Patch
// ******************************

// Compute the surface mapping with normalized domain
void
compute_normalized_surface_mapping(
  const Matrix6x3r& surface_mapping_coeffs,
  const ConvexPolygon& domain,
  Matrix6x3r& normalized_surface_mapping_coeffs)
{
  // Normalize the surface coefficients
  MatrixXr domain_vertices = domain.get_vertices();
  Eigen::Matrix<double, 6, 6> change_of_basis_matrix;
  generate_quadratic_coordinate_domain_triangle_normalization_matrix<double>(
    domain_vertices.row(0),
    domain_vertices.row(1),
    domain_vertices.row(2),
    change_of_basis_matrix);
  normalized_surface_mapping_coeffs =
    change_of_basis_matrix * surface_mapping_coeffs;
}

void
compute_bezier_points(const Matrix6x3r& normalized_surface_mapping_coeffs,
                      Matrix6x3r& bezier_points)
{
  // Compute the bezier points
  Eigen::Matrix<double, 6, 6> monomial_to_bezier_matrix;
  generate_monomial_to_bezier_matrix(monomial_to_bezier_matrix);

  bezier_points = monomial_to_bezier_matrix * normalized_surface_mapping_coeffs;
}

// **************
// Public Methods
// **************

QuadraticSplineSurfacePatch::QuadraticSplineSurfacePatch()
{
  m_surface_mapping_coeffs.setZero();
  m_normal_mapping_coeffs.setZero();
  m_normalized_surface_mapping_coeffs.setZero();
  m_bezier_points.setZero();
  m_min_point.setZero();
  m_max_point.setZero();
  m_cone_index = -1;
}

QuadraticSplineSurfacePatch::QuadraticSplineSurfacePatch(
  const Matrix6x3r& surface_mapping_coeffs,
  const ConvexPolygon& domain)
  : m_surface_mapping_coeffs(surface_mapping_coeffs)
  , m_domain(domain)
{
  // Compute derived mapping information from the surface mapping and domain
  generate_quadratic_surface_normal_coeffs(surface_mapping_coeffs,
                                           m_normal_mapping_coeffs);
  compute_normalized_surface_mapping(
    surface_mapping_coeffs, domain, m_normalized_surface_mapping_coeffs);
  compute_bezier_points(m_normalized_surface_mapping_coeffs, m_bezier_points);
  compute_point_cloud_bounding_box(m_bezier_points, m_min_point, m_max_point);

  // Do not mark a cone by default
  m_cone_index = -1;
}

int
QuadraticSplineSurfacePatch::dimension() const
{
  return m_surface_mapping_coeffs.cols();
}

void
QuadraticSplineSurfacePatch::mark_cone(int cone_index)
{
  m_cone_index = cone_index;
}

bool
QuadraticSplineSurfacePatch::has_cone() const
{
  return ((m_cone_index >= 0) && (m_cone_index < 3));
}

int
QuadraticSplineSurfacePatch::get_cone() const
{
  return m_cone_index;
}

Matrix6x3r const&
QuadraticSplineSurfacePatch::get_surface_mapping() const
{
  return m_surface_mapping_coeffs;
}

Matrix6x3r const&
QuadraticSplineSurfacePatch::get_normal_mapping() const
{
  return m_normal_mapping_coeffs;
}

Matrix6x3r const&
QuadraticSplineSurfacePatch::get_normalized_surface_mapping() const
{
  return m_normalized_surface_mapping_coeffs;
}

Matrix6x3r const&
QuadraticSplineSurfacePatch::get_bezier_points() const
{
  return m_bezier_points;
}

ConvexPolygon const&
QuadraticSplineSurfacePatch::get_domain() const
{
  return m_domain;
}

void
QuadraticSplineSurfacePatch::get_patch_boundaries(
  std::array<RationalFunction<4, 3>, 3>& patch_boundaries) const
{
  // Get parametrized domain boundaries
  std::array<LineSegment, 3> domain_boundaries;
  get_domain().parametrize_patch_boundaries(domain_boundaries);

  // Lift the domain boundaries to the surface
  Matrix6x3r const& surface_mapping_coeffs = get_surface_mapping();
  for (size_t i = 0; i < domain_boundaries.size(); ++i) {
    domain_boundaries[i].pullback_quadratic_function<3>(surface_mapping_coeffs,
                                                        patch_boundaries[i]);
  }
}

void
QuadraticSplineSurfacePatch::normalize_patch_domain(
  QuadraticSplineSurfacePatch& normalized_spline_surface_patch) const
{
  // Generate the standard u + v <= 1 triangle
  Matrix3x2r normalized_domain_vertices;
  normalized_domain_vertices << 0, 0, 1, 0, 0, 1;
  ConvexPolygon normalized_domain(normalized_domain_vertices);

  // Build the normalized surface patch
  Matrix6x3r const& normalized_surface_mapping_coeffs =
    get_normalized_surface_mapping();
  normalized_spline_surface_patch = QuadraticSplineSurfacePatch(
    normalized_surface_mapping_coeffs, normalized_domain);
}

PlanarPoint
QuadraticSplineSurfacePatch::denormalize_domain_point(
  PlanarPoint& normalized_domain_point) const
{
  // Get domain triangle vertices
  ConvexPolygon const& domain = get_domain();
  MatrixXr domain_vertices = domain.get_vertices();
  PlanarPoint v0 = domain_vertices.row(0);
  PlanarPoint v1 = domain_vertices.row(1);
  PlanarPoint v2 = domain_vertices.row(2);

  // Generate affine transformation mapping the standard triangle to the domain
  // triangle
  Eigen::Matrix<double, 2, 2> linear_transformation;
  PlanarPoint translation;
  linear_transformation.row(0) = v1 - v0;
  linear_transformation.row(1) = v2 - v0;
  translation = v0;

  // Denormalize the domain point
  return normalized_domain_point * linear_transformation + translation;
}

void
QuadraticSplineSurfacePatch::evaluate(const PlanarPoint& domain_point,
                                      SpatialVector& surface_point) const
{
  evaluate_quadratic_mapping<3>(
    m_surface_mapping_coeffs, domain_point, surface_point);
}

void
QuadraticSplineSurfacePatch::evaluate_normal(
  const PlanarPoint& domain_point,
  SpatialVector& surface_normal) const
{
  evaluate_quadratic_mapping<3>(
    m_normal_mapping_coeffs, domain_point, surface_normal);
}

void
QuadraticSplineSurfacePatch::sample(
  size_t sampling_density,
  std::vector<SpatialVector>& spline_surface_patch_points) const
{
  // Sample the convex domain
  std::vector<PlanarPoint> domain_points;
  m_domain.sample(sampling_density, domain_points);

  // Lift the domain points to the surface
  size_t num_points = domain_points.size();
  spline_surface_patch_points.resize(num_points);
  for (size_t i = 0; i < num_points; ++i) {
    evaluate(domain_points[i], spline_surface_patch_points[i]);
  }
}

void
QuadraticSplineSurfacePatch::triangulate(size_t num_refinements,
                                         Eigen::MatrixXd& V,
                                         Eigen::MatrixXi& F,
                                         Eigen::MatrixXd& N) const
{
  // Triangulate the domain
  Eigen::MatrixXd V_domain;
  m_domain.triangulate(num_refinements, V_domain, F);

  // Lift the domain vertices to the surface and also compute the normals
  V.resize(V_domain.rows(), dimension());
  N.resize(V_domain.rows(), dimension());
  for (Eigen::Index i = 0; i < V_domain.rows(); ++i) {
    SpatialVector surface_point;
    SpatialVector surface_normal;
    evaluate(V_domain.row(i), surface_point);
    evaluate_normal(V_domain.row(i), surface_normal);
    V.row(i) = surface_point;
    N.row(i) = surface_normal;
  }
}

void
QuadraticSplineSurfacePatch::add_patch_to_viewer(std::string patch_name) const
{
  // Generate mesh discretization
  int num_refinements = 2;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd N;
  triangulate(num_refinements, V, F, N);

  // Add patch mesh
  polyscope::init();
  polyscope::registerSurfaceMesh(patch_name, V, F);
}

std::ostream&
operator<<(std::ostream& out,
           const QuadraticSplineSurfacePatch& spline_surface_patch)
{
  out << spline_surface_patch.formatted_patch();
  return out;
}

void
QuadraticSplineSurfacePatch::serialize(std::ostream& out) const
{
  int prec = 17;
  out << "patch" << std::endl;

  // Serialize x coordinate coefficients
  out << "cx";
  for (Eigen::Index i = 0; i < m_surface_mapping_coeffs.rows(); ++i) {
    out << std::setprecision(prec) << " " << m_surface_mapping_coeffs(i, 0);
  }
  out << std::endl;

  // Serialize y coordinate coefficients
  out << "cy";
  for (Eigen::Index i = 0; i < m_surface_mapping_coeffs.rows(); ++i) {
    out << std::setprecision(prec) << " " << m_surface_mapping_coeffs(i, 1);
  }
  out << std::endl;

  // Serialize z coordinate coefficients
  out << "cz";
  for (Eigen::Index i = 0; i < m_surface_mapping_coeffs.rows(); ++i) {
    out << std::setprecision(prec) << " " << m_surface_mapping_coeffs(i, 2);
  }
  out << std::endl;

  // Serialize domain boundary
  MatrixXr const& vertices = m_domain.get_vertices();
  if (get_cone() == 0)
    out << "cp1 ";
  else
    out << "p1 ";
  out << std::setprecision(prec) << vertices(0, 0) << " " << vertices(0, 1)
      << std::endl;
  if (get_cone() == 1)
    out << "cp2 ";
  else
    out << "p2 ";
  out << std::setprecision(prec) << vertices(1, 0) << " " << vertices(1, 1)
      << std::endl;
  if (get_cone() == 2)
    out << "cp3 ";
  else
    out << "p3 ";
  out << std::setprecision(prec) << vertices(2, 0) << " " << vertices(2, 1)
      << std::endl;
}

void
QuadraticSplineSurfacePatch::write_patch(const std::string& filename) const
{
  spdlog::info("Writing spline patch to {}", filename);
  std::ofstream output_file(filename, std::ios::out | std::ios::trunc);
  serialize(output_file);
  output_file.close();
}

// ***************
// Private Methods
// ***************

// Print the formatted surface mapping
std::string
QuadraticSplineSurfacePatch::formatted_patch() const
{
  std::stringstream spline_surface_patch_string;
  spline_surface_patch_string
    << formatted_bivariate_quadratic_mapping<3>(m_surface_mapping_coeffs, 16);

  return spline_surface_patch_string.str();
}
