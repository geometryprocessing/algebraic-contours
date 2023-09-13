#include "conformal_ideal_delaunay/ConformalInterface.hh"
#include "common.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <CLI/CLI.hpp>

/*
std::filesystem::path
join_path(const std::filesystem::path& first_path,
          const std::filesystem::path& second_path)
{
  return first_path / second_path;
}
*/

template <typename T>
void read_vector_from_file(
  const std::string &filename,
  std::vector<T> &vec
) {
  vec.clear();

  // Open file
  std::ifstream input_file(filename);
  if (!input_file) return;

  // Read file
  std::string line;
  while (std::getline(input_file, line))
  {
    std::istringstream iss(line);
    T value;
    iss >> value;
    vec.push_back(value);
  }

  // Close file
  input_file.close();
}

template <typename VectorScalar, typename MatrixScalar>
void
convert_std_to_eigen_vector(const std::vector<VectorScalar>& vector_std,
                            Eigen::Matrix<MatrixScalar, Eigen::Dynamic, 1>& vector_eigen)
{
  size_t vector_size = vector_std.size();
  vector_eigen.resize(vector_size);
  for (size_t i = 0; i < vector_size; ++i) {
    vector_eigen[i] = MatrixScalar(vector_std[i]);
  }
}

template <typename VectorScalar, typename MatrixScalar>
void
convert_std_to_eigen_matrix(
  const std::vector<std::vector<VectorScalar>>& matrix_vec,
  Eigen::Matrix<MatrixScalar, Eigen::Dynamic, Eigen::Dynamic>& matrix
) {
  matrix.setZero(0, 0);
  if (matrix_vec.empty()) return;

  // Get dimensions of matrix
  int rows = matrix_vec.size();
  int cols = matrix_vec[0].size();
	matrix.resize(rows, cols);
  
  // Copy matrix by row
	for (int i = 0; i < rows; ++i)
	{
    // Check size validity
    if (static_cast<int>(matrix_vec[i].size()) != cols)
    {
      spdlog::error("Cannot copy vector of vectors of inconsistent sizes to a matrix");
      matrix.setZero(0, 0);
      return;
    }

    // Copy row
		for (int j = 0; j < cols; ++j)
		{
			matrix(i, j) = MatrixScalar(matrix_vec[i][j]);
		}
	}
}

void
parameterize_component(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	const std::vector<double>& Th_hat,
	Eigen::MatrixXd& V_o,
	Eigen::MatrixXi& F_o,
	Eigen::MatrixXd& uv_o,
	Eigen::MatrixXi& FT_o
) {
	// Get conformal mesh
	std::vector<int> pt_fids(0);
	std::vector<Eigen::Matrix<double, 3, 1>> pt_bcs(0);
	auto conformal_res = conformal_metric<double>(
		V,
		F,
		Th_hat,
		pt_fids,
		pt_bcs
	);
	OverlayMesh<double> m_o = std::get<0>(conformal_res);
	std::vector<double> u = std::get<1>(conformal_res);
	std::vector<int> vtx_reindex = std::get<4>(conformal_res);
  std::vector<std::vector<double>> V_overlay = std::get<5>(conformal_res);
  std::vector<std::pair<int,int>> endpoints = std::get<6>(conformal_res);

	// Build output mesh
	auto parametrize_res = overlay_mesh_to_VL<double>(
		V,
		F,
		Th_hat,
		m_o,
		u,
		V_overlay,
		vtx_reindex,
		endpoints,
		-1
	);
	std::vector<std::vector<double>> V_o_vec = std::get<0>(parametrize_res);
	std::vector<std::vector<int>> F_o_vec = std::get<1>(parametrize_res);
	std::vector<double> u_o_vec = std::get<2>(parametrize_res);
	std::vector<double> v_o_vec = std::get<3>(parametrize_res);
	std::vector<std::vector<int>> FT_o_vec = std::get<4>(parametrize_res);

	// Convert vector formats to matrices
	Eigen::VectorXd u_o, v_o;
	convert_std_to_eigen_matrix(V_o_vec, V_o);
	convert_std_to_eigen_matrix(F_o_vec, F_o);
	convert_std_to_eigen_matrix(FT_o_vec, FT_o);
	convert_std_to_eigen_vector(u_o_vec, u_o);
	convert_std_to_eigen_vector(v_o_vec, v_o);
	uv_o.resize(u_o.size(), 2);
	uv_o.col(0) = u_o;
	uv_o.col(1) = v_o;
}

int main(int argc, char *argv[])
{
  CLI::App app{"Parameterize a mesh."};
  std::string input_dir = "./";
  std::string mesh_name = "";
  std::string output_dir = "./";
  app.add_option("-i,--input", input_dir, "Mesh input directory")
    ->check(CLI::ExistingDirectory)
    ->required();
  app.add_option("-f,--fname", mesh_name, "Mesh filename")
    ->required();
  app.add_option("-o,--output", output_dir, "Output directory")
    ->check(CLI::ExistingDirectory);
  CLI11_PARSE(app, argc, argv);
  spdlog::set_level(spdlog::level::info);
  std::filesystem::create_directories(output_dir);
  std::string input_filename = join_path(input_dir, mesh_name + ".obj");
  std::string Th_hat_filename = join_path(input_dir, mesh_name + "_Th_hat");

  // Get input mesh
  Eigen::MatrixXd V, uv, N;
  Eigen::MatrixXi F, FT, FN;
	spdlog::info("Optimizing mesh at {}", input_filename);
  igl::readOBJ(input_filename, V, uv, N, F, FT, FN);
	
	// Get input angles
	std::vector<double> Th_hat;
	spdlog::info("Using cone angles at {}", Th_hat_filename);
	read_vector_from_file(Th_hat_filename, Th_hat);

  // Split mesh into components
  std::vector<Eigen::MatrixXd> V_components;
  std::vector<Eigen::MatrixXi> F_components;
  std::vector<Eigen::VectorXi> I_components, J_components;
  separate_mesh_components(V, F, V_components, F_components, I_components, J_components);
  int num_components = V_components.size();

  // Split cone angles
  std::vector<std::vector<double>> Th_hat_components(num_components);
  for (int i = 0; i < num_components; ++i)
  {
    int num_component_vertices = V_components[i].rows();
    Th_hat_components[i].resize(num_component_vertices);
    for (int j = 0; j < num_component_vertices; ++j)
    {
      Th_hat_components[i][j] = Th_hat[J_components[i][j]];
    }
  }

	// Get cone vertices and angles for each component
  std::vector<Eigen::MatrixXd> V_o_components(num_components);
  std::vector<Eigen::MatrixXi> F_o_components(num_components);
  std::vector<Eigen::MatrixXd> uv_o_components(num_components);
  std::vector<Eigen::MatrixXi> FT_o_components(num_components);
  for (int i = 0; i < num_components; ++i)
  {
    parameterize_component(
      V_components[i],
      F_components[i],
      Th_hat_components[i],
      V_o_components[i],
      F_o_components[i],
      uv_o_components[i],
      FT_o_components[i]
    );
  }

  // Recombine faces and vertices
  // WARNING: This will generally not be the same as the original mesh due
  // to reindexing and reordering
	Eigen::MatrixXd V_o;
	Eigen::MatrixXi F_o;
	Eigen::MatrixXd uv_o;
	Eigen::MatrixXi FT_o;
  combine_mesh_components(V_o_components, F_o_components, V_o, F_o);
  combine_mesh_components(uv_o_components, FT_o_components, uv_o, FT_o);

	// Write the output
	std::string output_filename = join_path(output_dir, mesh_name + ".obj");
	Eigen::MatrixXd N_o;
  Eigen::MatrixXi FN_o;
  igl::writeOBJ(output_filename, V_o, F_o, N_o, FN_o, uv_o, FT_o);
}

