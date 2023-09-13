#include "generate_optimal_cones.h"
#include "common.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include "igl/readOBJ.h"
#include "igl/writeOBJ.h"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include <CLI/CLI.hpp>

std::filesystem::path
join_path(const std::filesystem::path& first_path,
          const std::filesystem::path& second_path)
{
  return first_path / second_path;
}

/*
template <typename T>
void write_vector(
	const std::vector<T> &vec,
	const std::string &filename
) {
	std::ofstream output_file(filename, std::ios::out | std::ios::trunc);
  for (size_t i = 0; i < vec.size(); ++i)
  {
    output_file << vec[i] << std::endl;
  }
  output_file.close();
}
*/

void write_vector(
	const std::vector<double> &vec,
	const std::string &filename,
	int precision
) {
	std::ofstream output_file(filename, std::ios::out | std::ios::trunc);
  for (size_t i = 0; i < vec.size(); ++i)
  {
    output_file << std::setprecision(precision) << vec[i] << std::endl;
  }
  output_file.close();
}


int main(int argc, char *argv[])
{
  CLI::App app{"Generate cones for a mesh."};
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

  // Get input mesh
  Eigen::MatrixXd V, uv, N;
  Eigen::MatrixXi F, FT, FN;
  std::string input_filename = join_path(input_dir, mesh_name + ".obj");
	spdlog::info("Generating cones for mesh at {}", input_filename);
  igl::readOBJ(input_filename, V, uv, N, F, FT, FN);

  // Split mesh into components
  std::vector<Eigen::MatrixXd> V_components;
  std::vector<Eigen::MatrixXi> F_components;
  std::vector<Eigen::VectorXi> I_components, J_components;
  separate_mesh_components(V, F, V_components, F_components, I_components, J_components);
  int num_components = V_components.size();

	// Get cone vertices and angles for each component
  int cone_index_offset = 0;
  std::vector<int> cone_vertices(0);
  std::vector<double> cone_angles(0);
  for (int i = 0; i < num_components; ++i)
  {
    std::vector<int> cone_vertices_component;
    std::vector<double> cone_angles_component;
    generate_optimal_cones(
      V_components[i],
      F_components[i],
      cone_vertices_component,
      cone_angles_component);
    for (auto ci : cone_vertices_component)
    {
      cone_vertices.push_back(ci + cone_index_offset);
    }
    append(cone_angles, cone_angles_component);
    cone_index_offset += V_components[i].rows();
  }

  // Recombine faces and vertices
  // WARNING: This will generally not be the same as the original mesh due
  // to reindexing and reordering
  combine_mesh_components(V_components, F_components, V, F);

	// Write the output
	std::string output_filename;
	output_filename = join_path(output_dir, mesh_name + ".obj");
  igl::writeOBJ(output_filename, V, F, N, FN, uv, FT);
	output_filename = join_path(output_dir, mesh_name + "_cones");
	write_vector(cone_vertices, output_filename);
	output_filename = join_path(output_dir, mesh_name + "_Th_hat");
	write_vector(cone_angles, output_filename, 17);
}

