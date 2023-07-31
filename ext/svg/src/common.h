#pragma once

#include <Eigen/Core>
#include <iostream>
#include <vector>

// Real number representations
namespace svg
{
typedef Eigen::VectorXd VectorXr;
typedef Eigen::Vector3d Vector3r;
typedef Eigen::RowVectorXd OneFormXr;
typedef Eigen::MatrixXd MatrixXr;


typedef double real_t;
using Color = Eigen::Vector4f;
}
