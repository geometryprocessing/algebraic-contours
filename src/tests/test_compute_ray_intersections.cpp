#include <catch2/catch_test_macros.hpp>

#include "common.h"
#include "compute_ray_intersections.h"
#include "generate_shapes.h"

#include "compute_ray_intersections_pencil_method.h"

TEST_CASE("The intersection of a line and a plane can be found",
          "[compute_ray_intersections]") {
  Matrix6x3r surface_mapping_coeffs;
  Matrix2x3r ray_mapping_coeffs;
  std::array<PlanarPoint, MAX_PATCH_RAY_INTERSECTIONS> surface_intersections;
  std::array<double, MAX_PATCH_RAY_INTERSECTIONS> ray_intersections;

//  SECTION("Orthogonal ray") {
//    surface_mapping_coeffs << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        0;
//    ray_mapping_coeffs << 0.5, 0.5, 1, 0, 0, -2;
//
//    compute_normalized_quadratic_surface_ray_intersections(
//        surface_mapping_coeffs, ray_mapping_coeffs, surface_intersections,
//        ray_intersections);
//
//    REQUIRE(surface_intersections.size() == 1);
//    REQUIRE(ray_intersections.size() == 1);
//    REQUIRE(vector_equal(surface_intersections[0], PlanarPoint(0.5, 0.5)));
//    REQUIRE(float_equal(ray_intersections[0], 0.5));
//  }
//
//  SECTION("Parallel ray") {
//    surface_mapping_coeffs << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        0;
//    ray_mapping_coeffs << 0.5, 0.5, 1, 0, 1, 0;
//
//    compute_normalized_quadratic_surface_ray_intersections(
//        surface_mapping_coeffs, ray_mapping_coeffs, surface_intersections,
//        ray_intersections);
//
//    REQUIRE(surface_intersections.size() == 0);
//    REQUIRE(ray_intersections.size() == 0);
//  }
//
//  SECTION("Oblique ray") {
//    surface_mapping_coeffs << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        0;
//    ray_mapping_coeffs << 0.5, 0.5, 1, 0, 0.5, -4;
//
//    compute_normalized_quadratic_surface_ray_intersections(
//        surface_mapping_coeffs, ray_mapping_coeffs, surface_intersections,
//        ray_intersections);
//
//    REQUIRE(surface_intersections.size() == 1);
//    REQUIRE(ray_intersections.size() == 1);
//    REQUIRE(vector_equal(surface_intersections[0], PlanarPoint(0.5, 0.625)));
//    REQUIRE(float_equal(ray_intersections[0], 0.25));
//  }
//
//  SECTION("Folded surface") {
//    surface_mapping_coeffs << 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
//        0;
//    ray_mapping_coeffs << 0.6, 0.5, 1, 0, 0, -2;
//
//    compute_normalized_quadratic_surface_ray_intersections(
//        surface_mapping_coeffs, ray_mapping_coeffs, surface_intersections,
//        ray_intersections);
//
//    REQUIRE(surface_intersections.size() == 1);
//    REQUIRE(ray_intersections.size() == 1);
//  }
//
//  
//
  SECTION("Nonstandard domain") {
    // surface_mapping_coeffs << 0, 0, 0, // 1
    //     1, 0, 0,                       // u
    //     0, 1, 0,                       // v
    //     0, 0, 0,                       // uv
    //     0, 0, 0,                       // uu
    //     0, 0, 0;                       // vv

    surface_mapping_coeffs << 0.1, 0.2, 0, // 1
        1, 0.1, 0.2,                       // u
        0.3, 0.2, 0.6,                     // v
        -0.2, 0.3, 2,                      // uv
        -0.2, 0.1, 1.2,                    // uu
        1, 0.2, 2;                         // vv

    surface_mapping_coeffs << -1.7644781234132338, 2.3739354387396845,
        151.40962125112011,                                                   //
        0.064700996167022434, 0.060355200875465254, 0.11004004230059464,      //
        0.096340254611671616, 0.20934622175013734, -0.037910458133162371,     //
        -0.046780436914869238, 0.03653445519442261, 0.0067183808714837631,    //
        -0.019222335822841597, -0.0027165803799707611, 0.0060945333950915663, //
        -0.01945302670205272, 0.012961148925535924, -0.0050044769478544589;

    std::cout << surface_mapping_coeffs.transpose() << std::endl;
    MatrixXr uv(3, 2);
    // uv << 0.5, 0.5, 1.0, 0.0, 1.0, 1.0;
    uv << 0, 0, 1, 0, 0, 1;

    ConvexPolygon normalized_domain(uv);

    QuadraticSplineSurfacePatch spline_surface_patch(surface_mapping_coeffs,
                                                     normalized_domain);
    // ray_mapping_coeffs << 0.75, 0.6, -0.5, 0, 0, 1;
    // ray_mapping_coeffs << 8.75 / 2.0, 1.5 / 2.0, -0.5, 0, 0, 2;
    ray_mapping_coeffs << -1.69671806538393, 2.5609728762486372,
        51.377018112355586, 0, 0, 200;

    //compute_spline_surface_patch_ray_intersections(
    //    spline_surface_patch, ray_mapping_coeffs, surface_intersections,
    //    ray_intersections);

    std::cout << "here" << std::endl;
    std::cout << "original u v:" << std::endl;
    for (int i = 0; i < surface_intersections.size(); i++) {
      std::cout << surface_intersections[i].transpose() << std::endl;
    }
    std::cout << "original t:" << std::endl;
    for (int i = 0; i < ray_intersections.size(); i++) {
      std::cout << ray_intersections[i] << std::endl;
    }
    int num_intersections;
    long long ray_intersections_call;
    long long ray_bounding_box_call;
    compute_spline_surface_patch_ray_intersections_pencil_method(
        spline_surface_patch, ray_mapping_coeffs, num_intersections, surface_intersections,
        ray_intersections, ray_intersections_call, ray_bounding_box_call);

    std::cout << "here2" << std::endl;

    std::cout << "pencil u v:" << std::endl;
    for (int i = 0; i < num_intersections; i++) {
      std::cout << surface_intersections[i].transpose() << std::endl;
    }
    std::cout << "pencil t:" << std::endl;
    for (int i = 0; i < num_intersections; i++) {
      std::cout << ray_intersections[i] << std::endl;
    }

    // REQUIRE(surface_intersections.size() == 1);
    // REQUIRE(ray_intersections.size() == 1);
    // REQUIRE(vector_equal(surface_intersections[0], PlanarPoint(0.75, 0.6)));
    // REQUIRE(float_equal(ray_intersections[0], 0.5));
  }

  SECTION("Buggy case") {
    //surface_mapping_coeffs <<
    //-0.2693712613262523, -0.5153358428079721, -5.0861375988037265,
    //-0.6495837675648397, -0.4482439315535407,  0.0651557329235988,
     //0.0052419258883517, -0.0314664588921571,  0.0314649989041941,
     //0.0063504149683548,  0.0327342275532451,  0.0179998107390204,
     //0.0521452691525865,  0.3891106316368191,  0.0164684727730394,
     //0.0077627374615072,  0.0097248648585380,  0.0052741903599595;
     
    surface_mapping_coeffs <<
     0.0306990265776488, -0.0083980370575310, -5.0953307552169225,
    -0.0085178318208030, -1.1377644072476367,  0.0007980868147486,
     0.0988863835140743, -1.1312633598946833, -0.0935388401944144,
    -2.6232348137245647,  0.7342847594183667,  0.0073067668313452,
    -1.1976146269412042,  0.4707137286546729, -0.0035880162103296,
    -1.3980539186437533,  0.2347754967876012,  0.2238503060416120;

    std::cout << surface_mapping_coeffs.transpose() << std::endl;
    MatrixXr uv(3, 2);
    uv <<
        -0,  0.5,
         0,   -0,
      0.25, 0.25;

    ConvexPolygon normalized_domain(uv);

    QuadraticSplineSurfacePatch spline_surface_patch(surface_mapping_coeffs,
                                                     normalized_domain);
    ray_mapping_coeffs << 
    -0.4639130718485140, -0.6313829128108502, -25.0169047238889739,
                    0.0,                 0.0,  40.0000000000000000;
    int num_intersections;
    long long ray_intersections_call;
    long long ray_bounding_box_call;
    compute_spline_surface_patch_ray_intersections_pencil_method(
        spline_surface_patch, ray_mapping_coeffs, num_intersections, surface_intersections,
        ray_intersections, ray_intersections_call, ray_bounding_box_call);


    REQUIRE(num_intersections == 0);
    // REQUIRE(vector_equal(surface_intersections[0], PlanarPoint(0.75, 0.6)));
    // REQUIRE(float_equal(ray_intersections[0], 0.5));
  }


}