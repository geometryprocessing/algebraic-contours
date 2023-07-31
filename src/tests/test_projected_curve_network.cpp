#include <catch2/catch_test_macros.hpp>
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "common.h"
#include "projected_curve_network.h"
#include "generate_shapes.h"


TEST_CASE ( "A curve network can be built from topology information", "[projected_curve_network]")
{
  SECTION ( "Closed loop")
  {
    std::vector<AbstractCurveNetwork::SegmentIndex> out_array = { 0, 1, 2};
    std::vector<AbstractCurveNetwork::NodeIndex> to_array = { 1, 2, 0};
    std::vector<AbstractCurveNetwork::NodeIndex> intersections = { -1, -1, -1 };
    AbstractCurveNetwork curve_network(to_array, out_array, intersections);

    REQUIRE( curve_network.next(0) == 1 );
    REQUIRE( curve_network.next(1) == 2 );
    REQUIRE( curve_network.next(2) == 0 );
    REQUIRE( curve_network.prev(0) == 2 );
    REQUIRE( curve_network.prev(1) == 0 );
    REQUIRE( curve_network.prev(2) == 1 );
    REQUIRE( curve_network.to(0) == 1 );
    REQUIRE( curve_network.to(1) == 2 );
    REQUIRE( curve_network.to(2) == 0 );
    REQUIRE( curve_network.from(0) == 0 );
    REQUIRE( curve_network.from(1) == 1 );
    REQUIRE( curve_network.from(2) == 2 );
    REQUIRE( curve_network.out(0) == 0 );
    REQUIRE( curve_network.out(1) == 1 );
    REQUIRE( curve_network.out(2) == 2 );
    REQUIRE( curve_network.in(0) == 2 );
    REQUIRE( curve_network.in(1) == 0 );
    REQUIRE( curve_network.in(2) == 1 );
  }
}