//
// Created by konrad_guest on 27/07/2024.
//
#include "gtest/gtest.h"
#include "Headers/simulation/boundary_conditions.h"

TEST(BoundaryConditionsTest,WrapComponent2D){
    // Test case 1: component >= 0.5 * regionComponent
    double component1 = 6.0;
    double regionComponent1 = 10.0;
    wrap_component_2d(component1, regionComponent1);
    EXPECT_DOUBLE_EQ(component1, -4.0);

    // Test case 2: component < -0.5 * regionComponent
    double component2 = -6.0;
    double regionComponent2 = 10.0;
    wrap_component_2d(component2, regionComponent2);
    EXPECT_DOUBLE_EQ(component2, 4.0);

    // Test case 3: -0.5 * regionComponent <= component < 0.5 * regionComponent
    double component3 = 3.0;
    double regionComponent3 = 10.0;
    wrap_component_2d(component3, regionComponent3);
    EXPECT_DOUBLE_EQ(component3, 3.0);
}
