#include <iostream>

// Define structures for vectors
typedef struct {
    double x, y, z; // Assuming three-dimensional vectors
} VecR;

// Function to perform element-wise division of vectors
void VDiv(VecR &result, const VecR &numerator, const VecR &denominator) {
    result.x = numerator.x / denominator.x;
    result.y = numerator.y / denominator.y;
    result.z = numerator.z / denominator.z;
}

int main() {
    // Example usage of VDiv
    VecR region = {10.0, 10.0, 10.0}; // Example region size vector
    VecR initUcell = {2.0, 2.0, 2.0}; // Example unit cell dimensions vector
    VecR gap;

    // Perform element-wise division using VDiv
    VDiv(gap, region, initUcell);

    // Output the result
    std::cout << "Gap between unit cells: (" << gap.x << ", " << gap.y << ", " << gap.z << ")" << std::endl;

    return 0;
}
