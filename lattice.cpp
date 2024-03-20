#include <iostream>
#include <vector>

using namespace std;

struct Point {
    int x, y;
    Point(int _x, int _y) : x(_x), y(_y) {}
};

vector<Point> generateSquareLattice(int width, int height, int spacing) {
    vector<Point> latticePoints;
    for (int x = 0; x < width; x += spacing) {
        for (int y = 0; y < height; y += spacing) {
            latticePoints.emplace_back(x, y);
        }
    }
    return latticePoints;
}

void visualizeLattice(const vector<Point>& latticePoints, int width, int height) {
    vector<vector<bool>> grid(height, vector<bool>(width, false));

    // Mark lattice points on the grid
    for (const auto& point : latticePoints) {
        grid[point.y][point.x] = true;
    }

    // Print the grid
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (grid[y][x]) {
                cout << "X "; // Print lattice point
            } else {
                cout << ". "; // Print empty space
            }
        }
        cout << endl;
    }
}

int main() {
    int width = 10;  // Width of the lattice
    int height = 10; // Height of the lattice
    int spacing = 2; // Spacing between lattice points

    vector<Point> latticePoints = generateSquareLattice(width, height, spacing);

    visualizeLattice(latticePoints, width, height);

    return 0;
}
