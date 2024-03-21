#define DO_MOL for(n = 0; n< nMol;n++)

#include <iostream>

int main() {
    int nMol = 10;
    int n;
    DO_MOL {
        std::cout << "Molecule " << n << std::endl;
    }
    return 0;
}