#ifndef PHYSICAL_H
#define PHYSICAL_H
#include "types.h"

void VWrapComponent(real& component, const real& region_component) {
    if (component >= 0.5 * region_component) {
        component -= region_component;
    } else {
        if (component < -0.5 * region_component) {
            component += region_component;
        }
    }
}

// Function to wrap all components of a vector
void VWrapAll(VecR& v, const VecR& region) {
    VWrapComponent(v.x, region.x);
    VWrapComponent(v.y, region.y);
    // Add more lines for additional components if needed
}


#endif // PHYSICAL_H