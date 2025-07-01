#include <iostream>
#include <vector>
#include <complex>
#include <omp.h>
#include "champCplx_functions.h"
#include "omp.h"
double computeEnergy(const std::vector<std::complex<double>>& img) {
    double energy = 0.0;

    #pragma omp parallel for reduction(+:energy)//"reduction" avoid that each threads write in energy, instead create mulple copy and calculate local sums on a limited number of pixel. No conflict in calculation
    for (size_t i = 0; i < img.size(); ++i) {
        energy += std::norm(img[i]);
    }

    return energy;
}
