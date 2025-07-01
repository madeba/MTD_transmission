#ifndef __PHASE_FUNCTIONS__
#define __PHASE_FUNCTIONS__

#include <vector>
#include <complex>

void calcPhase_mpi_pi_atan2(std::vector<std::complex<double>> const &cplxField, std::vector<double> &phaseMod2pi);

#endif
