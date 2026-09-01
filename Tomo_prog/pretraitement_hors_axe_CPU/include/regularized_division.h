#ifndef REGULARIZED_DIVISION_H
#define REGULARIZED_DIVISION_H

#include <vector>
#include <complex>
#include <iostream>

std::vector <double> gradU_U_dampedDiv(std::vector<std::complex<double>> &UBorn,std::vector<std::complex<double>> &gradUBorn, double alpha);
std::vector<double> global_damp(std::vector<std::complex<double>> &gradUBorn, std::vector<std::complex<double>> UBorn, double alpha_global);
std::vector<double> local_damp(std::vector<std::complex<double>> &gradUBorn, std::vector<std::complex<double>> UBorn, double alpha_global, double beta_local);
#endif
