#ifndef REGULARIZED_DIVISION_H
#define REGULARIZED_DIVISION_H

#include <vector>
#include <complex>
#include <iostream>

std::vector <double> gradU_U_dampedDiv(std::vector<std::complex<double>> &UBorn,std::vector<std::complex<double>> &gradUBorn, double alpha);
std::vector<double> global_damp(std::vector<std::complex<double>> &gradUBorn, std::vector<std::complex<double>> UBorn, double alpha_global);

//std::vector<double> local_damp(std::vector<std::complex<double>> &gradUBorn, std::vector<std::complex<double>> UBorn, double alpha_global, double beta_local);
std::vector<double> local_damp(const std::vector<std::complex<double>> &gradUBorn,
                          const std::vector<std::complex<double>> &UBorn,
                          double alpha_global,
                          double beta_local=0.1,
                          int margin=2,
                          bool masque_jointure = true);

#endif
