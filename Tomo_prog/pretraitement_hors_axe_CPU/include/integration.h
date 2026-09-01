
#ifndef __INTEGRATION__
#define __INTEGRATION_
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <complex>
#include "vecteur.h"
#include "FFTW_init.h"
void integ_grad4(std::vector<double> const &gradx, std::vector<double> const& grady, std::vector<std::complex<double>> &sortie,std::vector<vecteur> &kvect_shift, FFTW_init &param_c2c);

struct GradientResult {
    std::vector<double> phase;
    std::vector<double> gradx;
    std::vector<double> grady;
};

Eigen::VectorXd integrateGradientEigen(
    const Eigen::MatrixXd& gradX,
    const Eigen::MatrixXd& gradY,
    int dimX,
    int dimY);

GradientResult generateEllipseWithGradient(
    Var2D centre,
    Var2D rayon,
    Var2D dim,
    double phase_max);
#endif
