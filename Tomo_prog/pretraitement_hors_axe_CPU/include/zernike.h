#ifndef __ZERNIKE__
#define __ZERNIKE__
#pragma once

#include <vector>
#include <complex>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>
class Zernike
{
public:

    Zernike(int width,
            int height,
            double radius,
            double cx,
            double cy);

    // Fit phase -> coefficients
    Eigen::VectorXd fitPhase(const Eigen::MatrixXd& phase,
                             int maxMode);

    // Rebuild phase from coeffs
    Eigen::MatrixXd reconstruct(const Eigen::VectorXd& coeffs);

    // Correct complex field
    Eigen::MatrixXcd correctWavefront(
        const Eigen::MatrixXcd& field,
        const Eigen::VectorXd& coeffs);

private:

    int W,H;
    double R;
    double CX,CY;

    double radial(int n, int m, double rho);

    double evaluateZernike(int n, int m,
                   double rho,
                   double theta);

    std::vector<std::pair<int,int>>
    generateModes(int maxMode);
};
#endif

