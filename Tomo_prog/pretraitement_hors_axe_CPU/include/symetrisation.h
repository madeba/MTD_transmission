#ifndef SYMETRISATION_H
#define SYMETRISATION_H



#include <vector>
#include <complex>
#include <iostream>
//symetrisation pour scalaire
void Symetrise_mirror2(const std::vector<double>& img,  std::vector<double>& imgSym);
//symétrisation pour gradient ( avant intégration)
void Symetrise_mirror_gradient(const std::vector<double>& gradx, const std::vector<double>& grady, std::vector<double>& gradxSym, std::vector<double>& gradySym);

void Symetrise_mirror2(
    const std::vector<std::complex<double>>& img,
    std::vector<std::complex<double>>& imgSym);//récupérer le 1e cadrant//get back the first quadrant
std::vector<double> cut_quad4(const std::vector<double>& img4);
std::vector<double> cut_quad4(std::vector<std::complex<double>> const &monImg4Quad);
std::vector<std::complex<double>> cut_quad4_cplx(std::vector<std::complex<double>> const &monImg4Quad);
#endif
