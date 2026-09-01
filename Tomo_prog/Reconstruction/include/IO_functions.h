#ifndef __IO_FUNCTIONS__
#define __IO_FUNCTIONS__
#include <H5Cpp.h>
#include <vector>
#include <complex>
#include <fstream>
void save_complex_volume_hdf5(const std::vector<std::complex<double>>& volume,const std::string& filename,hsize_t Nx, hsize_t Ny, hsize_t Nz,
                              double taille_pixel);


void save_real_hdf5_volume(const std::vector<std::complex<double>>& volume,
                            const std::string& filename,
                            hsize_t Nx, hsize_t Ny, hsize_t Nz, const std::string& dataset_name);

std::vector<std::complex<double>> load_complex_volume_hdf5(const std::string &filename, hsize_t &Nx, hsize_t &Ny, hsize_t &Nz);

void lire_bin(std::string chemin, double resultat[], short int precision, const size_t NbParam);
int get_bin_file_size(std::string chemin);

#endif
