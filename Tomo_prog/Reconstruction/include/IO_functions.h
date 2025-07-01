
#ifndef __IO_FUNCTIONS__
#define __IO_FUNCTIONS__
#include <H5Cpp.h>
#include <vector>
#include <complex>

void save_complex_volume_hdf5(const std::vector<std::complex<double>>& volume,const std::string& filename,hsize_t Nx, hsize_t Ny, hsize_t Nz,
                              double taille_pixel);


void save_real_hdf5_volume(const std::vector<std::complex<double>>& volume,
                            const std::string& filename,
                            hsize_t Nx, hsize_t Ny, hsize_t Nz, const std::string& dataset_name);

std::vector<std::complex<double>> load_complex_volume_hdf5(const std::string &filename,
                                                        hsize_t &Nx, hsize_t &Ny, hsize_t &Nz);


#endif
