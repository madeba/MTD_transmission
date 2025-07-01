#include <iostream>

#include "IO_functions.h"
using namespace H5;

#include <H5Cpp.h>

std::vector<std::complex<double>> load_complex_volume_hdf5(const std::string &filename,
                                                        hsize_t &Nx, hsize_t &Ny, hsize_t &Nz)
{
    // 1. Ouverture du fichier
    H5File file(filename, H5F_ACC_RDONLY);

    // 2. Ouverture du groupe
    Group grp = file.openGroup("/indice_complexe");

    // 3. Ouverture des datasets
    DataSet dset_re = grp.openDataSet("Re");
    DataSet dset_im = grp.openDataSet("Im");

    // 4. Récupération des dimensions depuis le dataspace
    DataSpace dataspace = dset_re.getSpace();
    hsize_t dims[3];
    dataspace.getSimpleExtentDims(dims, nullptr);

    // Attribution des dimensions (rappel: stockage Z, Y, X)
    Nz = dims[0];
    Ny = dims[1];
    Nx = dims[2];

    size_t total_size = Nx * Ny * Nz;

    // 5. Lecture des données Re et Im
    std::vector<float> data_re(total_size);
    std::vector<float> data_im(total_size);

    dset_re.read(data_re.data(), PredType::NATIVE_FLOAT);
    dset_im.read(data_im.data(), PredType::NATIVE_FLOAT);
    // dset_re.read(data_re.data(), PredType::NATIVE_DOUBLE);
   // dset_im.read(data_im.data(), PredType::NATIVE_DOUBLE);

    // 6. Reconstruction du volume complexe
    std::vector<std::complex<double>> volume(total_size);
    for (size_t i = 0; i < total_size; ++i) {
        volume[i] = std::complex<double>(data_re[i], data_im[i]);
    }
     return volume;
}



void save_complex_volume_hdf5(const std::vector<std::complex<double>> &volume,
                              const std::string &filename,
                              hsize_t Nx, hsize_t Ny, hsize_t Nz,
                              double pixel_size_m){
    // 1. Création du fichier
    H5File file(filename, H5F_ACC_TRUNC);

    // 2. Création du groupe si nécessaire
    Group grp = file.createGroup("/indice_complexe");

    // 3. Création des DataSpaces
    hsize_t dims[3] = {Nz, Ny, Nx};  // Z, Y, X
    DataSpace dataspace(3, dims);

    // 4. Extraction Re et Im dans buffers
    std::vector<float> data_re(Nx * Ny * Nz);
    std::vector<float> data_im(Nx * Ny * Nz);

    for (size_t i = 0; i < data_re.size(); ++i) {
        data_re[i] = volume[i].real();
        data_im[i] = volume[i].imag();
    }

    // 5. Création des datasets
    DataSet dset_re = grp.createDataSet("Re", PredType::NATIVE_FLOAT, dataspace);
    DataSet dset_im = grp.createDataSet("Im", PredType::NATIVE_FLOAT, dataspace);

    // 6. Écriture
    dset_re.write(data_re.data(), PredType::NATIVE_FLOAT);
    dset_im.write(data_im.data(), PredType::NATIVE_FLOAT);
}



void save_real_hdf5_volume(const std::vector<std::complex<double>>& volume,
                            const std::string& filename,
                            hsize_t Nx, hsize_t Ny, hsize_t Nz,
                            const std::string& dataset_name = "/volume_real")
{
    // Vérification de la taille
    if (volume.size() != Nx * Ny * Nz) {
        throw std::runtime_error("Dimensions incompatibles avec le volume fourni.");
    }

    // Création du tableau temporaire en float
    std::vector<float> data_real(Nx * Ny * Nz);
    for (size_t i = 0; i < volume.size(); ++i) {
        data_real[i] = static_cast<float>(volume[i].real());
    }

    // Dimensions du dataset
    hsize_t dims[3] = { Nz, Ny, Nx };  // Fiji attend ZYX

    try {
        H5::H5File file(filename, H5F_ACC_TRUNC);
        H5::DataSpace dataspace(3, dims);
        H5::DataSet dataset = file.createDataSet(dataset_name,
                                                 H5::PredType::NATIVE_FLOAT,
                                                 dataspace);
        dataset.write(data_real.data(), H5::PredType::NATIVE_FLOAT);
    } catch (H5::Exception& e) {
        std::cerr << "Erreur HDF5 : " << e.getCDetailMsg() << std::endl;
    }
}
