#include <iostream>

#include "IO_functions.h"
using namespace H5;
using namespace std;
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


///#########lecture  d'un fichier binaire 3D, connaissant sa taille et son type de données
///read a  binary file. (path, 3D table, data format (double=64), Nbpixels to be read). for data format, can use enum PRECISION, cf "projet.h"
int get_bin_file_size(string chemin)
{
    size_t lTaille, nb_elmnt_lu;//size and number of elements
    //unsigned short int dimData=precision/8;//taille en octet d'un element.
    FILE* pFichier = NULL;
    pFichier = fopen(chemin.c_str(), "r");  //ouverture de ce fichier en écriture binaire

    if(pFichier==NULL){
        fputs("Impossible d'ouvrir le fichier\n",stderr);
        cout<<chemin<<endl;
        exit (1);// obtenir la longueur du fichier, comparer avec donnée entrée.
    }
    else{
        fseek(pFichier,0,SEEK_END);//trouver la fin de fichier
        lTaille = ftell (pFichier);//retourne la position courante (en octet) du curseur de fichier : ici, position de la fin du fichier
        rewind(pFichier);
        fclose(pFichier);
    }
return lTaille;
}

///#########lecture  d'un fichire binaire 3D, connaissant sa taille et son type de données
///read a  binary file. (path, 3D table, data format (double=64), Nbpixels to be read). for data format, can use enum PRECISION, cf "projet.h"
void lire_bin(string chemin, double resultat[], short int precision, const size_t NbPix)
{
    size_t lTaille, nb_elmnt_lu;//size and number of elements
    unsigned short int dimData=precision/8;//taille en octet d'un element.
    FILE* pFichier = NULL;
    pFichier = fopen(chemin.c_str(), "r");  //ouverture de ce fichier en écriture binaire

    if(pFichier==NULL){
        fputs("Impossible d'ouvrir le fichier\n",stderr);
        cout<<chemin<<endl;
        exit (1);// obtenir la longueur du fichier, comparer avec donnée entrée.
    }
    else{
        fseek(pFichier,0,SEEK_END);//trouver la fin de fichier
        lTaille = ftell (pFichier);//retourne la position courante (en octet) du curseur de fichier : ici, position de la fin du fichier->taille du fichier
        cout<<"Fichier "<<chemin<<endl;
         printf("taille trouvée en octet par ftell %li, taille estimée : %i\n",lTaille, NbPix*dimData);//
        rewind(pFichier);

        if(NbPix*dimData!=lTaille)
            cout<<"Taille du fichier "<<chemin <<" incompatible avec les dimensions\n"<<endl;

        nb_elmnt_lu = fread (resultat,1,lTaille,pFichier);//lecture
        //nb_elmnt_lu = fread (&resultat_vector[0], 1,lTaille,pFichier);

        if(nb_elmnt_lu!=lTaille){
            cout<<"Problème lors de la lecture du fichier "<<chemin<<endl;
            cout<<"Nombre d'éléments lus="<<nb_elmnt_lu<<endl;
        }
        fclose(pFichier);
    }

}
///surcharge RAII
vector<double> lire_bin(const string& chemin, short int precision, size_t NbPix)
{
    ifstream fichier(chemin, ios::binary | ios::ate);

    if (!fichier) {
        throw runtime_error("Impossible d'ouvrir le fichier: " + chemin);
    }

    size_t tailleFichier = fichier.tellg();
    size_t tailleAttendue = NbPix * (precision / 8);

    if (tailleFichier != tailleAttendue) {
        throw runtime_error("Taille incompatible: " + to_string(tailleFichier)
                          + " vs " + to_string(tailleAttendue));
    }
    fichier.seekg(0);
    vector<unsigned char> buffer(tailleFichier);
    fichier.read(reinterpret_cast<char*>(buffer.data()), tailleFichier);

    vector<double> resultat(NbPix);

    // Conversion selon la précision
    for (size_t i = 0; i < NbPix; i++) {
        if (precision == 8) {
            resultat[i] = buffer[i];
        }
        else if (precision == 16) {
            resultat[i] = *reinterpret_cast<unsigned short*>(&buffer[i * 2]);
        }
        else if (precision == 32) {
            resultat[i] = *reinterpret_cast<float*>(&buffer[i * 4]);
        }
    }
    return resultat;
}
