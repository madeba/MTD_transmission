#ifndef TIFF_WRITER_H
#define TIFF_WRITER_H
#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <tiffio.h>
#include "projet.h"
#include "omp.h"
///SAV3Dtiff avec template
std::string formatFijiMsgMicron(double pixelSize);


template<typename T>
void write3D_Tiff(const std::vector<T>& var_sav,
                const Var3D& dim,
                const std::string& partie,
                const std::string& chemin,
                double taille_pixel_m,
                const std::string& description = "") {
    taille_pixel_m=1e6*taille_pixel_m;//conversion in micron for imageJ
    std::string msgCalibZ=formatFijiMsgMicron(taille_pixel_m);
    // Vérification de cohérence des données
    if (var_sav.size() != static_cast<size_t>(dim.x * dim.y * dim.z)) {
        throw std::invalid_argument("Taille du vecteur incompatible avec les dimensions");
    }

    if (taille_pixel_m <= 0.0) {
        throw std::invalid_argument("La taille de pixel doit être positive");
    }

    // Vérification de la partie pour les types complexes
    if constexpr (std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<float>>) {
        if (partie != "Re" && partie != "re" && partie != "Im" && partie != "im") {
            throw std::invalid_argument("Partie non reconnue pour type complexe : doit être 'Re' ou 'Im'");
        }
    } else {
        // Pour les types réels, on ignore le paramètre partie (ou on peut vérifier qu'il est vide)
        if (!partie.empty() && partie != "Re" && partie != "re") {
            throw std::invalid_argument("Pour les types réels, 'partie' doit être vide ou 'Re'");
        }
    }

    // Ouverture du fichier avec gestion d'erreur
    TIFF* out = TIFFOpen(chemin.c_str(), "w");
    if (!out) {
        throw std::runtime_error("Erreur à l'ouverture du fichier TIFF : " + chemin);
    }

    // Configuration des paramètres TIFF
    const uint32_t width = static_cast<uint32_t>(dim.x);
    const uint32_t height = static_cast<uint32_t>(dim.y);
    const uint32_t depth = static_cast<uint32_t>(dim.z);
    const uint16_t spp = 1;
    const uint16_t bpp = 32;

    // Résolution en pixels/cm (libtiff ne supporte par le metre)
    const float resolution = 1 / static_cast<float>(taille_pixel_m);//100 cm par metre->10 000/taille pixel pixel par cm
   // cout<<"resolution="<<resolution<<endl;
    // Déterminer si on veut la partie réelle (pour les complexes)
    const bool is_real = (partie.empty() || partie == "Re" || partie == "re");

    // Buffer 2D réutilisé pour chaque slice
    std::vector<float> buffer2D(width * height);

    try {
        for (uint32_t z = 0; z < depth; ++z) {
            const size_t offset_z = static_cast<size_t>(z) * width * height;

            // Extraction des données selon le type
            for (uint32_t y = 0; y < height; ++y) {
                const size_t offset_y = static_cast<size_t>(y) * width;
                for (uint32_t x = 0; x < width; ++x) {
                    const auto& val = var_sav[offset_z + offset_y + x];
                    const size_t idx = offset_y + x;

                    if constexpr (std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<float>>) {
                        // Type complexe
                        buffer2D[idx] = is_real ? static_cast<float>(val.real()) :
                                                 static_cast<float>(val.imag());
                    } else {
                        // Type réel (double, float, etc.)
                        buffer2D[idx] = static_cast<float>(val);
                    }
                }
            }

            // Configuration des tags TIFF pour cette page
            TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
            TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
            TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
            TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
            TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
            TIFFSetField(out, TIFFTAG_XRESOLUTION, resolution);
            TIFFSetField(out, TIFFTAG_YRESOLUTION, resolution);
            TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, RESUNIT_CENTIMETER);
            TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            TIFFSetField(out, TIFFTAG_PAGENUMBER, z, depth);
            TIFFSetField(out, TIFFTAG_DOCUMENTNAME, description);

            // Métadonnées optionnelles
            if (!msgCalibZ.empty()) {
                std::string page_desc = msgCalibZ + " (slice " + std::to_string(z + 1) +
                                       "/" + std::to_string(depth) + ")";
                TIFFSetField(out, TIFFTAG_IMAGEDESCRIPTION, page_desc.c_str());
            }

            // Métadonnées de pixel
            TIFFSetField(out, TIFFTAG_SOFTWARE, "Custom 3D TIFF Writer");

            // Écriture ligne par ligne
            for (uint32_t y = 0; y < height; ++y) {
                if (TIFFWriteScanline(out, &buffer2D[y * width], y, 0) < 0) {
                    throw std::runtime_error("Erreur lors de l'écriture de la ligne " +
                                           std::to_string(y) + " de la slice " + std::to_string(z));
                }
            }

            // Finaliser cette page
            if (!TIFFWriteDirectory(out)) {
                throw std::runtime_error("Erreur lors de la finalisation de la slice " + std::to_string(z));
            }
        }
    }
    catch (...) {
        TIFFClose(out);
        throw; // Re-lancer l'exception après fermeture du fichier
    }

    TIFFClose(out);
}

// Surcharges pratiques pour simplifier l'usage avec des types réels
template<typename T>
typename std::enable_if_t<!std::is_same_v<T, std::complex<double>> &&
                         !std::is_same_v<T, std::complex<float>>, void>
write3D_Tiff(const std::vector<T>& var_sav,
           const Var3D& dim,
           const std::string& chemin,
           double taille_pixel_m,
           const std::string& description = "") {
    write3D_Tiff(var_sav, dim, "", chemin, taille_pixel_m, description);
}

void Import3D_Tiff(std::vector<double>& imgTiff, const std::string& chemin);
int GetTiff3D_Dimensions(const std::string& path, Var3D &dim);
#endif
