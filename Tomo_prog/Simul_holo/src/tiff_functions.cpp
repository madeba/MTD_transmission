#include "tiff_functions.h"
string formatFijiMsgMicron(double pixelSize)
{ char spacing_line[64];
    snprintf(spacing_line, sizeof(spacing_line), "spacing=%.6f\n", pixelSize);

    string msg= "ImageJ=\n";
   msg += spacing_line;
    msg += "unit=micron\n";
    return msg;
}
///get Tiff Tags : width,height,depth=x,y,z
int GetTiff3D_Dimensions(const std::string& path, Var3D &dim)
 {
    uint32 width,height,depth;

    TIFF* tif = TIFFOpen(path.c_str(), "r");
    if (!tif) return -1;

    depth = 0;

    // Lire la première page pour width/height
    if (TIFFSetDirectory(tif, 0)) {
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    } else {
        TIFFClose(tif);
        return -1;
    }

    // Compter le nombre de pages = profondeur Z
    do {
        ++depth;
    } while (TIFFReadDirectory(tif));

    TIFFClose(tif);
    dim.x=width;
    dim.y=height;
    dim.z=depth;
    return 0;
}

void Import3D_Tiff(std::vector<double>& imgTiff, const std::string& chemin)
{
    TIFF* tif = TIFFOpen(chemin.c_str(), "r");
    if (!tif) {
        std::cerr << "Erreur : impossible d’ouvrir le fichier TIFF : " << chemin << std::endl;
        return;
    }

    uint32 width = 0, height = 0;
    uint16 bitsPerSample = 0, sampleFormat = 0, samplesPerPixel = 0;

    // Lire la première page pour obtenir les dimensions X et Y
    if (!TIFFSetDirectory(tif, 0)) {
        std::cerr << "Erreur : la première page TIFF est inaccessible." << std::endl;
        TIFFClose(tif);
        return;
    }

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &sampleFormat);

    if (bitsPerSample != 32 || sampleFormat != SAMPLEFORMAT_IEEEFP || samplesPerPixel != 1) {
        std::cerr << "Erreur : seules les images TIFF float32 monochromes sont prises en charge." << std::endl;
        TIFFClose(tif);
        return;
    }

    // Compter le nombre total de pages = profondeur Z
    uint32 depth = 0;
    do {
        ++depth;
    } while (TIFFReadDirectory(tif));

    // Réinitialiser la position à la première page
    TIFFSetDirectory(tif, 0);

    std::cout << "Dimensions détectées : " << width << " x " << height << " x " << depth << std::endl;

    // Allouer l’image 3D
    imgTiff.resize(static_cast<size_t>(width) * height * depth);

    std::vector<float> buffer(width * height);  // buffer temporaire en float
    for (uint32 z = 0; z < depth; ++z) {
        if (!TIFFSetDirectory(tif, z)) {
            std::cerr << "Erreur : impossible de lire la page " << z << std::endl;
            TIFFClose(tif);
            return;
        }

        for (uint32 y = 0; y < height; ++y) {
            float* rowPtr = &buffer[y * width];
            if (TIFFReadScanline(tif, rowPtr, y, 0) < 0) {
                std::cerr << "Erreur lecture ligne " << y << " de la page " << z << std::endl;
                TIFFClose(tif);
                return;
            }
        }

        // Copier les valeurs dans imgTiff avec conversion float -> double
        for (uint32 y = 0; y < height; ++y) {
            for (uint32 x = 0; x < width; ++x) {
                size_t index = static_cast<size_t>(z) * width * height + y * width + x;
                imgTiff[index] = static_cast<double>(buffer[y * width + x]);
            }
        }
    }

    TIFFClose(tif);
}
