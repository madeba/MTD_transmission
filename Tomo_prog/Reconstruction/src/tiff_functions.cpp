#include "tiff_functions.h"
using namespace std;
///fonction formattant les metadat imageJ personnalisé (dont la taille en z).
///spécifiquement la taille de pixels en z, non prise en charge en standard
std::string formatFijiMsgMicron(double pixelSize)
{ char spacing_line[64];
    snprintf(spacing_line, sizeof(spacing_line), "spacing=%.6f\n", pixelSize);

    string msg= "ImageJ=\n";
   msg += spacing_line;
    msg += "unit=micron\n";

    return msg;

    // Création des métadonnées ImageJ
/*std::ostringstream ij_metadata;
ij_metadata.precision(10);
ij_metadata << std::fixed;

ij_metadata << "ImageJ=1.53c\n"
            << "images=" << dim << "\n"
            << "channels=1\n"
            << "slices=" << dim << "\n"
            << "frames=1\n"
            << "hyperstack=true\n"
            << "mode=grayscale\n"
            << "unit=um\n"
            << "spacing=" << (taille_pixel_z * 1e6) << "\n"
            << "loop=false\n";
        return ij_metadata; */
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

void Import3D_Tiff(vector<double> &imgTiff, string chemin, double taille_pixel)
{
    const size_t dim=round(std::pow(imgTiff.size(), 1.0/3.0));
    uint32 image_width, image_height, dimz;
//    float xres, yres;
    uint16 spp; //autres arguments : photo, res_unit, zpage, bpp
    TIFF *Tiff_id;
    size_t x, y; //z;
    float *buffer2D=new float[dim * dim];
    Tiff_id = TIFFOpen(chemin.c_str(), "r");
    if (!Tiff_id)
        fprintf (stderr, "Can't open  for writing\n");

    image_width = dim;
    image_height = dim;
    dimz=dim;
    spp = 1; /* Samples per pixel */
    size_t num_page=0;

    for(num_page = 0; num_page < dim; num_page++) //z=page
    {
//z=page
        TIFFGetField(Tiff_id, TIFFTAG_IMAGEWIDTH, image_width / spp);
        TIFFGetField(Tiff_id, TIFFTAG_IMAGELENGTH, image_height);

        for (y = 0; y < image_height; y++) //écriture d'une page numérotée num_page, ligne par ligne (y).
        {
            TIFFReadScanline(Tiff_id, &buffer2D[y * image_width], y, 0);
        }
        int nbPix_plan=num_page*dim*dim;
        for (y = 0; y < dim; y++)
        {
            size_t num_lgn=y*dim;
            for(x = 0; x < dim; x++)
            {
                imgTiff[num_lgn + x+nbPix_plan]=buffer2D[num_lgn + x];
            }
        }
    }
    delete[] buffer2D;
    TIFFClose(Tiff_id);
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
