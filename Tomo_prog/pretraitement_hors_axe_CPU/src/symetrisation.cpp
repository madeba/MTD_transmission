#include "symetrisation.h"
#include <assert.h>
using namespace std;
// symétrie miroir assurant la continuité fft aux bords
//mirror symetry, to ensure fft continuity at the borders
void Symetrise_mirror2(
    const std::vector<double>& img,
    std::vector<double>& imgSym)
{
    size_t dim = static_cast<size_t>(std::sqrt(img.size()));//eviter warning type
    size_t dim2 = 2 * dim;

    assert(dim * dim == img.size());
    assert(imgSym.size() == dim2 * dim2);

    for (size_t y = 0; y < dim; y++) {
        for (size_t x = 0; x < dim; x++) {

            double val = img[x + y * dim];

            // coordonnées miroir
            size_t xm = dim - 1 - x;
            size_t ym = dim - 1 - y;

            // 4 quadrants
            // bas droite (original)
            imgSym[(x + dim) + (y + dim) * dim2] = val;

            // bas gauche (miroir x)
            imgSym[xm + (y + dim) * dim2] = val;

            // haut droite (miroir y)
            imgSym[(x + dim) + ym * dim2] = val;

            // haut gauche (miroir x+y)
            imgSym[xm + ym * dim2] = val;
        }
    }
}
///overload for complex number
void Symetrise_mirror2(
    const std::vector<complex<double>>& img,
    std::vector<complex<double>>& imgSym)
{
    size_t dim = static_cast<size_t>(std::sqrt(img.size()));//eviter warning type
    size_t dim2 = 2 * dim;

    assert(dim * dim == img.size());
    assert(imgSym.size() == dim2 * dim2);

    for (size_t y = 0; y < dim; y++) {
        for (size_t x = 0; x < dim; x++) {

            complex<double> val = img[x + y * dim];

            // coordonnées miroir
            size_t xm = dim - 1 - x;
            size_t ym = dim - 1 - y;

            // 4 quadrants
            // bas droite (original)
            imgSym[(x + dim) + (y + dim) * dim2] = val;

            // bas gauche (miroir x)
            imgSym[xm + (y + dim) * dim2] = val;

            // haut droite (miroir y)
            imgSym[(x + dim) + ym * dim2] = val;

            // haut gauche (miroir x+y)
            imgSym[xm + ym * dim2] = val;
        }
    }
}

///"antisymétrie", uniquement pour symétriser un gradient avant integration par fft
///"antisymtry", only to integrate gradient
void Symetrise_mirror_gradient(
    const std::vector<double>& gradx,
    const std::vector<double>& grady,
    std::vector<double>& gradxSym,
    std::vector<double>& gradySym)
{
    size_t dim = static_cast<size_t>(std::sqrt(gradx.size()));
    size_t dim2 = 2 * dim;

    assert(dim * dim == gradx.size());
    assert(grady.size() == gradx.size());
    assert(gradxSym.size() == dim2 * dim2);
    assert(gradySym.size() == dim2 * dim2);

    for (size_t y = 0; y < dim; y++) {
        for (size_t x = 0; x < dim; x++) {

            size_t id = x + y * dim;

            double gx = gradx[id];
            double gy = grady[id];

            size_t xm = dim - 1 - x;
            size_t ym = dim - 1 - y;

            // =========================
            // 1. Bas droite (original)
            // =========================
            gradxSym[(x + dim) + (y + dim) * dim2] = gx;
            gradySym[(x + dim) + (y + dim) * dim2] = gy;

            // =========================
            // 2. Bas gauche (miroir X)
            // =========================
            gradxSym[xm + (y + dim) * dim2] = -gx; // inversion !
            gradySym[xm + (y + dim) * dim2] =  gy;

            // =========================
            // 3. Haut droite (miroir Y)
            // =========================
            gradxSym[(x + dim) + ym * dim2] =  gx;
            gradySym[(x + dim) + ym * dim2] = -gy; // inversion !

            // =========================
            // 4. Haut gauche (miroir X+Y)
            // =========================
            gradxSym[xm + ym * dim2] = -gx;
            gradySym[xm + ym * dim2] = -gy;
        }
    }
}
std::vector<double> cut_quad4(std::vector<double>& img4)
{
    size_t dim2 = static_cast<size_t>(std::sqrt(img4.size()));
    assert(dim2 * dim2 == img4.size());
    assert(dim2 % 2 == 0);

    size_t dim = dim2 / 2;

    std::vector<double> img(dim * dim);

    for (size_t y = 0; y < dim; y++) {
        for (size_t x = 0; x < dim; x++) {

            size_t src = (x + dim) + (y + dim) * dim2;
            size_t dst = x + y * dim;

            img[dst] = img4[src];
        }
    }

    return img;
}

//découpe le 4eme cadrant pour récupérer l'image dans une image symétrisée (surcharge complexe)
///ATTENTION! RETOURNE LA PARTIE REELLE
vector<double> cut_quad4(vector<complex<double>> const &monImg4Quad)
{
    size_t cpt_orig,cpt4_final;
    size_t dim=sqrt(monImg4Quad.size()/4);//on balaye avec les coef de la petite image finale
    vector<double> monImgCut(dim*dim);
    for(size_t x=0; x<dim; x++)
        for(size_t y=0; y<dim; y++)
            {
                cpt_orig=x+y*dim;
                cpt4_final=dim+x+(dim+y)*(2*dim);
                monImgCut[cpt_orig]=monImg4Quad[cpt4_final].real();//dxW(x,y) //image originale
            }
    return monImgCut;
}
//découpe le 4eme cadrant pour récupérer l'image dans une image symétrisée (surcharge complexe)
vector<complex<double>> cut_quad4_cplx(vector<complex<double>> const &monImg4Quad)
{
    size_t cpt_orig,cpt4_final;
    size_t dim=sqrt(monImg4Quad.size()/4);//on balaye avec les coef de la petite image finale
    vector<complex<double>> monImgCut(dim*dim);
    for(size_t x=0; x<dim; x++)
        for(size_t y=0; y<dim; y++)
            {
                cpt_orig=x+y*dim;
                cpt4_final=dim+x+(dim+y)*(2*dim);
                monImgCut[cpt_orig]=monImg4Quad[cpt4_final];
            }
    return monImgCut;
}
