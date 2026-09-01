#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <vector>
#include "struct.h"
#include "integration.h"
#include "vecteur.h"
#include "FFTW_init.h"
#include "FFT_fonctions.h"
#define M_2PI 2*M_PI
using namespace std;
using namespace Eigen;

// index column-major (comme MATLAB)
inline int idx(int x, int y, int ny) {
    return y + x * ny;
}

//integration de 2 gradients par moindre carré
VectorXd integrateGradientEigenLS(
    const MatrixXd& fx,
    const MatrixXd& fy,
    int nx,
    int ny)
{
    int N = nx * ny;
    int M = 2 * N;

    typedef Triplet<double> T;//un triplet a .col, .row() et .value(), ils servent à construire des matrices quasi vide (sparse)
    vector<T> triplets;
    VectorXd b = VectorXd::Zero(M);

    int row = 0;

    // --- X gradients ---matrice creuse qui encode l'opérateur gradient en X par différences finies.
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {

            int i = idx(x, y, ny);

            if (x == 0) {//bord gauche difference forward : ∂f/∂x ≈ f(x+1,y) - f(x,y)
                int i2 = idx(x+1, y, ny);
                triplets.emplace_back(row, i, -1);
                triplets.emplace_back(row, i2, 1);
            }
            else if (x == nx-1) { //bord droit difference backward ∂f/∂x ≈ f(x,y) - f(x-1,y)
                int i1 = idx(x-1, y, ny);
                triplets.emplace_back(row, i1, -1);
                triplets.emplace_back(row, i, 1);
            }
            else {//intérieur, différence centrée : ∂f/∂x ≈ (f(x+1,y) - f(x-1,y)) / 2
                int i1 = idx(x-1, y, ny);
                int i2 = idx(x+1, y, ny);
                triplets.emplace_back(row, i1, -0.5);
                triplets.emplace_back(row, i2, 0.5);
            }

            b(row) = fx(y,x);
            row++;
        }
    }

    // --- Y gradients ---
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {

            int i = idx(x, y, ny);

            if (y == 0) {
                int i2 = idx(x, y+1, ny);
                triplets.emplace_back(row, i, -1);
                triplets.emplace_back(row, i2, 1);
            }
            else if (y == ny-1) {
                int i1 = idx(x, y-1, ny);
                triplets.emplace_back(row, i1, -1);
                triplets.emplace_back(row, i, 1);
            }
            else {
                int i1 = idx(x, y-1, ny);
                int i2 = idx(x, y+1, ny);
                triplets.emplace_back(row, i1, -0.5);
                triplets.emplace_back(row, i2, 0.5);
            }

            b(row) = fy(y,x);
            row++;
        }
    }
    // --- Construire matrice sparse ---
    SparseMatrix<double> A(M, N); //matrice creuse M x N
    A.setFromTriplets(triplets.begin(), triplets.end());
    // --- Fixer constante (f(0)=0) ---
    // enlever colonne 0
    vector<T> triplets2;
    triplets2.reserve(triplets.size());
    for (const auto& t : triplets) {
        if (t.col() == 0) {
            b(t.row()) -= t.value() * 0.0;
        } else {
            triplets2.emplace_back(t.row(), t.col()-1, t.value());
        }
    }//Le système A * u = b a une infinité de solutions car le gradient ne détermine f qu'à une constante près.
    //Pour lever cette ambiguïté, on impose f(0) = 0.

    SparseMatrix<double> A2(M, N-1);
    A2.setFromTriplets(triplets2.begin(), triplets2.end());
//a nouvelle matrice A2 est de taille M x (N-1), et résoudre A2 * u2 = b donnera directement
//la solution avec la contrainte f(0) = 0 implicitement intégrée.

    // --- Solve LS ---


    LeastSquaresConjugateGradient<SparseMatrix<double>> solver;
    solver.setMaxIterations(1000);
       solver.setTolerance(1e-6);
    solver.compute(A2);

    VectorXd x = solver.solve(b);
//Le système A2 * x = b est sur-déterminé (plus d'équations que d'inconnues, car il y a une ligne par pixel pour chaque gradient).
//Il n'a donc pas de solution exacte en général.
//LeastSquaresConjugateGradient trouve la solution qui minimise ‖A2*x - b‖²,
//c'est-à-dire le champ f dont les gradients sont le plus proche possible des gradients cibles fx et fy.
//C'est la méthode itérative recommandée par Eigen pour les systèmes creux surdéterminés.

    // --- reconstruire f ---
    VectorXd f = VectorXd::Zero(N);//Reconstruction de f avec la contrainte f(0) = 0
    for (int i = 1; i < N; i++)
        f(i) = x(i-1);

    return f;
}


///gradient analytique de l'ellipse, pour tester l'intégration.

GradientResult generateEllipseWithGradient(
    Var2D centre,
    Var2D rayon,
    Var2D dim,
    double phase_max)
{
    GradientResult res;
    int N = dim.x * dim.y;

    res.phase.resize(N, 0.0);
    res.gradx.resize(N, 0.0);
    res.grady.resize(N, 0.0);

    double cx = dim.x / 2.0 + centre.x;
    double cy = dim.y / 2.0 + centre.y;

    double sigma_x = rayon.x / 4.0;
    double sigma_y = rayon.y / 4.0;

    double inv_sigma_x2 = 1.0 / (sigma_x * sigma_x);
    double inv_sigma_y2 = 1.0 / (sigma_y * sigma_y);

    for (int y = 0; y < dim.y; y++) {
        for (int x = 0; x < dim.x; x++) {

            int id = x + y * dim.x;
            double X = x - cx;
            double Y = y - cy;

            // masque ellipse
            double ellipse_eq = (X*X)/(rayon.x*rayon.x)
                              + (Y*Y)/(rayon.y*rayon.y);

            if (ellipse_eq <= 1.0) {

                double f = phase_max *
                    exp(-0.5 * (X*X*inv_sigma_x2 + Y*Y*inv_sigma_y2));

                res.phase[id] = f;

                // gradient analytique
                res.gradx[id] = f * (-X * inv_sigma_x2);
                res.grady[id] = f * (-Y * inv_sigma_y2);
            }
        }
    }

    return res;
}

///intégration du gradient par fft, entrées réelles
void integ_grad4(vector<double> const &gradx, vector<double> const& grady, vector<complex<double>> &sortie,std::vector<vecteur> &kvect_shift, FFTW_init &param_c2c)
{
  complex<double> I(0,1);
  size_t nbPix=gradx.size();
  vector<complex<double>> TF_gradx(nbPix), TF_grady(nbPix), tampon(nbPix);
  vector<double> kvect_mod_sq(nbPix);
  TF2Dcplx(gradx, TF_gradx, param_c2c);
  TF2Dcplx(grady, TF_grady, param_c2c);
  //SAV_Tiff2D(TF_gradx,"Re","/home/mat/tomo_test/TF_gradx.tiff",1);
  //SAV_vec3D(kvect_shift,"x","/home/mat/tomo_test/kvect_shift_x.bin","wb",nbPix);
  for(size_t cpt=0;cpt<nbPix;cpt++){
        double kx = kvect_shift[cpt].getx();
        double ky = kvect_shift[cpt].gety();
        // kvect_mod_sq[cpt]=pow(kvect_shift[cpt].getx(),2)+pow(kvect_shift[cpt].gety(),2);
        kvect_mod_sq[cpt]=kx*kx+ky*ky;
        if(kvect_mod_sq[cpt]>1e-12){
        tampon[cpt]=-I*(TF_gradx[cpt]*kx+TF_grady[cpt]*ky)/(M_2PI*kvect_mod_sq[cpt]);
        }
        else{
        tampon[cpt]=0;
              // cout<<"division par zéro="<<kvect_mod_sq[cpt]<<endl;
        }
  }
  TF2Dcplx_INV(tampon,sortie,param_c2c);
}

