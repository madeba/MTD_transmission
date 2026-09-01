#include "regularized_division.h"
#include <assert.h>
using namespace std;
vector<double> gradU_U_dampedDiv(vector<complex<double>> &UBorn,vector<complex<double>> &gradUBorn, double alpha){

    size_t nbPix=UBorn.size();
    size_t dim=sqrt(nbPix);
    // for (auto& u : UBorn) max_norm2 = std::max(max_norm2, std::norm(u));
    double alpha_global=alpha;
   // int margin = 2;//exclude 2 pixels on the image borders.

    vector<double> result(nbPix);

    double beta_local=0.1;
    //result=global_damp(gradUBorn,UBorn, alpha_global);
    result=local_damp(gradUBorn,UBorn, alpha_global,beta_local);
    return result;

}

//global damp, based on max or |UBorn|^2
vector<double> global_damp(vector<complex<double>> &gradUBorn, vector<complex<double>> UBorn, double alpha_global){
    double max_norm2 = 0.0;
    size_t nbPix=UBorn.size();
    size_t dim=sqrt(nbPix);
    for (const auto& ChpCplx : UBorn)//find max |Uborn|
        max_norm2 = std::max(max_norm2, std::norm(ChpCplx));

    double eps2 = alpha_global * max_norm2; // ε² = (alpha * max|U|)²
    vector<double> result(UBorn.size());
    size_t margin =2;
    int mid = dim / 2;
        for (size_t cpt = 0; cpt < nbPix; ++cpt) {
        size_t x = cpt % dim;
        size_t y = cpt / dim;

//exclure les bords et la croix centrale dans l'espace symétrisé
//exclude borders which can have artifact due to fft windows
if (x < margin || x >= dim - margin ||
    y < margin || y >= dim - margin ||
    std::abs((int)x - mid) < margin ||
    std::abs((int)y - mid) < margin)
        {
        result[cpt] = 0.0;
        }
        else{
                        result[cpt] = imag(gradUBorn[cpt] * std::conj(UBorn[cpt]) / (std::norm(UBorn[cpt]) + eps2));
        }
    }
    return result;
}

//global AND local damp, based on max or |UBorn|^2 and local UBorn^2. this one should be prefered in  gradU_U_dampedDiv
vector<double> local_damp(vector<complex<double>> &gradUBorn, vector<complex<double>> UBorn, double alpha_global, double beta_local)
{
      double max_norm2 = 0.0;
      vector<double> result(UBorn.size());
      size_t nbPix=UBorn.size();
      size_t dim=sqrt(nbPix);
      size_t margin =2;
      int mid = dim / 2;
      for (const auto& ChpCplx : UBorn){//find max |Uborn|
        max_norm2 = std::max(max_norm2, std::norm(ChpCplx));
      }

     double eps2 = alpha_global * max_norm2; // ε² = (alpha * max|U|)²
            for (size_t cpt = 0; cpt < nbPix; ++cpt) {
        size_t x = cpt % dim;
        size_t y = cpt / dim;
        if (x < margin || x >= dim - margin ||
        y < margin || y >= dim - margin ||
        std::abs((int)x - mid) < margin ||
        std::abs((int)y - mid) < margin)
        {
        result[cpt] = 0.0;
        }
        else{

          double norm2 = std::norm(UBorn[cpt]);
          double eps2_local = beta_local * norm2 + alpha_global * max_norm2;
          double denom = norm2 + eps2_local;
          //pondération à la fois locale et globale
          result[cpt] =  imag(gradUBorn[cpt] * conj(UBorn[cpt])/denom);
          }
      }

      return result;
}
