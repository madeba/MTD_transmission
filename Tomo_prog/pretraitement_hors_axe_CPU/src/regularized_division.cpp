#include "regularized_division.h"
#include <assert.h>
using namespace std;
vector<double> gradU_U_dampedDiv(vector<complex<double>> &UBorn,vector<complex<double>> &gradUBorn, double alpha){

    size_t nbPix=UBorn.size();
    // for (auto& u : UBorn) max_norm2 = std::max(max_norm2, std::norm(u));
    double alpha_global=alpha;
   // int margin = 2;//exclude 2 pixels on the image borders.

    vector<double> result(nbPix);

    double beta_local=0.1;
    //result=global_damp(gradUBorn,UBorn, alpha_global);
    int margin=2;
    result=local_damp(gradUBorn,UBorn, alpha_global,beta_local,margin,true);

    return result;

}

//global damp, based on max or |UBorn|^2 : it's only "local damp" with beta=0

vector <double> global_damp(vector<complex<double>> &UBorn,vector<complex<double>> &gradUBorn, double alpha_global){
    return local_damp(gradUBorn, UBorn, alpha_global, 0.0, 2, true);
}
// Division régularisée : Im( grad(U) * conj(U) / (|U|^2 + eps^2) )
// Prévue pour l'image symétrisée 2N x 2N (dim = 2N).
// eps^2 = alpha_global * max|U|^2 / (1 + beta_local)
//   -> même seuil global que dans la version d'origine, mais sans le facteur
//      1/(1+beta) sur l'amplitude du gradient de phase.
// Le masque (bords + jointure centrale) est symétrique par rapport au miroir
// (x <-> dim-1-x), ce qui préserve l'imparité exacte du champ de gradient.
vector<double> local_damp(const vector<complex<double>> &gradUBorn,
                          const vector<complex<double>> &UBorn,
                          double alpha_global,
                          double beta_local ,
                          int margin,
                          bool masque_jointure)
{
    const size_t nbPix = UBorn.size();
    const int dim = static_cast<int>(std::lround(std::sqrt(static_cast<double>(nbPix))));
    assert(static_cast<size_t>(dim) * dim == nbPix);
    assert(gradUBorn.size() == nbPix);

    const int mid = dim / 2;   // la jointure est entre mid-1 et mid

    double max_norm2 = 0.0;
    for (const auto &c : UBorn)
        max_norm2 = std::max(max_norm2, std::norm(c));

    const double eps2 = alpha_global * max_norm2 / (1.0 + beta_local);

    vector<double> result(nbPix, 0.0);

    for (size_t cpt = 0; cpt < nbPix; ++cpt) {
        const int x = static_cast<int>(cpt % dim);
        const int y = static_cast<int>(cpt / dim);

        // bords : pixels 0..margin-1 et dim-margin..dim-1 (appariés par le miroir)
        const bool bord = (x < margin) || (x >= dim - margin) ||
                          (y < margin) || (y >= dim - margin);

        // jointure : pixels mid-margin .. mid+margin-1 (symétrique autour de mid-0.5)
        const bool jointure = masque_jointure &&
                              ((x >= mid - margin && x < mid + margin) ||
                               (y >= mid - margin && y < mid + margin));

        if (bord || jointure) {
            result[cpt] = 0.0;
        } else {
            const double denom = std::norm(UBorn[cpt]) + eps2;
            result[cpt] = std::imag(gradUBorn[cpt] * std::conj(UBorn[cpt])) / denom;
        }
    }
    return result;
}
//global AND local damp, based on max or |UBorn|^2 and local UBorn^2. this one should be prefered in  gradU_U_dampedDiv
/*vector<double> local_damp(vector<complex<double>> &gradUBorn, vector<complex<double>> UBorn, double alpha_global, double beta_local)
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
*/
