
#include "math_functions.h"
#include <cmath>
using namespace std;



// percentile dans [0,1] (ex: 0.05 pour le 5e percentile)
// retourne la valeur de visibilité en dessous de laquelle se trouve cette fraction des hologrammes
double return_percentile_value(std::vector<double> values, double percentile)
{
    std::sort(values.begin(), values.end());

    double position = percentile * (values.size() - 1);
    size_t idx_bas = static_cast<size_t>(std::floor(position));
    size_t idx_haut = static_cast<size_t>(std::ceil(position));

    if (idx_bas == idx_haut)
        return values[idx_bas];
       // interpolation linéaire de la valeur (un peu exagérée ? )
    double frac = position - idx_bas;
    return values[idx_bas] * (1.0 - frac) + values[idx_haut] * frac;
}
// MAD (Median Absolute Deviation) : mesure robuste de dispersion autour de la médiane.
// Remplit med et mad par référence ; values est passé par valeur (copié), donc
// le tri interne ne modifie pas le tableau de l'appelant.
//MAD = médiane des valeurs absolues des écarts à la médiane.
//Formellement : MAD = median(|xᵢ − median(x)|)
void computeMedianAndMAD(std::vector<double> values, double &med, double &mad)
{
    std::sort(values.begin(), values.end());
    med = values[values.size() / 2];

    std::vector<double> absdev(values.size());
    for (size_t i = 0; i < values.size(); i++)
        absdev[i] = std::abs(values[i] - med);

    std::sort(absdev.begin(), absdev.end());
    mad = absdev[absdev.size() / 2];
}

// seuil de détection des hologrammes à faible visibilité, basé sur la médiane et la MAD
// facteur_mad : typiquement 3.0 (large tolérance) à 2.0 (plus strict)
//Plus facteur_mad est grand, plus on s'autorise à être loin de la médiane avant de considérer une valeur comme anormale
//— donc plus le seuil est permissif (moins d'hologrammes exclus).
double computeThresholdMAD(std::vector<double> values, double facteur_mad)
{
    double med, mad;
    computeMedianAndMAD(values, med, mad);
    return std::max(0.0, med - facteur_mad * mad);
}
