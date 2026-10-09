
#include <algorithm>
#include <vector>
#include "math_functions.h"
//calculate median value of an array (used toi estime epsilon in sup_redon division
double median(std::vector<double> valeurs)  // copie volontaire : nth_element réordonne le vecteur
{
    if (valeurs.empty()) return 0.0;

    size_t indiceMilieu = valeurs.size() / 2;
    std::nth_element(valeurs.begin(), valeurs.begin() + indiceMilieu, valeurs.end());

    if (valeurs.size() % 2 == 1) {
        return valeurs[indiceMilieu];  // taille impaire : élément du milieu
    }

    double elementHaut = valeurs[indiceMilieu];
    std::nth_element(valeurs.begin(), valeurs.begin() + indiceMilieu - 1, valeurs.end());
    double elementBas = valeurs[indiceMilieu - 1];

    return (elementBas + elementHaut) / 2.0;  // taille paire : moyenne des deux éléments centraux
}
//calculate median value of an array, zero excluded
double medianeSansZeros(std::vector<double> const& poids)
{
    std::vector<double> poidsNonNuls;
    poidsNonNuls.reserve(poids.size());
    for (double valeur : poids) {
        if (valeur > 0.0) poidsNonNuls.push_back(valeur);
    }

    if (poidsNonNuls.empty()) return 0.0;

    size_t indiceMilieu = poidsNonNuls.size() / 2;
    std::nth_element(poidsNonNuls.begin(), poidsNonNuls.begin() + indiceMilieu, poidsNonNuls.end());

    if (poidsNonNuls.size() % 2 == 1) {
        return poidsNonNuls[indiceMilieu];
    }

    double elementHaut = poidsNonNuls[indiceMilieu];
    std::nth_element(poidsNonNuls.begin(), poidsNonNuls.begin() + indiceMilieu - 1, poidsNonNuls.end());
    double elementBas = poidsNonNuls[indiceMilieu - 1];

    return (elementBas + elementHaut) / 2.0;
}
