#include "phase_functions.h"
//#include "projet.h"
using namespace std;

///calculate the wrapped phase from -pi to pi
void calcPhase_mpi_pi_atan2(vector<complex<double>> const &cplxField, vector<double> &phaseMod2pi){///calcul phase -PI-PI
for(size_t cpt=0;cpt<(size_t)cplxField.size();cpt++)
phaseMod2pi[cpt]=atan2(cplxField[cpt].imag(),cplxField[cpt].real());
}



