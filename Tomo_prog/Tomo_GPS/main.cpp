#include <iostream>
#include <fstream>
#include <vector>
#include "manip.h"
#include "IO_fonctions.h"
#include "FFT_fonctions.h"
#include <chrono>
#include "fonctions.h"
using namespace std;

int main()
{
  manip m1;
  ///Init complex refractive index
  vector<float> indice(readTiff3D(m1.chemin_result+"/indice.tif"));  //load refracive index
  vector<float> absorption(readTiff3D(m1.chemin_result+"/absorption.tif"));//load absorption
  int dim=round(pow(indice.size(),1.0/3.0));
  cout<<"dimension="<<dim<<endl;
  double quad_err=0;
  vector<complex<double>> indiceCplx(indice.size());//init complex refractive index
  vector<complex<double>> indiceCplx_previous(indice.size(),0.0);//init complex refractive index
  vector<complex<double>> SpectreIndiceCplx(indice.size());//init complex refractive index
  size_t nbPix3D=indice.size();

  string chemin_blanc_indice=m1.chemin_result+"/indice_blanc.tif";
  if(chemin_blanc_indice.c_str()){
   ifstream ifile;
   ifile.open(chemin_blanc_indice.c_str());
   if(ifile){
   cout<<"soustraction du blanc"<<endl;
   vector<float> indice_blanc(readTiff3D(chemin_blanc_indice.c_str()));//init complex refractive index
   for(size_t cpt=0;cpt<indice.size();cpt++) indice[cpt]=indice[cpt]-indice_blanc[cpt];
   }
   else       cout<<"Pas de  blanc pour l'indice. ";
  }



  vector<float> OTF3D(readTiff3D(m1.chemin_result+"/OTF3D.tif"));//init complex refractive index

  cout<<"indice.size()"<<indice.size()<<endl;
  ///measured refractive index
  for(size_t cpt=0;cpt<indice.size();cpt++){
    indiceCplx[cpt].real((double)indice[cpt]);
    indiceCplx[cpt].imag((double)absorption[cpt]);
  }
  vector<float>().swap(indice);//forcer la libération mémoire de sup_redon
  vector<float>().swap(absorption);//forcer la libération mémoire de sup_redon


  const unsigned int dim_final=dim;//;m1.dim_final mais lu dans le tiff d'origine
  Point3D dim3D(dim_final,dim_final,dim_final,dim_final);
  cout<<"dim_final="<<dim_final<<endl;
//   FFT_encaps tf3D_out(dim3D,m1.nbThreads);
  bool b_inplace=true;
  FFTW_init tf3D_INP(dim3D,6, b_inplace);

  ///calculate 'measured' spectrum
  FFT3Dcplx_Inplace(fftshift3D(indiceCplx), SpectreIndiceCplx,tf3D_INP);
  SpectreIndiceCplx=fftshift3D(SpectreIndiceCplx);

 // SAV3D_Tiff(SpectreIndiceCplx,"Re",m1.chemin_result+"SpectreIndiceCplxRE.tif",m1.tailleTheoPixelTomo);
  //SAV3D_Tiff(SpectreIndiceCplx,"Im",m1.chemin_result+"SpectreIndiceCplxIm.tif",m1.tailleTheoPixelTomo);
  vector<complex<double>> SpectreIndiceCplx_estim(nbPix3D);//init complex refractive index
  //for(cpt=0;cpt<nbPix3D;cpt++) SpectreIndiceCplx_estim[cpt]=SpectreIndiceCplx[cpt];


  cout<<"Itérations"<<endl;
  unsigned short int num_iter, nb_iter=m1.nbIterGPS;
  size_t cpt=0;

  for(num_iter=0;num_iter<nb_iter;num_iter++){
    auto start_part2= std::chrono::system_clock::now();
    //constraints on refractive index
     for(cpt=0;cpt<nbPix3D;cpt++){
        if(indiceCplx[cpt].real()>m1.delta_nMax) indiceCplx[cpt].real(m1.delta_nMax);//imposer la valeur MAX d'indice
            else if(indiceCplx[cpt].real()<m1.delta_nMin) indiceCplx[cpt].real(m1.delta_nMin);//imposer la valeur MIN d'indice

        if(indiceCplx[cpt].imag()<m1.kappa_Min) indiceCplx[cpt].imag(m1.kappa_Min);//imposer la valeur MIN de coefficient d'extinction
            else if(indiceCplx[cpt].imag()>m1.kappa_Max) indiceCplx[cpt].imag(m1.kappa_Max);//imposer la valeur MAX de coefficient d'extinction
        }

     FFT3Dcplx_Inplace(fftshift3D(indiceCplx), SpectreIndiceCplx_estim,tf3D_INP);
     SpectreIndiceCplx_estim=fftshift3D(SpectreIndiceCplx_estim);

     for(cpt=0;cpt<nbPix3D;cpt++)
       if(OTF3D[cpt]==1) SpectreIndiceCplx_estim[cpt]=SpectreIndiceCplx[cpt];//imposer les données spectrales expérimentales


     FFT3Dcplx_Inplace_INV(fftshift3D(SpectreIndiceCplx_estim),indiceCplx,tf3D_INP);//recalculer l'indice à partir du nouveau spectre, et recommencer
     indiceCplx=fftshift3D(indiceCplx);

 cout<<"quad err="<<calc_quad_err(indiceCplx,indiceCplx_previous)<<endl;

      indiceCplx_previous=indiceCplx;
     cout<<"num_iter="<<num_iter<<endl;
     auto end_part2= std::chrono::system_clock::now();
     auto elapsed_part2 = end_part2 - start_part2;
     std::cout <<"Temps pour 1 iter= "<< elapsed_part2.count()/(pow(10,9)) << '\n';
  }
  /// Remarque Fabien : on applique une dernière fois la positivité/ la contrainte
   for(cpt=0;cpt<nbPix3D;cpt++){
        if(indiceCplx[cpt].real()>m1.delta_nMax) indiceCplx[cpt].real(m1.delta_nMax);//imposer la valeur MAX d'indice
            else if(indiceCplx[cpt].real()<m1.delta_nMin) indiceCplx[cpt].real(m1.delta_nMin);//imposer la valeur MIN d'indice

        if(indiceCplx[cpt].imag()<m1.kappa_Min) indiceCplx[cpt].imag(m1.kappa_Min);//imposer la valeur MIN de coefficient d'extinction
            else if(indiceCplx[cpt].imag()>m1.kappa_Max) indiceCplx[cpt].imag(m1.kappa_Max);//imposer la valeur MAX de coefficient d'extinction
        }
cout<<"m1.tailleTheopixelTomo"<<m1.tailleTheoPixelTomo<<endl;
SAV3D_Tiff(SpectreIndiceCplx_estim,"Re",m1.chemin_result+"spectre_RE_GPS.tif",m1.tailleTheoPixelTomo*pow(10,-5));
SAV3D_Tiff(indiceCplx,"Re",m1.chemin_result+"indiceGPS.tif",m1.tailleTheoPixelTomo*pow(10,-5));
SAV3D_Tiff(indiceCplx,"Im",m1.chemin_result+"absorptionGPS.tif",m1.tailleTheoPixelTomo*pow(10,-5));

}
