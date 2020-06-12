#include <iostream>
#include <vector>
#include "manip.h"
#include "IO_fonctions.h"
#include "FFT_fonctions.h"
#include <chrono>
using namespace std;

int main()
{
    manip m1;
  ///Init complex refractive index
  vector<float> indice(readTiff3D(m1.chemin_result+"/indice.tif"));  //load refracive index
  vector<float> absorption(readTiff3D(m1.chemin_result+"/absorption.tif"));//load absorption
  vector<complex<double>> indiceCplx(indice.size());//init complex refractive index
  vector<complex<double>> SpectreIndiceCplx(indice.size());//init complex refractive index
  size_t nbPix3D=indice.size();


  vector<float> OTF3D(readTiff3D(m1.chemin_result+"/OTF3D.tif"));//init complex refractive index
  cout<<"indice.size()"<<indice.size()<<endl;
  ///measured refractive index
  for(size_t cpt=0;cpt<indice.size();cpt++){
    indiceCplx[cpt].real((double)indice[cpt]);
    indiceCplx[cpt].imag((double)absorption[cpt]);
  }
  vector<float>().swap(indice);//forcer la libération mémoire de sup_redon
  vector<float>().swap(absorption);//forcer la libération mémoire de sup_redon


  const unsigned int dim_final=512;//;m1.dim_final;
  Point3D dim3D(dim_final,dim_final,dim_final,dim_final);
  cout<<"dim_final="<<dim_final<<endl;
  // FFT_encaps tf3D_out(dim3D,m1.nbThreads);
  bool b_inplace=true;
  FFTW_init tf3D_INP(dim3D,m1.nbThreads, b_inplace);

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
        if(indiceCplx[cpt].real()>m1.delta_nMax) indiceCplx[cpt].real(m1.delta_nMax);
            else if(indiceCplx[cpt].real()<m1.delta_nMin) indiceCplx[cpt].real(m1.delta_nMin);
        if(indiceCplx[cpt].imag()<m1.kappa_Min) indiceCplx[cpt].imag(m1.delta_nMin);
            else if(indiceCplx[cpt].imag()>m1.kappa_Max) indiceCplx[cpt].real(m1.kappa_Max);
        }

     FFT3Dcplx_Inplace(fftshift3D(indiceCplx), SpectreIndiceCplx_estim,tf3D_INP);
     SpectreIndiceCplx_estim=fftshift3D(SpectreIndiceCplx_estim);

     for(cpt=0;cpt<nbPix3D;cpt++)
       if(OTF3D[cpt]==1) SpectreIndiceCplx_estim[cpt]=SpectreIndiceCplx[cpt];


     FFT3Dcplx_Inplace_INV(fftshift3D(SpectreIndiceCplx_estim),indiceCplx,tf3D_INP);
     indiceCplx=fftshift3D(indiceCplx);

     cout<<"num_iter="<<num_iter<<endl;
      auto end_part2= std::chrono::system_clock::now();
  auto elapsed_part2 = end_part2 - start_part2;
  std::cout <<"Temps pour 1 iter= "<< elapsed_part2.count()/(pow(10,9)) << '\n';
  }
SAV3D_Tiff(indiceCplx,"Re",m1.chemin_result+"indiceGPS.tif",m1.tailleTheoPixelTomo);
SAV3D_Tiff(indiceCplx,"Im",m1.chemin_result+"absorptionGPS.tif",m1.tailleTheoPixelTomo);

}
