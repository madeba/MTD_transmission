#include "OTF.h"
#include "fonctions.h"
#include <fstream>

using namespace std;
///class OTF : use to generate different OTF, corresponding to different scanning pattern : rosace ("fleur"), spiral, Annular etc.

//OTF::OTF(manip m1):manipOTF(m1.dim_final),Obj3D::Obj3D(m1.dim_final)
//constructor : init 3D OTF.
OTF::OTF(manip m1):manipOTF(m1.dimROI_Cam),Valeur(pow(m1.dim_final,3))
{
    int nbPix=pow(m1.dim_final,3);
    cout<<"OTf dimfinal="<<m1.dim_final<<endl;
    cout<<"nbpix="<<nbPix<<endl;
    for(size_t cpt=0;cpt<nbPix;cpt++)
    {
      Valeur[cpt].real(0);
      Valeur[cpt].imag(0);
    }
}

OTF::~OTF()
{
    //dtor
}
///Fill the 3D  OTF values. Need 2D spec values from a 2D scanning.
void OTF::retropropag(Point2D spec)
{
    int Nmax=manipOTF.NXMAX;
    int dim_final=manipOTF.dim_final;
    double fmcarre=pow(Nmax,2);
   // cout<<"fmcarre="<<fmcarre<<endl;
    double rcarre=pow(manipOTF.R_EwaldPix,2);
    //cout<<"rcarre"<<rcarre<<endl;
    if(rcarre-spec.x*spec.x-spec.y*spec.y<0)
          cout<<"problème : ki imaginaire"<<endl;
   // Point3D ki(spec,round(sqrt(rcarre-spec.x*spec.x-spec.y*spec.y)),dim_final);
    Point3D ki(spec,sqrt(rcarre-spec.x*spec.x-spec.y*spec.y),dim_final);
    Point3D kobj(0,0,0,dim_final);
    Point3D kd(0,0,0,dim_final);//dim espace erronée mais sinon problème soustraction

        for(kd.x=-Nmax; kd.x<Nmax; kd.x++){
            for(kd.y=-Nmax; kd.y<Nmax; kd.y++){

                //kd.z=round(sqrt(rcarre-(kd.x)*(kd.x)-(kd.y)*(kd.y)));
                if((kd.x*kd.x)+(kd.y*kd.y)<fmcarre){//le spectre est dans un disque de rayon NXMAX
                    kd.z=sqrt(rcarre-(kd.x)*(kd.x)-(kd.y)*(kd.y));
                        //kobj=kd-ki;
                    kobj.z=round(kd.z-ki.z);
                    kobj.y=kd.y-ki.y;
                    kobj.x=kd.x-ki.x;

                    if(Valeur[kobj.coordI().cpt3D()].real()==0){
                    Valeur[kobj.coordI().cpt3D()].real(1);
                    Valeur[kobj.coordI().cpt3D()].imag(1);
                    nbPixEff++;
                    }
                    else{

                       // Valeur[kobj.coordI().cpt3D()].real(1+Valeur[kobj.coordI().cpt3D()].real());//supredon
                        nbPixRedon++;
                    }
                }
            }
        }
}


void OTF::symetrize_xoy()
{   int dimcarre=manipOTF.dim_final*manipOTF.dim_final;

    for(int cpt=0;cpt<Valeur.size();cpt++){
    int zi=cpt/(dimcarre), cpt2D=cpt-zi*dimcarre, yi=cpt2D/manipOTF.dim_final,xi=cpt2D%manipOTF.dim_final;

   // cout<<"("<<xi<<","<<yi<<","<<zi<<")"<<endl;

    Point3D k_orig(xi-manipOTF.dim_final/2,yi-manipOTF.dim_final/2,zi-manipOTF.dim_final/2,manipOTF.dim_final),
    k_sym(xi-manipOTF.dim_final/2,yi-manipOTF.dim_final/2,-zi+manipOTF.dim_final/2,manipOTF.dim_final);
    int cpt3D=k_sym.coordI().cpt3D();

    if(cpt3D>Valeur.size())
     {
    /* cout<<"cpt="<<cpt<<endl;
     cout<<"(xi,yi,zi)=("<<xi<<","<<yi<<","<<zi<<")"<<endl;
     cout<<"(x,y,z)=("<<xi-dim_finale3D/2<<","<<yi<<","<<zi<<")"<<endl;
     cout<<"(xsym,ysym,zsym)=("<<xi-dim_finale3D/2<<","<<yi-dim_finale3D/2<<","<<-zi+dim_finale3D/2<<")"<<endl;
 cout<<"cpt3D="<<cpt3D<<endl;*/
    }
    if(Valeur[cpt].real()!=0)
        Valeur[k_sym.coordI().cpt3D()].real(1);
    }
   // SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_sym__Re.bin",t_float,"wb");
    SAV3D_Tiff(Valeur,"Re","/home/mat/tomo_test/otf_sym__Re.tif",1);
}


void OTF::bFermat(int nbHolo)
{
    double theta=0,theta_max=nbHolo;//number of holograms = theta max in radian !!
    double angle_dor=(3-sqrt(5))*M_PI;//nombre d'or
    double delta_theta=theta_max/nbHolo;
    int Nmax_cond=manipOTF.NXMAX_cond, dim_Uborn=manipOTF.dim_Uborn;///set-up paramaters
    vector<double> centre(dim_Uborn*dim_Uborn,0);///save the 2D pattern
    Point2D spec(0,0,dim_Uborn);

    for(int numHolo=1; numHolo<=nbHolo; numHolo++)      ///1) génération de l'image du pattern dans centre[]
    {//on démarre à 1 pour distinguer le point zéro du fond égal à 0
        spec.x=sqrt((theta)/theta_max)*round(dim_Uborn/2)*cos((theta+0.5)*angle_dor);
        spec.y=sqrt((theta)/theta_max)*round(dim_Uborn/2)*sin((theta+0.5)*angle_dor);
        centre[spec.coordI().cpt2D()]=1;
        retropropag(spec);
        theta=theta+delta_theta;
    }
}
///Annular scanning with multiple circles
void OTF::bMultiCercleUNI(int nb_cercle)
{
    int Nmax_cond=manipOTF.NXMAX_cond;
    int dim_Uborn=manipOTF.dim_Uborn;
    double rcarre=pow(Nmax_cond,2);
    vector<double> centre(dim_Uborn*dim_Uborn,0);
    Point2D spec(0,0,dim_Uborn);
    double longTot=0;
    int nbSpec=0;

    //calculer la longueur totale des périmètres des cercles
    for(double R_cercle=0; R_cercle<Nmax_cond+1; R_cercle=R_cercle+Nmax_cond/nb_cercle)
    {
        cout<<"R_cercle="<<R_cercle<<endl;
        longTot=longTot+2*M_PI*R_cercle;
    }
    cout<<"longueur totale="<<longTot<<endl;
    cout<<"Nmax_cond="<<Nmax_cond<<endl;
    for(double R_cercle=Nmax_cond-5; R_cercle>0; R_cercle=R_cercle-Nmax_cond/nb_cercle)
    {
        double perimetre=2*M_PI*R_cercle;

        int nb_ki=round(manipOTF.nbHolo*perimetre/longTot);
        cout<<"nb_ki=="<<nb_ki<<endl;
        for(double theta=0; theta<=2*M_PI; theta=theta+2*M_PI/nb_ki)
        {
            spec.x=round(R_cercle*cos(theta));
            spec.y=round((R_cercle*sin(theta)));
            if(spec.x*spec.x+spec.y*spec.y<rcarre)
            {
                centre[spec.coordI().cpt2D()]=1;
                retropropag(spec);
                nbSpec++;
            }
        }
    }
    spec.x=0,spec.y=0;
    centre[spec.coordI().cpt2D()]=1;
    retropropag(spec);

//SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_Re.bin",t_float,"wb");
    SAV2(centre,"/home/mat/tomo_test/centre.bin",t_float,"wb");
}
///annular scanning with one circle and one paramter to limit the scanning NA (100=100%=max NA)
void OTF::bCercle(int pourcentage_NA)
{
    int Nmax=manipOTF.NXMAX;
    int dim_Uborn=manipOTF.dim_Uborn;
    double rcarre=Nmax*Nmax;
    vector<double> centre(dim_Uborn*dim_Uborn,0);
    Point2D spec(0,0,dim_Uborn);
    int nbSpec=0;
    double R_cercle=round(Nmax*pourcentage_NA/100);
    double perimetre=2*M_PI*R_cercle;


    for(double theta=0;theta<=2*M_PI;theta=theta+2*M_PI/manipOTF.nbHolo)
       {
        spec.x=round(R_cercle*cos(theta));
        spec.y=round((R_cercle*sin(theta)));
        if(spec.x*spec.x+spec.y*spec.y<rcarre){
        centre[spec.coordI().cpt2D()]=1;
        retropropag(spec);
        nbSpec++;
        }
       }
SAV3D_Tiff(Valeur,"Re","/home/mat/tomo_test/otf.tif",1);
//SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_Re.bin",t_float,"wb");
SAV2(centre,"/home/mat/tomo_test/centre.bin",t_float,"wb");
}

///double spiral
void OTF::bDblSpiral()
{
    int Nmax=manipOTF.NXMAX;
    int dim_Uborn=manipOTF.dim_Uborn;
 vector<double> centre(dim_Uborn*dim_Uborn,0);
 double a=4,rho=0;///amplification du rayon polaire rho par raport à l'angle, rayon polaire
 Point2D spec(0,0,dim_Uborn),spec_anti(0,0,dim_Uborn);
 double rayon_float=Nmax;
 double theta_m=Nmax/a;
 double rcarre=Nmax*Nmax;
 double L=a/2*(log(theta_m+sqrt(theta_m*theta_m+1))+theta_m*sqrt(theta_m*theta_m+1));
int nbSpec=0;

 // for(double theta=0;theta<=theta_m;theta=theta+12*M_PI/(round(nbHolo/2)))
 for(double theta=0;theta<=theta_m;theta=theta+2*L/(a*manipOTF.nbHolo*sqrt(theta*theta+1)))
  {
      rho=a*theta;
     // cout<<"rho="<<rho<<endl;
      spec.x=round(rho*cos(theta));
      spec.y=round((rho*sin(theta)));
      spec_anti.x=-round(rho*cos(theta));
      spec_anti.y=-round(rho*sin(theta));

    if(spec.x*spec.x+spec.y*spec.y<rcarre){
        centre[spec.coordI().cpt2D()]=1;
        centre[spec_anti.coordI().cpt2D()]=1;
        retropropag(spec);
        retropropag(spec_anti);
        nbSpec++;
    }
  }
//cout<<"nbspec="<<2*nbSpec<<endl;
SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_Re.bin",t_float,"wb");
SAV2(centre,"/home/mat/tomo_test/centre.bin",t_float,"wb");
}
///simple archimede spiral; point are uniformly distributed in curvilinear abscissa
void OTF::bSpiral(){
    int Nmax=manipOTF.NXMAX;
    int dim_Uborn=manipOTF.dim_Uborn;
 vector<double> centre(dim_Uborn*dim_Uborn,0);
 double a=2,rho=0,theta_m=Nmax/a;///amplification du rayon polaire rho par rapport à l'angle, rayon polaire
 Point2D spec(0,0,dim_Uborn);
double fmcarre=Nmax*Nmax;
 int nbSpec=0;

 double L=a/2*(log(theta_m+sqrt(theta_m*theta_m+1))+theta_m*sqrt(theta_m*theta_m+1));
  //for(double theta=0;theta<rayon/a;theta=theta+Nmax/(a*nbHolo)){
cout<<"Longueur spirale="<<L<<endl;
cout<<"nombre de spires="<<Nmax/(2*a*M_PI)<<endl;
nbPixEff=0;
nbPixRedon=0;
  for(double theta=0.01;theta<theta_m;theta=theta+L/(a*manipOTF.nbHolo*sqrt(theta*theta+1))){
   // cout<<"theta="<<theta<<endl;
   // cout<<"delta_theta="<<Nmax*Nmax/(a*a*nbHolo*theta)<<endl;
    rho=a*theta;
    spec.x=round(rho*cos(theta));
    spec.y=round((rho*sin(theta)));

    if(spec.x*spec.x+spec.y*spec.y<fmcarre){
        centre[spec.coordI().cpt2D()]=1;
        retropropag(spec);
        nbSpec++;
        }
  }

cout<<"nbSpec="<<nbSpec<<endl;
SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_Re.bin",t_float,"wb");
SAV2(centre,"/home/mat/tomo_test/centre.bin",t_float,"wb");
}

//spirale uniforme en angle, donc non uniforme en chemin curviligne
void OTF::bSpiralNU(){
int Nmax=manipOTF.NXMAX;
    int dim_Uborn=manipOTF.dim_Uborn;
 vector<double> centre(dim_Uborn*dim_Uborn,0);
 double a=2,rho=0,theta_m=Nmax/a;///amplification du rayon polaire rho par raport à l'angle, rayon polaire
 Point2D spec(0,0,dim_Uborn);
 double rcarre=Nmax*Nmax;
 int nbSpec=0;

// double L=a/2*(log(theta_m+sqrt(theta_m*theta_m+1))+theta_m*sqrt(theta_m*theta_m+1));
  //for(double theta=0;theta<rayon/a;theta=theta+Nmax/(a*nbHolo)){
//cout<<"Longueur spirale="<<L<<endl;
cout<<"nombre de spires="<<Nmax/(2*a*M_PI)<<endl;
  ofstream myfile;
  myfile.open ("spiral_uni_400.txt");
  for(double theta=0.01;theta<theta_m;theta=theta+theta_m/manipOTF.nbHolo){
   // cout<<"theta="<<theta<<endl;
   // cout<<"delta_theta="<<Nmax*Nmax/(a*a*nbHolo*theta)<<endl;
    rho=a*theta;
    spec.x=(rho*cos(theta));
    spec.y=((rho*sin(theta)));

    if(spec.x*spec.x+spec.y*spec.y<rcarre){
    myfile<<"x "<<spec.x<<", y "<<spec.y<<endl;
        centre[spec.coordI().cpt2D()]=nbSpec;
        retropropag(spec);
        nbSpec++;
        }
  }
myfile.close();
//cout<<"nbSpec="<<nbSpec<<endl;
//SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_Re.bin",t_float,"wb");
//SAV2(centre,"/home/mat/tomo_test/centre.bin",t_float,"wb");
}

/*void OTF::bFleur(){
int Nmax=manipOTF.NXMAX;
int dim_Uborn=manipOTF.dim_Uborn;
vector<double> centre(dim_Uborn*dim_Uborn,0);
size_t nb=4;///controle du nombre de branches
Point2D spec(0,0,dim_Uborn);
//cout<<"Nxmax====="<<Nmax<<endl;
double rcarre=Nmax*Nmax;
int nbSpec=0;

for(double theta=0;theta<2*M_PI;theta=theta+2*M_PI/manipOTF.nbHolo){
         spec.x=round(Nmax*cos(nb*theta)*cos(theta));
         spec.y=round(Nmax*cos(nb*theta)*sin(theta));//arrondi trop tot?

        if(spec.x*spec.x+spec.y*spec.y<=rcarre){
            nbSpec++;
            centre[spec.coordI().cpt2D()]=1;
            retropropag(spec);
        }
    }
SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_Re.bin",t_float,"wb");
SAV2(centre,"/home/mat/tomo_test/centre.bin",t_float,"wb");
//cout<<"nbspec="<<nbSpec<<endl;
}*/
/*
vector<Var2D> OTF::bFleur(){
int Nmax=manipOTF.NXMAX;
int dim_Uborn=manipOTF.dim_Uborn;
vector<double> centre(dim_Uborn*dim_Uborn,0);

vector<Var2D> CoordSpec(manipOTF.nbHolo);
size_t nb=4;///controle du nombre de branches
Point2D spec(0,0,dim_Uborn);

//cout<<"Nxmax====="<<Nmax<<endl;
double rcarre=Nmax*Nmax;
int nbSpec=0;
short unsigned int num_holo=0;
for(double theta=0;theta<2*M_PI;theta=theta+2*M_PI/manipOTF.nbHolo){
         spec.x=round(Nmax*cos(nb*theta)*cos(theta));
         spec.y=round(Nmax*cos(nb*theta)*sin(theta));//arrondi trop tot?

        if(spec.x*spec.x+spec.y*spec.y<=rcarre){
            nbSpec++;
            centre[spec.coordI().cpt2D()]=1;
            CoordSpec[num_holo].x=spec.x;
            CoordSpec[num_holo].y=spec.y;
            retropropag(spec);
        }
        num_holo++;
    }
SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_Re.bin",t_float,"wb");
SAV2(centre,"/home/mat/tomo_test/centre.bin",t_float,"wb");
return CoordSpec;
//cout<<"nbspec="<<nbSpec<<endl;
}
*/
/*
vector<Var2D> OTF::bFleur(){
int Nmax=manipOTF.NXMAX;
int dim_Uborn=manipOTF.dim_Uborn;
vector<double> centre(dim_Uborn*dim_Uborn,0);

vector<Var2D> CoordSpec(manipOTF.nbHolo);

size_t nb=4;///controle du nombre de branches
Point2D spec(0,0,dim_Uborn);

//cout<<"Nxmax====="<<Nmax<<endl;
double rcarre=Nmax*Nmax;
int nbSpec=0;
short unsigned int num_holo=0;
for(double theta=0;theta<2*M_PI;theta=theta+2*M_PI/manipOTF.nbHolo){
        //cout<<"num_holo="<<num_holo<<endl;
      //  cout<<"theta="<<theta<<endl;
         spec.x=(int)round(Nmax*cos(nb*theta)*cos(theta));
         spec.y=(int)round(Nmax*cos(nb*theta)*sin(theta));//arrondi trop tot?

      //  if(spec.x*spec.x+spec.y*spec.y<=rcarre){
            nbSpec++;
            centre[spec.coordI().cpt2D()]=1;
            CoordSpec[num_holo].x=1;//(int)round(spec.x);
            CoordSpec[num_holo].y=1;//(int)round(spec.y);
            retropropag(spec);
            if(num_holo>98 && num_holo<110){
           //cout<<"num_holo="<<num_holo<<" : specOTF=("<<spec.x<<","<<spec.y<<")"<<endl;
           cout<<"num_holo="<<num_holo<<" : CoordOTF=("<<CoordSpec[num_holo].x<<","<<CoordSpec[num_holo].y<<")"<<endl;
            }
          //  cout<<"num_holoOTF="<<num_holo<<endl;
      //  }
        num_holo++;
    }

//SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_Re.bin",t_float,"wb");
SAV2(centre,manipOTF.chemin_result+"/centre_dans_otf.bin",t_float,"wb");
return CoordSpec;
//cout<<"nbspec="<<nbSpec<<endl;
}
*/


vector<Point2D> OTF::bFleur(short unsigned int const nbAxes){
int Nmax=manipOTF.NXMAX;
int dim_Uborn=manipOTF.dim_Uborn;
vector<double> centre(dim_Uborn*dim_Uborn,0);
Point2D ptInit(0,0,dim_Uborn);
vector<Point2D> CoordSpec(manipOTF.nbHolo,ptInit);
//vector<Point2D> *CoordSpec2=new vector<Point2D>(manipOTF.nbHolo);
short unsigned int const nb=nbAxes;///controle du nombre de branches
Point2D spec(0,0,dim_Uborn);

//cout<<"Nxmax====="<<Nmax<<endl;
double rcarre=Nmax*Nmax;
int nbSpec=0;
 int num_holo=0;
// double const delta_theta=2*M_PI/(manipOTF.nbHolo);
//for(double theta=0;theta<2*M_PI;theta=theta+delta_theta){
 double const delta_theta=2*M_PI/(manipOTF.nbHolo);
for(double theta=0;theta<2*M_PI;theta=theta+delta_theta){
        //cout<<"num_holo="<<num_holo<<endl;
      //  cout<<"theta="<<theta<<endl;
         spec.x=(int)round(Nmax*cos(nb*theta)*cos(theta));
         spec.y=(int)round(Nmax*cos(nb*theta)*sin(theta));//arrondi trop tot?

        if(spec.x*spec.x+spec.y*spec.y<=rcarre){
            nbSpec++;
            //centre[spec.coordI().cpt2D()]=1;//erreur exportation Coorspec si on laisse cette ligne
            CoordSpec[num_holo].x=round(spec.x);
            CoordSpec[num_holo].y=round(spec.y);
          //  cout<<"num_holo"<< num_holo<<", "<<CoordSpec[num_holo].dim2D<<endl;
            retropropag(spec);///project 2D pattern into a 3D OTF
          /*  if(num_holo>20 && num_holo<30){
           cout<<"num_holo="<<num_holo<<" : specOTF=("<<spec.x<<","<<spec.y<<")"<<endl;
           cout<<" : CoordOTF=("<<CoordSpec[num_holo].x<<","<<CoordSpec[num_holo].y<<")"<<endl;
            }*/
          //  cout<<"num_holoOTF="<<num_holo<<endl;
        }
        num_holo++;
    }

//SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_Re.bin",t_float,"wb");
//SAV2(centre,manipOTF.chemin_result+"/centre_dans_otf.bin",t_float,"wb");
return CoordSpec;
//cout<<"nbspec="<<nbSpec<<endl;
}
/*
void OTF::bFleur(vector<Point2D> &CoordSpec){
int Nmax=manipOTF.NXMAX;
int dim_Uborn=manipOTF.dim_Uborn;
vector<double> centre(dim_Uborn*dim_Uborn,0);
Point2D ptInit(0,0,dim_Uborn);
//vector<Point2D> CoordSpec(manipOTF.nbHolo,ptInit);

size_t nb=4;///controle du nombre de branches
Point2D spec_H(0,0,dim_Uborn);

//cout<<"Nxmax====="<<Nmax<<endl;
double rcarre=Nmax*Nmax;
int nbSpec=0;
 int num_holo=0;
 int K_attenuation=1;
  double const delta_theta=2*M_PI/(K_attenuation*manipOTF.nbHolo);
for(double theta=0;theta<2*M_PI/K_attenuation;theta=theta+delta_theta){
//for(double theta=0;theta<2*M_PI;theta=theta+2*M_PI/manipOTF.nbHolo){
//cout<<"num_holo="<<num_holo<<endl;
//cout<<"delta_theta"<<delta_theta<<endl;
  //      cout<<"theta="<<theta<<endl;
         spec_H.x=(Nmax-1)*cos(nb*theta)*cos(theta);
         spec_H.y=(Nmax-1)*cos(nb*theta)*sin(theta);//arrondi trop tot?

       // if(spec_H.x*spec_H.x+spec_H.y*spec_H.y<=rcarre){
            nbSpec++;
            centre[spec_H.coordI().cpt2D()]=1;
            CoordSpec[num_holo].x=round(spec_H.x);
            CoordSpec[num_holo].y=round(spec_H.y);

            retropropag(spec_H);

         //  cout<<"num_holo="<<num_holo<<" : specOTF=("<<spec_H.x<<","<<spec_H.y<<")"<<endl;
          //cout<<" : CoordOTF=("<<CoordSpec[num_holo].x<<","<<CoordSpec[num_holo].y<<")"<<endl;

          //  cout<<"num_holoOTF="<<num_holo<<endl;
      //  }
        num_holo++;
    }

//SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_Re.bin",t_float,"wb");
//SAV2(centre,manipOTF.chemin_result+"/centre_dans_otf.bin",t_float,"wb");

//cout<<"nbspec="<<nbSpec<<endl;
}
*/

void OTF::bFleur(vector<Point2D> &CoordSpec, size_t const nbAxes){
size_t Nmax=manipOTF.NXMAX;
size_t dim_Uborn=manipOTF.dim_Uborn;
vector<double> centre(dim_Uborn*dim_Uborn,0);
Point2D ptInit(0,0,dim_Uborn);
Point2D spec_H(0,0,dim_Uborn);

double rcarre=Nmax*Nmax;
size_t num_holo=0;
int K_attenuation=1;
double theta=0;
double const delta_theta=2*M_PI/(K_attenuation*manipOTF.nbHolo);;

for(num_holo=0;num_holo<manipOTF.nbHolo;num_holo++){
    theta=delta_theta*num_holo;
    spec_H.x=(Nmax-1)*cos(nbAxes*theta)*cos(theta);
    spec_H.y=(Nmax-1)*cos(nbAxes*theta)*sin(theta);

    CoordSpec[num_holo].x=round(spec_H.x);
    CoordSpec[num_holo].y=round(spec_H.y);
    retropropag(spec_H);
    }
}
