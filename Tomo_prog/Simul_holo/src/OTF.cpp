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
    this->b_Reflex=m1.b_Reflex;
}
/// coonstructor overload, allowing to override manip. b_reflex->dangerous ?
/*OTF::OTF(const  manip &my_m1, bool my_b_reflex)
    : manipOTF(my_m1), b_Reflex(my_b_reflex),Valeur(pow(my_m1.dim_final,3))
{
 int nbPix=pow(my_m1.dim_final,3);
    cout<<"OTf dimfinal="<<my_m1.dim_final<<endl;
    cout<<"nbpix="<<nbPix<<endl;
    for(size_t cpt=0;cpt<nbPix;cpt++)
    {
      Valeur[cpt].real(0);
      Valeur[cpt].imag(0);
    }
}*/


OTF::~OTF()
{
    //dtor
}
///Fill the 3D  OTF values. Need 2D spec values from a 2D scanning.
void OTF::retropropag(Point2D spec)
{
    int sign=1;
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
    if(manipOTF.b_Reflex==true){
        sign=-1;
    }
        for(kd.x=-Nmax; kd.x<Nmax; kd.x++){
            for(kd.y=-Nmax; kd.y<Nmax; kd.y++){

                //kd.z=round(sqrt(rcarre-(kd.x)*(kd.x)-(kd.y)*(kd.y)));
                if((kd.x*kd.x)+(kd.y*kd.y)<fmcarre){//le spectre est dans un disque de rayon NXMAX

                    kd.z=sqrt(rcarre-(kd.x)*(kd.x)-(kd.y)*(kd.y));
                        //kobj=kd-ki;
                    kobj.z=round(sign*kd.z-ki.z);
                    kobj.y=sign*kd.y-ki.y;
                    kobj.x=sign*kd.x-ki.x;

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


void OTF::symetrize_central()
{   int dimcarre=manipOTF.dim_final*manipOTF.dim_final;

    for(int cpt=0;cpt<Valeur.size();cpt++){
    int zi=cpt/(dimcarre), cpt2D=cpt-zi*dimcarre, yi=cpt2D/manipOTF.dim_final,xi=cpt2D%manipOTF.dim_final;

   // cout<<"("<<xi<<","<<yi<<","<<zi<<")"<<endl;
    int xH=xi-manipOTF.dim_final/2, yH=yi-manipOTF.dim_final/2, zH=zi-manipOTF.dim_final/2;
   // Point3D k_orig(xH,yH, zH,manipOTF.dim_final);
    Point3D  k_sym(xH,yH,-zH,manipOTF.dim_final);
    int cpt3D_sym=k_sym.coordI().cpt3D();

    if(Valeur[cpt].real()!=0)
        Valeur[k_sym.coordI().cpt3D()].real(1);
        Valeur[k_sym.coordI().cpt3D()].imag(1);
    }
   // SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_sym__Re.bin",t_float,"wb");
   // SAV3D_Tiff(Valeur,"Re","/home/mat/tomo_test/otf_sym__Re.tif",1);
}

void OTF::symetrize_xoy()
{   int dimcarre=manipOTF.dim_final*manipOTF.dim_final;
        cout<<"dim final="<<manipOTF.dim_final<<endl;
    for(int cpt=0;cpt<Valeur.size();cpt++){
        int zi=cpt/(dimcarre), cpt2D=cpt-zi*dimcarre, yi=cpt2D/manipOTF.dim_final,xi=cpt2D%manipOTF.dim_final;///get back 3D coordinate in 3D computer space

        int xH=xi-manipOTF.dim_final/2, yH=yi-manipOTF.dim_final/2, zH=zi-manipOTF.dim_final/2;///convert computer space coordinate into "human" coordinate system
        Point3D  k_sym(xH,yH,-zH,manipOTF.dim_final);
        int cpt3D_sym=k_sym.coordI().cpt3D();///calculate the 3D counter in 1D vector

        if(Valeur[cpt].real()!=0){
            Valeur[k_sym.coordI().cpt3D()].real(1);
            Valeur[k_sym.coordI().cpt3D()].imag(1);
        }
   // SAVCplx(Valeur,"Re","/home/mat/tomo_test/otf_sym__Re.bin",t_float,"wb");
   // SAV3D_Tiff(Valeur,"Re","/home/mat/tomo_test/otf_sym__Re.tif",1);
    }
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
cout<<"nombre de spires="<<Nmax/(2*a*M_PI)<<endl;
  ofstream myfile;
  myfile.open ("spiral_uni_400.txt");
  for(double theta=0.01;theta<theta_m;theta=theta+theta_m/manipOTF.nbHolo){
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


void OTF::scan_uniform3D(vector<Point2D> &CoordSpec,float coef_limit)//coef  limit is used toi kimit the effective NA (scanning can never reach exactly the max NA and scanning at max NA can cause sampling problem)
{   const int Nxmax=manipOTF.NXMAX;
    const int nbHolo=manipOTF.nbHolo;
    // Validate input
    if (nbHolo <= 0 || Nxmax <= 0) {
        cerr << "Error: Invalid parameters nbHolo=" << nbHolo
             << " Nxmax=" << Nxmax << endl;
        return;
    }
    ///calculate sampling paramters
    const double surface_moyenne_pt=2*M_PI/nbHolo;
    const double distance_moyenne=sqrt(surface_moyenne_pt);
    const double theta_max=M_PI/2.0*coef_limit;//UDHS. if Nmax is limited by coef limit, Theta must also be limited by this coef.
    const double theta_min=0;//i
   // double theta_max=asin(Nxmax*coef_limit/manipOTF.R_EwaldPix);UDCS
    const int nbCercle=static_cast<int>(round((theta_max-theta_min)/distance_moyenne));//MTHETA=nbCercle concentrique
    cout<<"nb_total_cercle="<<nbCercle<<endl;
    double delta_theta=(theta_max-theta_min)/nbCercle;///angular sampling between circles
    double theta=0,phi=0,delta_curv_phi=0;

    //Pre-calculate total circle length for uniform sampling
    double total_circle_length=0;
    for(int num_cercle=0;num_cercle<nbCercle;num_cercle++){//scan polar half diameter
        theta=(num_cercle)*delta_theta;//theta=0, -> exclusion zéro fréquence
        delta_curv_phi=surface_moyenne_pt/delta_theta; //donc d ?
        total_circle_length=total_circle_length+2.0*M_PI*sin(theta);
    }

    const double dist_recalc=total_circle_length/nbHolo;
    cout<<"nb_cercle="<<nbCercle<<endl;
    size_t  nbTotalPoint=0, numHolo=0;
    Point2D spec_H(0,0,2*Nxmax);//spec_H is a 2D point in a space of size (2*NXMAx,2*NXMAX)
    vector<double> centres(4*Nxmax*Nxmax,0.0);
    for(int num_cercle=0;num_cercle<nbCercle;num_cercle++){//scan polar half diameter
        theta=(num_cercle)*delta_theta;//theta=0, singularité donc on exclut num_cercle==0
        //delta_phi=surface_moyenne_pt/delta_theta; //donc d ?
        delta_curv_phi=dist_recalc;//delta abscisse curviligne sur le cercle concentrique
         size_t nbPtPhi=round(2.0*M_PI*sin(theta)/(delta_curv_phi));//M_PHI=nbPoint dans le cercle concentrique actuel=perimetre/delta_phi
         nbTotalPoint+=nbPtPhi;
         if(nbTotalPoint>nbHolo){///avoid overflow due to rounding, which cause table index overflow and segfault (numHolo>nbHolo->index overflox in Coordspec)
            nbPtPhi=nbPtPhi-(nbTotalPoint-nbHolo);
         cout<<"nbTotalPoint="<<nbTotalPoint<<" exceeding nbHolo : removing "<<nbTotalPoint-nbHolo<<" point(s) in the biggest circle"<<endl;
         }


        for(int numPt=0;numPt<nbPtPhi;numPt++){//scan  concentric circles
           phi=2.0*M_PI*(numPt)/nbPtPhi;//sampling circle with phi
           spec_H.x=coef_limit*Nxmax*sin(theta)*cos(phi);
           spec_H.y=coef_limit*Nxmax*sin(theta)*sin(phi);
           CoordSpec[numHolo].x=round(spec_H.x);
           CoordSpec[numHolo].y=round(spec_H.y);
           centres[spec_H.coordI().cpt2D()]=numHolo+1;
           retropropag(spec_H);
           numHolo++;
        }
    }
///0 frequency specular
CoordSpec[0].x=0;
CoordSpec[0].y=0;
centres[CoordSpec[0].coordI().cpt2D()]=1;
retropropag(CoordSpec[0]);
SAV_Tiff2D(centres,manipOTF.chemin_result+"centres_dans_otf_uni3d.tif",manipOTF.Tp_holo);

}
