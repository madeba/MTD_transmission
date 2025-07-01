#include "fonctions.h"
#define PI 3.14159
#include <chrono>
#include <vector>
#include <algorithm>
#include "FFT_encaps.h"
#include "FFT_fonctions.h"
/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */
using namespace std::chrono;
using namespace std;


vector<complex<double>> extractSliceZ(vector<complex<double>> &index3D, string axis, size_t height)
{
    int dim=cbrt(index3D.size());///cbrt calculate cube root
    size_t cpt2D,cpt3D;
    int nbPixPlan=dim*dim;//nb pixel in a plane, to avoid calculation inside the loop
    vector<complex<double>> slice2D(nbPixPlan);
    if(axis=="z"){
        for(int y=0;y<dim;y++)
        for(int x=0;x<dim;x++){
                cpt2D=x+y*dim;
                cpt3D=x+y*dim+height*nbPixPlan;
                slice2D[cpt2D]=index3D[cpt3D];
        }
    }
    else if(axis=="y"){
        for(int z=0;z<dim;z++)
        for(int x=0;x<dim;x++){
                cpt2D=x+z*dim;
                cpt3D=x+height*dim+z*nbPixPlan;
                slice2D[cpt2D]=index3D[cpt3D];
        }
    }
    else if(axis=="x"){
        for(int z=0;z<dim;z++)
        for(int y=0;y<dim;y++){
                cpt2D=y+z*dim;
                cpt3D=height+y*dim+z*nbPixPlan;
                slice2D[cpt2D]=index3D[cpt3D];
        }
    }
    else {
        throw invalid_argument("Invalid axis: must be 'x', 'y', or 'z'");
    }
    return slice2D;
}

string extract_string(std::string const &token,  std::string chemin_fic)
{
    ifstream fichier(chemin_fic.c_str(), ios::in);  // on ouvre en lecture
    string ligne,motcle,valeurMot,separ=" ";
    string valeur;
    vector<std::string> tokens;

    if(fichier)  // si l'ouverture a fonctionné
    {
        while(!fichier.eof())
        {
            getline(fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire
                tokens.push_back(ligne);
        }
    }
    else
        cerr << "Impossible d'ouvrir le fichier !"<< chemin_fic<< endl;

    size_t nb_tok=tokens.size();
    for(size_t cpt=0; cpt<nb_tok; cpt++)
    {
        ligne=tokens[cpt];
        if(ligne!="")
        {
            size_t pos_separ=ligne.find(separ);
            size_t long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin)
            if(motcle==token)
            {
                valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));
                cout<<motcle<<"="<<valeurMot<<endl;
                valeur=valeurMot.c_str();
            }
        }
    }
    if(valeur.empty())
        cout<<"mot_clé "<<token<<" inexistant dans le fichier "<<chemin_fic<<endl;
    fichier.close();
    return valeur;
}
float extract_val(string token,  string chemin_fic)
{
    ifstream fichier(chemin_fic.c_str(), ios::in);  // on ouvre en lecture
    string ligne,motcle,valeurMot,separ=" ";
    float valeur=0;
    vector<std::string> tokens;

    if(fichier)  // si l'ouverture a fonctionné
    {
        while(!fichier.eof())
        {
            getline(fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire
                tokens.push_back(ligne);
        }
    }
    else
        cerr << "Impossible d'ouvrir le fichier !"<< chemin_fic<< endl;

    int nb_tok=tokens.size();
    for(int cpt=0; cpt<nb_tok; cpt++)
    {
        ligne=tokens[cpt];
        if(ligne!="")
        {
            int pos_separ=ligne.find(separ);
            int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin)
            if(motcle==token)
            {
                valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));
                cout<<motcle<<"="<<valeurMot<<endl;
                valeur=atof(valeurMot.c_str());
            }
        }
    }
    if(valeurMot.empty())
        cout<<"mot_clé "<<token<<" inexistant dans le fichier "<<chemin_fic<<endl;
    fichier.close();
    return valeur;
}

float extract_val(string token,  string chemin_fic, double defaut)
{
    ifstream fichier(chemin_fic.c_str(), ios::in);  // on ouvre en lecture
    string ligne,motcle,valeurMot,separ=" ";
    float valeur=0;
    vector<std::string> tokens;

    if(fichier)  // si l'ouverture a fonctionné
    {
        while(!fichier.eof())
        {
            getline(fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire
                tokens.push_back(ligne);
        }
    }
    else
        cerr << "Impossible d'ouvrir le fichier !"<< chemin_fic<< endl;

    int nb_tok=tokens.size();
    for(int cpt=0; cpt<nb_tok; cpt++)
    {
        ligne=tokens[cpt];
        if(ligne!="")
        {
            int pos_separ=ligne.find(separ);
            int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin)
            if(motcle==token)
            {
                valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));
                cout<<motcle<<"="<<valeurMot<<endl;
                valeur=atof(valeurMot.c_str());
            }
        }
    }
    if(valeurMot.empty()){
        cout<<"mot_clé "<<token<<" inexistant dans le fichier "<<chemin_fic<<endl;
    cout<<"Utlisation valeur par défaut: "<<defaut<<endl;
    valeur=defaut;
    }
    fichier.close();
    return valeur;
}







void prepare_wisdom2D(Var2D dim, const char *chemin)
{
    fftw_plan_with_nthreads(4);
    int N=dim.x*dim.y;

    fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
    fftw_plan p;
    //Réservation memoire
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_BACKWARD, FFTW_EXHAUSTIVE);
    fftw_export_wisdom_to_filename(chemin);
    fftw_destroy_plan(p);
}
void prepare_wisdom3D(Var3D dim, char *chemin)
{
    fftw_plan_with_nthreads(4);
    int N=dim.x*dim.y*dim.z;

    fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
    fftw_plan p;
    //Réservation memoire
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p=fftw_plan_dft_3d( dim.x,  dim.y, dim.z,in, out,FFTW_BACKWARD, FFTW_EXHAUSTIVE);
    fftw_export_wisdom_to_filename(chemin);
    fftw_destroy_plan(p);
}
int coordSpec(nbCplx* TF_UBorn, double *TF_champMod,Var2D NMAX)
{
    int cpt_max=0;
    TF_champMod[0]=pow(TF_UBorn[0].Re,2)+pow(TF_UBorn[0].Im,2);

    for(int cpt=1; cpt<(4*NMAX.x*NMAX.y); cpt++)
    {
        TF_champMod[cpt]=pow(TF_UBorn[cpt].Re,2)+pow(TF_UBorn[cpt].Im,2);
        if(TF_champMod[cpt]>TF_champMod[cpt_max])
        {
            cpt_max=cpt;
        }
    }
    return cpt_max;
}


double calc_mediane(double entree[], size_t nb_elem)
{
    int pos_med=round(nb_elem/2);
    //float a[] = {9, 8, 7, 6, 5, 0, 1, 2.5, 3, 4,11,125, -1};
    std::nth_element(entree, entree + pos_med, entree + nb_elem);
    double mediane=entree[pos_med];
    // cout<<"mediane : "<<mediane<<endl;
    return mediane;
}



void calcPhase2pi(nbCplx* obj, Var2D taille,double* phaseMod2pi)///calcul phase -PI-PI
{
    ///calcul phase de 0 à 2 PI
    size_t image_size=taille.x*taille.y;

    for(size_t pixel=0; pixel<image_size; pixel++)
    {
        double cos_phase=obj[pixel].Re;
        double sin_phase=obj[pixel].Im;

        if(sin_phase>0)
        {
            if(cos_phase>0)
            {
                double argTanPhase=sin_phase/cos_phase;
                phaseMod2pi[pixel]=atan(argTanPhase);
            }
            else if(cos_phase<0)
            {
                double argTanPhase=sin_phase/cos_phase;
                phaseMod2pi[pixel]=atan(argTanPhase)+PI;
            }
            else if(cos_phase==0)
            {
                phaseMod2pi[pixel]=PI/2;
            }
        }
        else if(sin_phase<0)
        {
            if(cos_phase<0)
            {
                double argTanPhase=sin_phase/cos_phase;
                phaseMod2pi[pixel]=atan(argTanPhase)+PI;
            }
            else if(cos_phase>0)
            {
                double argTanPhase=sin_phase/cos_phase;
                phaseMod2pi[pixel]=atan(argTanPhase)+2*PI;
            }
            else if(cos_phase==0)
            {
                phaseMod2pi[pixel]=3*PI/2;
            }
        }
        else if(sin_phase==0)
        {
            if(cos_phase>0)
            {
                phaseMod2pi[pixel]=0;
            }
            else if(cos_phase<0)
            {
                phaseMod2pi[pixel]=PI;
            }

        }
    }
}

/*void calc_Uborn(nbCplx *TF_UBorn,nbCplx *UBorn,Var2D dim2DHA,Var2D PosSpec)
{
    Var2D recalUBorn={-PosSpec.x,-PosSpec.y},DecalU_Born={dim2DHA.x/2,dim2DHA.y/2};
    size_t NbPixU_Born=dim2DHA.x*dim2DHA.y;
    nbCplx TF_UBorn_I[NbPixU_Born];
    nbCplx UBorn_I[NbPixU_Born];
    decal2DCplxGen(TF_UBorn,TF_UBorn_I,dim2DHA,recalUBorn);
    TF2Dcplx_INV(TF_UBorn_I,UBorn_I,dim2DHA);
    circshift2DCplx(UBorn_I,UBorn,dim2DHA,DecalU_Born);

}*/

void Chrono(temps *t, string message)
{
    t->fin = clock ();
    float temps_cpu = (t->fin - t->init) * 1e-6;
    t->total= t->total+temps_cpu;
    // cout<<"coucou total="<<t->fin<<endl;
    t->init = clock();
    cout<<message<<endl;

}
int chargeBin(float *objet, string chemin,  int NbPix)
{
    size_t precision=sizeof(objet[0]);//attention si on n'indique pas "0", sizeof donnera la taille du pointeur (64bits)!
    FILE *fichier_ID = NULL;
    fichier_ID=fopen(chemin.c_str(),"rb");
    if(fichier_ID==NULL)
        cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;
    cout<<"precision en nb d'octet="<<precision<<endl;
    int isfileok=fread(objet,precision,NbPix,fichier_ID);
    if(isfileok==0)
        cout<<"erreur lecture : "<<chemin<<endl;
    fclose(fichier_ID);
    return 1;
}

/*
///recalage par correlation croisée
void recale(nbCplx* obj,nbCplx* objDecal,nbCplx *objRecal, Var3D dimVol)
{
    unsigned int Npix3D=dimVol.x*dimVol.y*dimVol.z, cpt=0;
    nbCplx *A=new nbCplx[Npix3D];
    nbCplx *B=new nbCplx[Npix3D];
    nbCplx *C=new nbCplx[Npix3D];
    nbCplx *R=new nbCplx[Npix3D];
    TF3DCplx(obj, A, dimVol);
     TF3DCplx(objDecal, B, dimVol);

///Calcul de R, correlation de phase dans Fourier=A*conj(B)/(modA.modB)
   for(cpt=0;cpt<Npix3D;cpt++)
        {
            double modA=sqrt( pow(A[cpt].Re,2)+pow(A[cpt].Im,2));
            double modB=sqrt( pow(B[cpt].Re,2)+pow(B[cpt].Im,2));
            if(modA==0)
                modA=1;
            if(modB==0)
                modB=1;
            R[cpt].Re= (A[cpt].Re*B[cpt].Re + A[cpt].Im*B[cpt].Im)/(modA*modB);
            R[cpt].Im=(-B[cpt].Im*A[cpt].Re + A[cpt].Im*B[cpt].Re)/(modA*modB);
        }

    delete[] A;

    ///Recalage en multipliant B Par R
   for(cpt=0;cpt<Npix3D;cpt++)
        {
            C[cpt].Re=B[cpt].Re*R[cpt].Re  -  R[cpt].Im*B[cpt].Im;
            C[cpt].Im=B[cpt].Im*R[cpt].Re  +  B[cpt].Re*R[cpt].Im;
        }
    TF3DCplx_INV(C, objRecal, dimVol);

    delete[] B;
    delete[] C;
    delete[] R;

}*/




/*Var2D corr_crois(double* obj2D_A, double* obj2D_B, Var2D dim)
{
    int NbPts2D=dim.x*dim.y;

//copie pour passer en complexe...
    nbCplx *copie_obj2D_A=new nbCplx[NbPts2D];
    nbCplx *copie_obj2D_B=new nbCplx[NbPts2D];
        for(int cpt=0;cpt<NbPts2D;cpt++)//copies pour passer en complexe (la fonction Tf2D demandent des cplx)
    {
        copie_obj2D_A[cpt].Re=obj2D_A[cpt];
        copie_obj2D_A[cpt].Im=0;
        copie_obj2D_B[cpt].Re=obj2D_B[cpt];
        copie_obj2D_B[cpt].Im=0;
    }
    //stockage des spectres
    nbCplx *spect2D_A= new nbCplx[NbPts2D];
    nbCplx *spect2D_B= new nbCplx[NbPts2D];
    //stocker la correlation et son spectre
    nbCplx *spectCorr=new nbCplx[NbPts2D];
    nbCplx *Corr=new nbCplx[NbPts2D];


    //SAV_Re(copie_obj2D_B,NbPts2D,"/home/mat/tomo_test/copie_objetB.bin",t_float,"wb");
    TF2Dcplx(copie_obj2D_A,spect2D_A,dim);
    TF2Dcplx(copie_obj2D_B,spect2D_B,dim);

//Elimination du module pour isoler  la phase : A*conj(B)/(modA*modB)

for(int cpt=0;cpt<NbPts2D;cpt++){
    //double modA=sqrt(pow(spect2D_A[cpt].Re,2)+pow(spect2D_A[cpt].Im,2));
    //double modB=sqrt(pow(spect2D_B[cpt].Re,2)+pow(spect2D_B[cpt].Im,2));
    spectCorr[cpt].Re=double(spect2D_A[cpt].Re*spect2D_B[cpt].Re+spect2D_A[cpt].Im*spect2D_B[cpt].Im);
    spectCorr[cpt].Im=double(spect2D_A[cpt].Im*spect2D_B[cpt].Re-spect2D_A[cpt].Re*spect2D_B[cpt].Im);
    }
    //TF de l'expoentielle contenant le dépahsage->translation dans l'espace direct
    TF2Dcplx_INV(spectCorr,Corr,dim);
    int cptMax=0;
    double valMax=0;
    double valModCorr=0;
    for(int cpt=0;cpt<NbPts2D;cpt++)
    {
     valModCorr=pow(Corr[cpt].Re,2)+pow(Corr[cpt].Im,2);

    if(valModCorr>valMax)
        {
        valMax=valModCorr;
            cptMax=cpt;
        }
    }
    Var2D decal2D_I={0,0};
    decal2D_I.x=cptMax%dim.x;
    decal2D_I.y=cptMax/dim.y;
    if(decal2D_I.x>15 && decal2D_I.x<dim.x-15)
    decal2D_I.x=0;
     if( decal2D_I.y>15 &&  decal2D_I.y<dim.y-15)
    decal2D_I.y=0;
    cout<<"---------------------------"<<endl;
    cout<<"X="<<decal2D_I.x<<endl;
    cout<<"Y="<<decal2D_I.y<<endl;
   // SAV_Re(spectCorr,NbPts2D,"/home/mat/tomo_test/SpectCorrRE.bin",t_float,"wb");
    //SAV_Re(Corr,NbPts2D,"/home/mat/tomo_test/Corr.bin",t_float,"wb");
    return decal2D_I;
}
*/
double bruit(int attenuation)
{
    double valBruit=0;
    valBruit=((double)rand()/(double)RAND_MAX);

    return valBruit;
}


void genere_rectang2D(double *objet,Var2D posI_Coin,Var2D dimRect,Var2D dim)
{
    for(int y=posI_Coin.y; y<posI_Coin.y+dimRect.y; y++)
    {
        int NbPixY=dim.x*y;
        for(int x=posI_Coin.x; x<posI_Coin.x+dimRect.x; x++)
        {
            int cpt3D=NbPixY+x;
            objet[cpt3D]=1.0;//+bruit(1);
        }
    }
}

void genere_rectang3D(nbCplx *objet,Var3D posI_Coin,Var3D dimRect,Var3D dim)
{
    for(int z=posI_Coin.z; z<posI_Coin.z+dimRect.z; z++)
    {
        int altitude=dim.x*dim.y*z;
        for(int y=posI_Coin.y; y<posI_Coin.y+dimRect.y; y++)
        {
            int NbPixY=dim.x*y;
            for(int x=posI_Coin.x; x<posI_Coin.x+dimRect.x; x++)
            {
                int cpt3D=altitude+NbPixY+x;
                objet[cpt3D].Re=1.0;//+bruit(1);
                objet[cpt3D].Im=7.0;//+bruit(1);
            }
        }
    }


}
void genere_OTF_T_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon)
{
    int k=0, dimPlanFinal=round(dim_final.x*dim_final.y), Kdx0=posSpec.x,Kdy0=posSpec.y, dimVolX=dim_final.x;
    ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
    double r2=rayon*rayon,Kdz0;
    double Kdz0_carre = rayon*rayon-Kdx0*Kdx0-Kdy0*Kdy0;
    if(round(Kdz0_carre)>-1)
    {
        Kdz0=sqrt(Kdz0_carre);
        int NXMAX_CARRE=NMAX.x*NMAX.x;
        for (int Kdy = -NMAX.y; Kdy < NMAX.y; Kdy++)   //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
        {
            int Kdy_carre=Kdy*Kdy;
            for (int Kdx = -NMAX.x; Kdx < NMAX.x; Kdx++)   //on balaye l'image 2D en y, centre au milieu
            {
                //int cpt=(Kdy+NMAX.y)*2*NMAX.x+Kdx+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D

                if(Kdx*Kdx+Kdy_carre<NXMAX_CARRE)
                {
                    //ne pas depasser l'ouverture numérique pour 1 hologramme
                    double Kdz_carre=r2-Kdx*Kdx-Kdy_carre; //altitude au carré des données
                    double Kdz=round(sqrt(Kdz_carre)-Kdz0);///-Kdz->réflexion
                    double altitude=(Kdz+decal.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur z!!
                    k=(-Kdx0+Kdx+decal.x)+(-Kdy0+Kdy+decal.y)*dimVolX+round(altitude);//indice du tableau 1D du volume 3D
                    OTFr[k].Re=1;//
                    OTFr[k].Im=1;//e
                }
            } //fin for y
        }
    }//fin if zm0>-1
}
void genere_OTF_RB_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon)
{
    int k=0, dimPlanFinal=round(dim_final.x*dim_final.y), Kdx0=posSpec.x,Kdy0=posSpec.y, dimVolX=dim_final.x;
    ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
    double r2=rayon*rayon, Kdz0;
    double Kdz0_carre = rayon*rayon-Kdx0*Kdx0-Kdy0*Kdy0;
    if(round(Kdz0_carre)>-1)
    {
        Kdz0=sqrt(Kdz0_carre);
        int NXMAX_CARRE=NMAX.x*NMAX.x;
        for (int Kdy = -NMAX.y; Kdy < NMAX.y; Kdy++)   //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
        {
            int Kdy_carre=Kdy*Kdy;
            for (int Kdx = -NMAX.x; Kdx < NMAX.x; Kdx++)   //on balaye l'image 2D en y, centre au milieu
            {
                //  int cpt=(Kdy+NMAX.y)*2*NMAX.x+Kdx+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D

                if(Kdx*Kdx+Kdy_carre<NXMAX_CARRE)
                {
                    //ne pas depasser l'ouverture numérique pour 1 hologramme
                    double Kdz_carre=r2-Kdx*Kdx-Kdy_carre; //altitude au carré des données
                    double Kdz=round(sqrt(Kdz_carre)+Kdz0);///
                    double altitude=(Kdz+decal.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur z!!
                    k=(+Kdx0-Kdx+decal.x)+(Kdy0-Kdy+decal.y)*dimVolX+round(altitude);//indice du tableau 1D du volume 3D
                    OTFr[k].Re=1;//
                    OTFr[k].Im=1;//e
                }
            } //fin for y
        }
    }//fin if zm0>-1
}
void genere_OTF_RH_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon)
{
    int k=0, dimPlanFinal=round(dim_final.x*dim_final.y), Kdx0=posSpec.x,Kdy0=posSpec.y, dimVolX=dim_final.x;
    ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
    double r2=rayon*rayon,Kdz0; //arg_z_arc=0,z_arc=0,Kdz0;
    double Kdz0_carre = rayon*rayon-Kdx0*Kdx0-Kdy0*Kdy0;
    if(round(Kdz0_carre)>-1)
    {
        Kdz0=sqrt(Kdz0_carre);
        int NXMAX_CARRE=NMAX.x*NMAX.x;
        for (int Kdy = -NMAX.y; Kdy < NMAX.y; Kdy++)   //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
        {
            int Kdy_carre=Kdy*Kdy;
            for (int Kdx = -NMAX.x; Kdx < NMAX.x; Kdx++)   //on balaye l'image 2D en y, centre au milieu
            {
               // int cpt=(Kdy+NMAX.y)*2*NMAX.x+Kdx+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D

                if(Kdx*Kdx+Kdy_carre<NXMAX_CARRE)
                {
                    //ne pas depasser l'ouverture numérique pour 1 hologramme
                    double Kdz_carre=r2-Kdx*Kdx-Kdy_carre; //altitude au carré des données
                    double Kdz=round(-sqrt(Kdz_carre)-Kdz0);///-Kdz->réflexion
                    double altitude=(Kdz+decal.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur z!!
                    k=(-Kdx0+Kdx+decal.x)+(-Kdy0+Kdy+decal.y)*dimVolX+round(altitude);//indice du tableau 1D du volume 3D
                    OTFr[k].Re=1;//
                    OTFr[k].Im=1;//e
                }
            } //fin for y
        }
    }//fin if zm0>-1
}

void multiplier_masqueCplx2(nbCplx *image, nbCplx *masque, int t_image, int t_mask, Var2D CentreI)
{

    int y=0;
    int x=0;
    int xMaskInf=round(t_mask/2);
    int yMaskInf=round(t_mask/2);
    //if((t_mask%2)!=0)
    //cout<<"fonction multiplier_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;
    if(CentreI.x<xMaskInf || CentreI.y<yMaskInf)
        cout<<"Débordement par le bord supérieur. ("<<CentreI.x<<","<<CentreI.y<<")"<<endl;
    if(t_image-CentreI.x>xMaskInf || t_image-CentreI.y<yMaskInf)
        cout<<"Débordement! Coordonnées supérieures du masque=("<<CentreI.x<<","<<CentreI.y<<")"<<endl;

    for(int xMask=-xMaskInf; xMask<xMaskInf; xMask++)
    {
        for(int yMask=-yMaskInf; yMask<yMaskInf; yMask++)
        {
            //si le masque déborde de l'image (attention le masque repasse du coté opposé
            ///vers coord info pour le masque
            int yiMask=yMask+yMaskInf;
            int xiMask=xMask+xMaskInf;
            ///compteur 1D des coordonnées
            y=yiMask+CentreI.y;
            x=xiMask+CentreI.x;
            int cptj=y*t_image+x;
            int cptMask=yiMask*t_mask+xiMask;
            cout<<cptMask<<endl;
            //masque[cptMask].Re=0;
            image[cptj].Re=image[cptj].Re*masque[cptMask].Re - image[cptj].Im*masque[cptMask].Im;
            image[cptj].Im=image[cptj].Re*masque[cptMask].Im + image[cptj].Im*masque[cptMask].Re;
        }
    }
}

void Plan_ds_VolCplx(nbCplx *Vol3D, nbCplx *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di)
{
    int x3Di=0, y3Di=0;
    // int taillePlan=dimPlan.x*dimPlan.y;
    int planVol3D=dimVol.x*dimVol.y;
    int ecartX=(dimVol.x-dimPlan.x)/2, ecartY=(dimVol.y-dimPlan.y)/2;
    int altitude_k=z3Di*planVol3D;
    int cpt2D=0, k=0;

    for (int yi = 0; yi < dimPlan.y; yi++)   //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
    {
        for (int xi = 0; xi < dimPlan.x; xi++)   //on balaye l'image 2D en y, centre au milieu
        {
            cpt2D=yi*dimPlan.x+xi;//calcul du cpt du tableau 1D de l'image 2D
            x3Di=xi+ecartX;
            y3Di=yi+ecartY;
            k=altitude_k+y3Di*dimVol.x+x3Di;
            Vol3D[k].Re=plan2D[cpt2D].Re;
            Vol3D[k].Im=plan2D[cpt2D].Im;
        }
    }
}
void Plan_ds_Vol(double *Vol3D, double *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di)
{
    int x3Di=0, y3Di=0;
    //int taillePlan=dimPlan.x*dimPlan.y;
    int planVol3D=dimVol.x*dimVol.y;
    int ecartX=(dimVol.x-dimPlan.x)/2, ecartY=(dimVol.y-dimPlan.y)/2;
    int altitude_k=z3Di*planVol3D;
    int cpt2D=0, k=0;

    for (int yi = 0; yi < dimPlan.y; yi++)   //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
    {
        for (int xi = 0; xi < dimPlan.x; xi++)   //on balaye l'image 2D en y, centre au milieu
        {
            cpt2D=yi*dimPlan.x+xi;//calcul du cpt du tableau 1D de l'image 2D
            x3Di=xi+ecartX;
            y3Di=yi+ecartY;
            k=altitude_k+y3Di*dimVol.x+x3Di;
            Vol3D[k]=plan2D[cpt2D];

        }
    }
}
void Vol_ds_Plan(double *Vol3D, double *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di)
{
    int x3Di=0, y3Di=0;
    // int taillePlan=dimPlan.x*dimPlan.y;
    int planVol3D=dimVol.x*dimVol.y;
    int ecartX=(dimVol.x-dimPlan.x)/2, ecartY=(dimVol.y-dimPlan.y)/2;
    int altitude_k=z3Di*planVol3D;
    int cpt2D=0, k=0;

    for (int yi = 0; yi < dimPlan.y; yi++)   //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
    {
        for (int xi = 0; xi < dimPlan.x; xi++)   //on balaye l'image 2D en y, centre au milieu
        {

            cpt2D=yi*dimPlan.x+xi;//calcul du cpt du tableau 1D de l'image 2D
            x3Di=xi+ecartX;
            y3Di=yi+ecartY;
            k=altitude_k+y3Di*dimVol.x+x3Di;
            plan2D[cpt2D]=Vol3D[k];
        }
    }
}


/*void  retroPropagSA(int deltaZ, nbCplx *fft_shift_norm, nbCplx * planObjet, Var3D decal, Var2D NMAX, double rayon)
 {
    double kz[2*NMAX.x*2*NMAX.y];
    nbCplx rephase[2*NMAX.x*2*NMAX.y];
    nbCplx spectrePropag[2*NMAX.x*2*NMAX.y];
    int dimPlan=round(2*NMAX.x*2*NMAX.y);

         ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
                                double r2=rayon*rayon;
                                double arg_z_arc=0,z_arc=0;
                                //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));
                                //double zm0_carre = rayon*rayon-xm0*xm0-ym0*ym0;

                                        int NXMAX_CARRE=NMAX.x*NMAX.x;
                                            ///CALCULER LES KZ-----------------------------------
                                        for (int y = -NMAX.y; y < NMAX.y; y++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
                                                int y_carre=y*y;
                                                for (int x = -NMAX.x; x < NMAX.x; x++) { //on balaye l'image 2D en y, centre au milieu
                                                        int cpt=(y+NMAX.y)*2*NMAX.x+x+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D
                                                        if(x*x+y_carre<NXMAX_CARRE)
                                                        { //ne pas depasser l'ouverture numérique pour 1 hologramme
                                                                double z_carre=r2-x*x-y_carre; //altitude au carré des données
                                                                int cpt=y*(NMAX.x)*2+x;
                                                                kz[cpt]=sqrt(z_carre);
                                                        }
                                                        else  {
                                                               int cpt=y*(NMAX.x)*2+x;
                                                                kz[cpt]=0;
                                                        }
                                                } //fin for y
                                        }
                                        ///---------------------CALCULER LE SPECTRE REPROPAGE-----------------------------------

                                                         for(int pix=0;pix<dimPlan;pix++)
                                                         {  int xc=0,yc=0;
                                                            rephase[pix].Re=cos(kz[pix]*deltaZ);
                                                            rephase[pix].Im=sin(kz[pix]*deltaZ);
                                                            spectrePropag[pix].Re=fft_shift_norm[pix].Re*rephase[pix].Re-fft_shift_norm[pix].Im*rephase[pix].Im;
                                                            spectrePropag[pix].Im=fft_shift_norm[pix].Re*rephase[pix].Im+fft_shift_norm[pix].Im*rephase[pix].Re;
                                                            //xc=pix%(2*NMAX.x)-NMAX.x, yc=pix/(2*NYMAX)-NMAX.y;
                                                         }
                                                          TF2Dcplx_INV(spectrePropag, planObjet, NMAX);

 }*/

void decalCoupeCplx(nbCplx *fft, nbCplx *fft_tmp, Var2D NMAX,Var2D dimCCD)
{
    for (int xi=0; xi<NMAX.x; xi++)
    {
        for (int yi=0; yi<NMAX.y; yi++)
        {
            int cpt1=yi*dimCCD.x+xi;
            int cpt2=yi*(NMAX.x)*2+xi;

            fft[cpt2].Re=fft_tmp[cpt1].Re;
            fft[cpt2].Im=fft_tmp[cpt1].Im;
        }
        for (int yi=dimCCD.y-NMAX.y; yi<dimCCD.y; yi++)
        {
            int cpt1=yi*dimCCD.x+xi;

            int cpt2 = xi+(yi-dimCCD.y+2*NMAX.y)*2*(NMAX.x);
            fft[cpt2].Re=fft_tmp[cpt1].Re;
            fft[cpt2].Im=fft_tmp[cpt1].Im;
        }
    }
    ///---////////////////////////////////////////////deuxieme demi-espace
    for(int xi=dimCCD.x-NMAX.x; xi<dimCCD.x; xi++)
    {
        for (int yi=0; yi<NMAX.y; yi++)
        {
            int cpt1=yi*dimCCD.x+xi;
            int cpt2=yi*(NMAX.x)*2+(xi-dimCCD.x+2*NMAX.x);

            fft[cpt2].Re=fft_tmp[cpt1].Re;
            fft[cpt2].Im=fft_tmp[cpt1].Im;
        }

        for (int yi=dimCCD.y-NMAX.y; yi<dimCCD.y; yi++)
        {
            int cpt1=yi*dimCCD.x+xi;

            int cpt2 = (xi-dimCCD.x+2*NMAX.x)+(yi-dimCCD.y+2*NMAX.y)*2*(NMAX.x);
            fft[cpt2].Re=fft_tmp[cpt1].Re;
            fft[cpt2].Im=fft_tmp[cpt1].Im;
        }
    }
}

void decalCoupe(double *fft_reel, double *fft_imag, double *fft_reel_tmp, double *fft_imag_tmp, Var2D NMAX,Var2D dimCCD)
{
    for (int xi=0; xi<NMAX.x; xi++)
    {
        for (int yi=0; yi<NMAX.y; yi++)
        {
            int cpt1=yi*dimCCD.x+xi;
            int cpt2=yi*(NMAX.x)*2+xi;

            fft_reel[cpt2]=fft_reel_tmp[cpt1];
            fft_imag[cpt2]=fft_imag_tmp[cpt1];
        }
        for (int yi=dimCCD.y-NMAX.y; yi<dimCCD.y; yi++)
        {
            int cpt1=yi*dimCCD.x+xi;

            int cpt2 = xi+(yi-dimCCD.y+2*NMAX.y)*2*(NMAX.x);
            fft_reel[cpt2]=fft_reel_tmp[cpt1];
            fft_imag[cpt2]=fft_imag_tmp[cpt1];
        }
    }
    ///---////////////////////////////////////////////deuxieme demi-espace
    for(int xi=dimCCD.x-NMAX.x; xi<dimCCD.x; xi++)
    {
        for (int yi=0; yi<NMAX.y; yi++)
        {
            int cpt1=yi*dimCCD.x+xi;
            int cpt2=yi*(NMAX.x)*2+(xi-dimCCD.x+2*NMAX.x);

            fft_reel[cpt2]=fft_reel_tmp[cpt1];
            fft_imag[cpt2]=fft_imag_tmp[cpt1];
        }

        for (int yi=dimCCD.y-NMAX.y; yi<dimCCD.y; yi++)
        {
            int cpt1=yi*dimCCD.x+xi;

            int cpt2 = (xi-dimCCD.x+2*NMAX.x)+(yi-dimCCD.y+2*NMAX.y)*2*(NMAX.x);
            fft_reel[cpt2]=fft_reel_tmp[cpt1];
            fft_imag[cpt2]=fft_imag_tmp[cpt1];
        }
    }
}
void InitTabCplx(nbCplx *z,int taille)//initailiser un tableau (taille totale="taille") de structure à zéro.
{
    for(int cpt=0; cpt<taille; cpt++)
    {
        z[cpt] = {0}; /* Tous les champs à zéro */
    }
}



/*int retroPropag_Born(nbCplx *TF3D_PotObj, nbCplx *TF_Uborn_norm, double * sup_redon, int dim_final, Var2D posSpec, Var3D decal3D, Var2D NMAX, double rayon)
{
                        int kxmi=posSpec.x, kymi=posSpec.y;
                        int kxm0=(kxmi-NMAX.x), kym0=(kymi-NMAX.x);//coordonnée dans l'image2D centrée (xm0,ym0)=(0,0)=au centre de l'image
                        int points_faux=0;
                        rayon=rayon;
                        int dimVolX=round(dim_final), dimPlanFinal=round(dim_final*dim_final);
                        float n0=1.515, lambda=pow(633,10^(-9)),kv=2*3.1416/lambda;

                                //création de variable pou-9r éviter N calculs dans la boucle sur le volume 3D
                                int cptPot=0; //indice tableau 1d des données du potentiel3D
                                double r2=rayon*rayon, kzm0, kzm0_carre = rayon*rayon-kxm0*kxm0-kym0*kym0;
                               // double r2=rayon*rayon, kzm0, kzm0_carre = n0*rayon*rayon-kxm0*kxm0-kym0*kym0;
                                //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));

                                if(round(kzm0_carre)>-1) {

                                        kzm0=sqrt(kzm0_carre);

                                        int NMAX_CARRE=NMAX.x*NMAX.x;

                                        float ctePotUb=1;//2/lambda;
                                        for (int fdy = -NMAX.y; fdy < NMAX.y; fdy++) { //on balaye le champ Uborn2D en x , origine (0,0) de l'image au milieu
                                                int fdy_carre=fdy*fdy;
                                                for (int kdx = -NMAX.x; kdx < NMAX.x; kdx++) { //on balaye l'image 2D en y, centre au milieu
                                                        int cpt=(fdy+NMAX.y)*2*NMAX.x+kdx+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D
                                                        if(kdx*kdx+fdy_carre<NMAX_CARRE) { //ne pas depasser l'ouverture numérique pour 1 hologramme
                                                                double kdz_carre=r2-kdx*kdx-fdy_carre; //altitude au carré des données
                                                                double koz=round(sqrt(kdz_carre)-kzm0);
                                                                 //double kz=n0*round(sqrt(kz_carre)-kzm0);
                                                                 double m=sqrt(rayon*rayon-kdx*kdx-fdy*fdy);
                                                                double altitude=(koz+decal3D.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur z!!

                                                                cptPot=(-kxm0+kdx+decal3D.x)+(-kym0+fdy+decal3D.y)*dimVolX+round(altitude);//indice du tableau 1D du volume 3D
                                                                //cout<<"k"<<k<<endl;
                                                                TF3D_PotObj[cptPot].Re+=-ctePotUb*m*TF_Uborn_norm[cpt].Im;//Inversion Re->Im à cause du coefficient i entre potentiel et Uborn
                                                                TF3D_PotObj[cptPot].Im+=ctePotUb*m*TF_Uborn_norm[cpt].Re;//
                                                                sup_redon[cptPot]+=1;//pour calculer le support
                                                        } else
                                                                points_faux++;
                                                } //fin for y
                                        }
                                }//fin if zm0>-1
}*/


void retroPropag_Born(vector <complex<double>> &TF3D_PotObj, vector<complex<double>> const &TF_Uborn_norm, vector<double>  &sup_redon, int dim_final, Var2D posSpec, Var3D decal3D, Var2D NMAX, double rayon, manip m1)
{ //int Nmax_obj=m1.NXMAX_OBJ;
//cout<<"Nmax_obj="<<Nmax_obj<<endl;
    int fxmi=posSpec.x, fymi=posSpec.y;
    //  cout<<"fi : "<<fxmi<<","<<fymi<<endl;
    int fxm0=(fxmi-NMAX.x), fym0=(fymi-NMAX.x);//coordonnée dans le repère humain (xm0,ym0)=(0,0)=au centre de l'image
    // cout<<"fm0 : "<<fxm0<<","<<fym0<<endl;
    int points_faux=0;

    int dimVolX=round(dim_final), dimPlanFinal=round(dim_final*dim_final);
    float n0=m1.n0, //refractive index of the background medium
    kv=2*PI/m1.lambda0, //wavevector in vaccuum
    k0=kv*n0; //wavevector in the background medium

    //création de variable pour éviter N calculs dans la boucle sur le volume 3D
    //create variable before loop
    int cptPot=0; //indice tableau 1d des données du potentiel3D
    double cteNorm=-2*PI;
    double r2=rayon*rayon, fzm0, fzm0_carre = rayon*rayon-fxm0*fxm0-fym0*fym0;
    double norm_altitude=1.0/rayon;//normaliser fdz pour passer en  sdz;
   // if(sqrt(fxm0*fxm0+fym0*fym0)<Nmax_obj){

    if(round(fzm0_carre)>=0)
    {
        fzm0=sqrt(fzm0_carre);
        int NMAX_CARRE=NMAX.x*NMAX.x;

        complex<double> cteUb2Pot(0,k0/PI);//

        //#pragma omp parallel for
        // cout<<"NMAX.y="<<NMAX.y<<endl;
        for(short int fdy = -NMAX.y; fdy < NMAX.y; fdy++)    //X cscan of the filed Uborn2D, origin (0,0) in the middle of the picture
        {

            //cout<<"-----------------------"<<fdy<<endl;
            //  cout<<"fdy haut="<<fdy<<endl;
            int fdy_carre=fdy*fdy;
            for (int fdx = -NMAX.x; fdx < NMAX.x; fdx++)    //on balaye l'image 2D en y, centre au milieu
            {
                int cpt=(fdy+NMAX.y)*2*NMAX.x+fdx+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D
                if(fdx*fdx+fdy_carre<m1.coef_NA_obj_limit*NMAX_CARRE)    //ne pas depasser l'ouverture numérique pour 1 hologramme
                {
                    // cout<<"---------------------"<<endl;
                   //  cout<<"fxm0,fym0="<<fxm0<<","<<fym0<<endl;
                    // cout<<"fdx,fdy="<<fdx<<","<<fdy<<endl;
                    double fdz_carre=r2-fdx*fdx-fdy_carre; //altitude au carré des données
                    double koz=round(sqrt(fdz_carre)-fzm0);
                    double sdz=sqrt(rayon*rayon-fdx*fdx-fdy*fdy)*norm_altitude;
                    double altitude=(koz+decal3D.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur koz!!
                   // cout<<"fxm0="<<fxm0<<endl;
                   // cout<<"dimPlanFinal="<<dimPlanFinal<<endl;
                    cptPot=(-fxm0+fdx+decal3D.x)+(-fym0+fdy+decal3D.y)*dimVolX+round(altitude);//indice du tableau 1D du volume 3D
                    TF3D_PotObj[cptPot]+=cteUb2Pot*cteNorm*sdz*TF_Uborn_norm[cpt];
                    sup_redon[cptPot]+=1;//redundacy in frequency support.
                }
                else points_faux++;
            }
        }
    }
   // }
}


double max(double *entree, int tailleTab)
{
    double valMax=0;
    for (int cpt=0; cpt<tailleTab; cpt++)
    {
        if(entree[cpt]>valMax)
            valMax=entree[cpt];
    }
    return valMax;
}
int coordMaxMod2D(nbCplx *entree, size_t const tailleTab)
{
    vector<double> moduleQ(tailleTab);
    for(size_t cpt=0; cpt<tailleTab; cpt++)
    {
        moduleQ[cpt]=pow(entree[cpt].Re,2)+pow(entree[cpt].Im,2);
    }
    int cptMax=0;
    double valMax=0;
    for (size_t cpt=0; cpt<tailleTab; cpt++)
    {
        if(moduleQ[cpt]>valMax)
        {
            cptMax=cpt;
            valMax=moduleQ[cpt];
        }
    }
    return cptMax;
}


void changeDim2D(double* tab, double* tabFinal, Var2D dimInit, Var2D dimFin)
{
    int diffX=round(dimFin.x-dimInit.x)/2;
    int diffY=round(dimFin.y-dimInit.y)/2;

    for(int x=0; x<dimInit.x; x++)
    {
        for( int y=0; y<dimInit.y; y++)
        {
            int cpt_init=y*dimInit.x+x;
            int cpt_final=(y+diffY)*dimFin.x+x+diffX;
            tabFinal[cpt_final]=tab[cpt_init];
        }
    }
}
void changeDim2DCplx(nbCplx *tab, nbCplx* tabFinal, Var2D dimInit, Var2D dimFin)
{
    int diffX=round(dimFin.x-dimInit.x)/2;
    int diffY=round(dimFin.y-dimInit.y)/2;

    for(int x=0; x<dimInit.x; x++)
    {
        for( int y=0; y<dimInit.y; y++)
        {
            int cpt_init=y*dimInit.x+x;
            int cpt_final=(y+diffY)*dimFin.x+x+diffX;
            tabFinal[cpt_final].Re=tab[cpt_init].Re;
            tabFinal[cpt_final].Im=tab[cpt_init].Im;
        }
    }
}



///Interpolation3D : attend un volume "troué" et les dimensions en x,y,et z
void interp_lin3D(vector <complex<double>> &volume_interp_3D)
{
    vector <complex<double>> zmin_max(volume_interp_3D.size());
    int dim=round(pow(volume_interp_3D.size(),1.0/3));
    cout<<"dim="<<dim<<endl;
    int dim_x=dim,dim_y=dim, dim_z=dim;
    int z=0;
    int z_min=0;//dim_z;
    int z_max=0;//zmin et zmax sont les 2 points définissant la droite d'interpolation

    for (int x=0; x < dim_x; x++){    //on balaye l'image, référentiel avec centre informatique
        for (int y=0; y<dim_y; y++){   //on balaye l'image,
            z=0;
            z_min=0,z_max=0;
            while(z<dim_z){  // pas de boucle for car l'indice est variable
                double valeur_test=volume_interp_3D[x+y*dim_x+z*dim_y*dim_x].real();
                if(volume_interp_3D[x+y*dim_x+z*dim_y*dim_x].real() !=0){   //si valeur !=0, alors borne inférieure
                    z_min=z; //on a trouvé z_min
                    z_max=z_min+1; //initialiser z_max au point suivant (sinon volume_interp_3D!=0 ->while jamais verifie)
                    ///find zmax=next point containing a data and !=dimZ
                    while( z_max<dim_z && volume_interp_3D[x+y*dim_x+z_max*dim_y*dim_x].real()==0){    //soit trouver prochain z_max, soit fin de colonne. Indifférent real() ou imag(), donc on prend arbitrairement real()
                        z_max=z_max+1;
                        z=z_max;
                    }
                    //cout<<"zmax="<<z_max<<endl;
                    ///zmin and zmax found, , if there is a "hole", we can interpolate : LINEAR INTERPOLATION
                    if(z_max<=dim_z && (z_max-z_min)>=2 && volume_interp_3D[x+y*dim_x+z_max*dim_y*dim_x].real() !=0){   //il faut au moins un trou pour interpoler
                       //y=ax+b-> interpolation linéaire
                        double a_real=(volume_interp_3D[x+y*dim_x+z_max*dim_y*dim_x].real()-volume_interp_3D[x+y*dim_x+z_min*dim_y*dim_x].real())/(z_max-z_min);
                        double b_real=volume_interp_3D[x+y*dim_x+z_max*dim_y*dim_x].real()-a_real*z_max;
                        double a_imag=(volume_interp_3D[x+y*dim_x+z_max*dim_y*dim_x].imag()-volume_interp_3D[x+y*dim_x+z_min*dim_y*dim_x].imag())/(z_max-z_min);
                        double b_imag=volume_interp_3D[x+y*dim_x+z_max*dim_y*dim_x].imag()-a_imag*z_max;

                        for(int cpt_z=z_min+1; cpt_z<z_max; cpt_z++){
                             volume_interp_3D[x+y*dim_x+cpt_z*dim_y*dim_x].real(a_real*cpt_z+b_real);
                             volume_interp_3D[x+y*dim_x+cpt_z*dim_y*dim_x].imag(a_imag*cpt_z+b_imag);
                        }
                    z=z_max-1; // nouveau compteur=ancienne borne sup (-1 car z++)
                    }
                }
                z++;
            }//fin z
        }//fin y
    } //fin x
}

//##############§FIN INTERP3D#####################################################


//#############################################################"
/*double *tukey2D(int dimx,int dimy, float alpha)
{
        int N=dimx;
        double *  tuk2D = new double[dimx*dimy];
        double tuk1Dx [dimx];
        double tuk1Dy [dimy];

        int borne1=round(alpha*(N-1)/2);
        int borne2=round((N-1)*(1-alpha/2));

        //memset(tuk2D, 0, dim_entree.x*dim_entree.y*8);
        //memset(tuk1Dx, 0, dim_entree.x*8);
        for(int cpt=0; cpt<borne1+1; cpt++)
                tuk1Dx[cpt]=0.5*(1+cos(3.1415*(2*cpt/(alpha*(N-1))-1)));
        for(int cpt=borne1+1; cpt<borne2+1; cpt++)
                tuk1Dx[cpt]=1;
        for(int cpt=borne2+1; cpt<N; cpt++)
                tuk1Dx[cpt]=0.5*(1+cos(3.1415*(2*cpt/(alpha*(N-1))-2/alpha+1)));


        for(int cpt=0; cpt<N*N; cpt++) {
                int cptx=cpt%(N);
                int cpty=cpt/(N);
                tuk2D[cpt]=tuk1Dx[cptx]*tuk1Dx[cpty];
        }
        return tuk2D;
}*/

void ecrire_rapport(int NXMAX,float rayon,float Rf,  int DIMX_CCD2,int coin_x, int coin_y,short int precision_exportation,string chemin,int nb_proj,float n1,float NA,float Tp, int G)
{
    time_t date;
    time(&date);
//    char nom_rapport[12];
    string nom_rapport2=chemin+"/rapport.txt";
//        concatener(chemin,"/rapport.txt",nom_rapport);
    FILE *fichier_rapport ;
    cout<<"nom_rapport:" <<nom_rapport2<<endl;
    /*  ouverture pour ecriture (w) en mode texte (t) */
    fichier_rapport = fopen (nom_rapport2.c_str(), "wt") ;
    if (fichier_rapport == NULL)
        printf ("impossible de créer le fichier rapport_calcul.txt\n");
    //fprintf(fichier_rapport,"Date     : %s\nNXMAX    : %i\n,Rayon : %i\nRf       : %f\nK : %f\ndimx_ccd : %d\ncoin_x   : %d\ncoin_y   : %d\nPrecision: %i bits\nSession  : %s\nnb_proj  : %i\nindice n1: %f\nNA       : %f\nT_pixel  : %e\nG        : %i",ctime(&date),NXMAX,rayon,Rf,K,DIMX_CCD2,coin_x,coin_y,8*sizeof(precision),chemin,nb_proj,n1,NA,Tp,G);
    fprintf(fichier_rapport,"Date : %s\n NXMAX=%i\n Rf=%f\n Precision: %i bits\n Session  : %s\n Nombre de projections : %i",ctime(&date),NXMAX,Rf,8*precision_exportation,chemin.c_str(),nb_proj);
    fclose (fichier_rapport);
}


void genereCache(double masque[], int t_image, int t_mask, int centreX, int centreY)
{
    int t_imageX=t_image;
    int t_imageY=t_image;
    if((t_mask%2)!=0)
        cout<<"fonction genere_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;
    for(int pixel=0; pixel<2*t_image*2*t_image; pixel++)
    {
        masque[pixel]=0;
    }
    for(int x=centreX-t_mask/2; x<centreX+t_mask/2; x++)
    {
        for(int y=centreY-t_mask/2; y<centreY+t_mask/2; y++)
        {
            //le masque repasse du coté gauche lorsque le cache touche le bord droit! A corriger (?)
            int x_2=x;
            int y_2=y;
            if(x>2*t_imageX)
                x_2=x-2*t_imageX;
            if(y>2*t_imageY)
                y_2=y-2*t_imageY;
            if(x<0)
                x_2=2*t_imageX+x;
            if(y<0)
                y_2=2*t_imageY+y;
            ///coordonnées 1D du centre
            int cptj=2*t_imageY*y_2+x_2;
            masque[cptj]=1;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////
//fonction chargeant une image 2D. Retourne un pointeur sur le tableau 1D qui la contient
//  Nécessite le nom du pointeur à remplir,
//le numéro du pas de phase shifting (1,2 3 ou 4), le chemin, et cpt_fichier
//pour boucler sur les 1000 projections
/////////////////////////////////////////////////////////////////////////////////////

/*
void charger_image2D(unsigned char* phasei, string imgFile, Var2D coin,Var2D taille)
{

        //rempli_tableau(phasei, imgFile, coin,taille);
        Image Monimage;
        //int i, currentImageWidth, currentImageHeight;
        Monimage.read(imgFile);////// chargement en memoire de l'image
        Monimage.type( GrayscaleType );	////// Mon image est N&B
        //Monimage.display();
        //Monimage.channelDepth();
        Monimage.compressType(NoCompression);
        //Monimage.crop(Geometry(dimx_ccd,dimy_ccd, 0, 0) );
        //Monimage.write("/home/mat/tomo_test/test.bmp");
        //currentImageWidth = Monimage.columns();////// extraction des caractéristiques de l'image
        //currentImageHeight = Monimage.rows();
        //	finalArray = new unsigned char[taille_x*taille_y];////// reservation de la taille memoire du tableau final
        ////// lecture de l'image
        Monimage.getPixels(coin.x,coin.y,taille.x,taille.y);
        ////// ecriture de l'image dans un tableau
        Monimage.writePixels(GrayQuantum, phasei);

}*/




/////////////////////////////////////////////////////////////////////////////////////
//fonction de circshift

/////////////////////////////////////////////////////////////////////////////////////

void circshift3(double* entree, double* result, Var2D dim,Var2D decal)
{
    //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
    decal.y=decal.y%dim.y;
    decal.x=decal.x%dim.x;

    for(int yi=0; yi<decal.y; yi++)
    {
        for(int xi=0; xi<decal.x; xi++)

        {
            int pixel=yi*dim.x+xi;
            int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
            //1er quadrant vers 4 eme
            result[pixel_shift]=entree[pixel];
            //4 eme quadrant vers 1er
            result[pixel]=entree[pixel_shift];
            //2eme vers 3eme
            result[(yi+decal.y)*dim.x+xi]=entree[pixel+decal.x];
            //3eme vers 2eme
            result[pixel+decal.x]=entree[(yi+decal.y)*dim.x+xi];
        }
    }
}
void decal2DGen(double* entree, double* result, Var2D dim,Var2D decal)
{
    //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
    decal.y=decal.y%dim.y;
    decal.x=decal.x%dim.x;
    if(decal.x<0)
        decal.x=dim.x+decal.x;
    if(decal.y<0)
        decal.y=dim.y+decal.y;
    for(int yi=0; yi<dim.y-decal.y; yi++)
    {
        for(int xi=0; xi<dim.x-decal.x; xi++)
        {
            //cout<<"xi,yi="<<xi<<","<<yi<<endl;
            int pixel=yi*dim.x+xi;
            int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
            result[pixel_shift]=entree[pixel];
        }
        for(int xi=dim.x-decal.x; xi<dim.x; xi++)
        {
            int pixel=yi*dim.x+xi;
            int pixel_shift=(yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
            result[pixel_shift]=entree[pixel];
        }
    }
    for(int yi=dim.y-decal.y; yi<dim.y; yi++)
    {
        for(int xi=0; xi<dim.x-decal.x; xi++)
        {
            int pixel=yi*dim.x+xi;
            int pixel_shift=(-dim.y+yi+decal.y)*dim.x+xi+decal.x;
            result[pixel_shift]=entree[pixel];
        }
        for(int xi=dim.x-decal.x; xi<dim.x; xi++)
        {
            int pixel=yi*dim.x+xi;
            int pixel_shift=(-dim.y+yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
            result[pixel_shift]=entree[pixel];
        }
    }

}
void circshift2DCplx(nbCplx* entree, nbCplx* result, Var2D dim,Var2D decal)///___/!\ ne fonctionne qu'avec des demi espace??
{
    //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo

    decal.y=decal.y%dim.y;
    decal.x=decal.x%dim.x;
    // cout<<"decal.x="<<decal.x<<"; decal.y="<<decal.y<<endl;
    for(int yi=0; yi<decal.y; yi++)
    {
        for(int xi=0; xi<decal.x; xi++)
        {
            int pixel=yi*dim.x+xi;
            //  cout<<"pixel="<<pixel<<endl;
            // cout<<"result[pixel]="<<result[pixel];
            int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
            //  cout<<"pixel_shift="<<pixel_shift<<endl;
            //1er quadrant vers 4 eme
            result[pixel_shift]=entree[pixel];

            //4 eme quadrant vers 1er
            result[pixel]=entree[pixel_shift];

            //2eme vers 3eme
            result[(yi+decal.y)*dim.x+xi]=entree[pixel+decal.x];

            //3eme vers 2eme
            result[pixel+decal.x]=entree[(yi+decal.y)*dim.x+xi];
        }
    }
}
void decal2DCplxGen(vector<complex<double>> const &entree, vector<complex<double>> &result, Var2D dim,Var2D decal)
{
    //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
    decal.y=decal.y%dim.y;
    decal.x=decal.x%dim.x;
    if(decal.x<0)
        decal.x=dim.x+decal.x;
    if(decal.y<0)
        decal.y=dim.y+decal.y;
    for(int yi=0; yi<dim.y-decal.y; yi++)
    {
        for(int xi=0; xi<dim.x-decal.x; xi++)
        {
            //cout<<"xi,yi="<<xi<<","<<yi<<endl;
            int pixel=yi*dim.x+xi;
            int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
            result[pixel_shift]=entree[pixel];
        }
        for(int xi=dim.x-decal.x; xi<dim.x; xi++)
        {
            int pixel=yi*dim.x+xi;
            int pixel_shift=(yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
            result[pixel_shift]=entree[pixel];
        }
    }
    for(int yi=dim.y-decal.y; yi<dim.y; yi++)
    {
        for(int xi=0; xi<dim.x-decal.x; xi++)
        {
            int pixel=yi*dim.x+xi;
            int pixel_shift=(-dim.y+yi+decal.y)*dim.x+xi+decal.x;
            result[pixel_shift]=entree[pixel];
        }
        for(int xi=dim.x-decal.x; xi<dim.x; xi++)
        {
            int pixel=yi*dim.x+xi;
            int pixel_shift=(-dim.y+yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
            result[pixel_shift]=entree[pixel];
        }
    }

}



void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D)
{
    decal3D.x=decal3D.x%dimFinal3D.x;//élmiiner les "modulos"
    decal3D.y=decal3D.y%dimFinal3D.y;
    decal3D.z=decal3D.z%dimFinal3D.y;

    unsigned short int xi,yi,zi=0;
    short int x2,y2,z2=0; //signé car une fois décalé, peuvent être négatifs!
    const unsigned int taille_plan =dimFinal3D.x*dimFinal3D.y;

    for(zi=0; zi<dimFinal3D.z; zi++)
    {
        if(zi+decal3D.z>dimFinal3D.z-1)   //dépassement à droite
        {
            z2=zi+decal3D.z-dimFinal3D.z;
        }
        else
        {
            if(zi+decal3D.z<0)   //dépassement à gauche
            {
                z2=dimFinal3D.z+(decal3D.z+zi);
            }
            else
            {
                z2=zi+decal3D.z;
            }
        }
        int nb_pixelz_decal=z2*taille_plan;
        unsigned int nb_pixelz=zi*taille_plan;
        for(yi=0; yi<dimFinal3D.y; yi++)
        {
            if(yi+decal3D.y>dimFinal3D.y-1)   //dépassement à droite
            {
                y2=yi+decal3D.y-dimFinal3D.y;
            }
            else
            {
                if(yi+decal3D.y<0)   //dépassement à gauche
                {
                    y2=dimFinal3D.y+(decal3D.y+yi);
                }
                else
                {
                    y2=yi+decal3D.y;
                }
            }
            int nb_lignes=yi*dimFinal3D.x;
            int nb_lignes_decal=y2*dimFinal3D.x;

            for(xi=0; xi<dimFinal3D.x; xi++)
            {
                if(xi+decal3D.x>dimFinal3D.x-1)   //dépassement à droite
                {
                    x2=xi+decal3D.x-dimFinal3D.x;
                }
                else
                {
                    if(xi+decal3D.x<0)   //dépassement à gauche
                    {
                        x2=dimFinal3D.x+(decal3D.x+xi);
                    }
                    else
                    {
                        x2=xi+decal3D.x;
                    }
                }

                volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2]=volume3D[nb_pixelz+nb_lignes+xi];
            }
        }
    }
}
void circshift3DCplx(vector <complex<double>> const &volume3D, vector <complex<double>> &volume3D_shift, Var3D dimFinal3D, Var3D decal3D)
{
    decal3D.x=decal3D.x%dimFinal3D.x;//éliminer les "modulos"
    decal3D.y=decal3D.y%dimFinal3D.y;
    decal3D.z=decal3D.z%dimFinal3D.y;

    unsigned short int xi,yi,zi=0;
    short int x2,y2,z2=0; //signé car une fois décalé, peuvent être négatifs!
    const unsigned int taille_plan =dimFinal3D.x*dimFinal3D.y;
    // #pragma omp parallel for private(zi)
    for(zi=0; zi<dimFinal3D.z; zi++)
    {
        if(zi+decal3D.z>dimFinal3D.z-1)   //dépassement à droite
        {
            z2=zi+decal3D.z-dimFinal3D.z;
        }
        else
        {
            if(zi+decal3D.z<0)   //dépassement à gauche
            {
                z2=dimFinal3D.z+(decal3D.z+zi);
            }
            else
            {
                z2=zi+decal3D.z;
            }
        }
        int nb_pixelz_decal=z2*taille_plan;
        unsigned int nb_pixelz=zi*taille_plan;
        for(yi=0; yi<dimFinal3D.y; yi++)
        {
            if(yi+decal3D.y>dimFinal3D.y-1)   //dépassement à droite
            {
                y2=yi+decal3D.y-dimFinal3D.y;
            }
            else
            {
                if(yi+decal3D.y<0)   //dépassement à gauche
                {
                    y2=dimFinal3D.y+(decal3D.y+yi);
                }
                else
                {
                    y2=yi+decal3D.y;
                }
            }
            int nb_lignes=yi*dimFinal3D.x;
            int nb_lignes_decal=y2*dimFinal3D.x;

            for(xi=0; xi<dimFinal3D.x; xi++)
            {
                if(xi+decal3D.x>dimFinal3D.x-1)   //dépassement à droite
                {
                    x2=xi+decal3D.x-dimFinal3D.x;
                }
                else
                {
                    if(xi+decal3D.x<0)   //dépassement à gauche
                    {
                        x2=dimFinal3D.x+(decal3D.x+xi);
                    }
                    else
                    {
                        x2=xi+decal3D.x;
                    }
                }
                volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2]=volume3D[nb_pixelz+nb_lignes+xi];
                // volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2].Re=volume3D[nb_pixelz+nb_lignes+xi].Re;
                // volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2].Im=volume3D[nb_pixelz+nb_lignes+xi].Im;
            }

        }
    }
}
///fonction de masquage
void multiplier_masque(double image[], unsigned char masque[], int t_image, int t_mask, int centreX, int centreY)
{
    int t_imageX=t_image;
    int t_imageY=t_image;
    //if((t_mask%2)!=0)
    //cout<<"fonction multiplier_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;

    for(int x=centreX-t_mask/2; x<centreX+t_mask/2; x++)
    {
        for(int y=centreY-t_mask/2; y<centreY+t_mask/2; y++)
        {
            //si le masque déborde de l'image (attention le masque repasse du coté gauche lorsque le jumeau touche le bord droit! A corriger
            int x_2=x;
            int y_2=y;
            if(x>2*t_imageX)
                x_2=x-2*t_imageX;
            if(y>2*t_imageY)
                y_2=y-2*t_imageY;
            if(x<0)
                x_2=2*t_imageX+x;
            if(y<0)
                y_2=2*t_imageY+y;
            ///coordonnées 1D du centre
            int cptj=2*t_imageY*y_2+x_2;
            //if((double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)!=1)
            ///cout<<"masque:"<<(double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)<<endl;
            image[cptj]=image[cptj]*masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)]/255;
        }
    }
}

///fonction de masquage
void multiplier_masque2(double image[], double masque[], int t_image, int t_mask, int centreX, int centreY)
{
    int t_imageX=t_image;
    int t_imageY=t_image;
    //if((t_mask%2)!=0)
    //cout<<"fonction multiplier_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;

    for(int x=centreX-t_mask/2; x<centreX+t_mask/2; x++)
    {
        for(int y=centreY-t_mask/2; y<centreY+t_mask/2; y++)
        {
            //si le masque déborde de l'image (attention le masque repasse du coté gauche lorsque le jumeau touche le bord droit! A corriger
            int x_2=x;
            int y_2=y;
            if(x>2*t_imageX)
                x_2=x-2*t_imageX;
            if(y>2*t_imageY)
                y_2=y-2*t_imageY;
            if(x<0)
                x_2=2*t_imageX+x;
            if(y<0)
                y_2=2*t_imageY+y;
            ///coordonnées 1D du centre
            int cptj=2*t_imageY*y_2+x_2;
            //if((double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)!=1)
            ///cout<<"masque:"<<(double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)<<endl;
            image[cptj]=image[cptj]*masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)];
        }
    }
}


/*
void recal_obj(nbCplx *a, nbCplx *b,nbCplx *objRecal, Var3D dimVol)
{
    int NPix3D=dimVol.x*dimVol.y*dimVol.z;
    nbCplx *A=new nbCplx[NPix3D];
    nbCplx *B=new nbCplx[NPix3D];



    TF3DCplx(a, A,dimVol);
    TF3DCplx(b, B,dimVol);


    nbCplx BConj[NPix3D];
    nbCplx prodAB[NPix3D];
    nbCplx decalPhi[NPix3D];
    conj_cplx(B, BConj, NPix3D);

    nbCplx produitABconj[NPix3D];

    AXB_cplx(A, BConj, prodAB, NPix3D);

for(int cpt=0;cpt<NPix3D;cpt++)
        {
            decalPhi[cpt].Re=prodAB[cpt].Re/sqrt(pow(A[cpt].Re,2))+sqrt(pow(BConj[cpt].Re,2));
            decalPhi[cpt].Im=prodAB[cpt].Im/sqrt(pow(A[cpt].Re,2))+sqrt(pow(BConj[cpt].Re,2));
        }

}
*/
void multiplier_masque2Cplx(nbCplx *image, double masque[], int t_image, int t_mask, Var2D Centre)
{
    int t_imageX=t_image;
    int t_imageY=t_image;
    //if((t_mask%2)!=0)
    //cout<<"fonction multiplier_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;

    for(int x=Centre.x-t_mask/2; x<Centre.x+t_mask/2; x++)
    {
        for(int y=Centre.y-t_mask/2; y<Centre.y+t_mask/2; y++)
        {
            //si le masque déborde de l'image (attention le masque repasse du coté gauche lorsque le jumeau touche le bord droit! A corriger
            int x_2=x;
            int y_2=y;
            if(x>2*t_imageX)
                x_2=x-2*t_imageX;
            if(y>2*t_imageY)
                y_2=y-2*t_imageY;
            if(x<0)
                x_2=2*t_imageX+x;
            if(y<0)
                y_2=2*t_imageY+y;
            ///coordonnées 1D du centre
            int cptj=2*t_imageY*y_2+x_2;
            //if((double(masque[t_mask*(y-Centre.y+t_mask/2)+(x-Centre.x+t_mask/2)])/255)!=1)
            ///cout<<"masque:"<<(double(masque[t_mask*(y-Centre.y+t_mask/2)+(x-Centre.x+t_mask/2)])/255)<<endl;
            image[cptj].Re=image[cptj].Re*masque[t_mask*(y-Centre.y+t_mask/2)+(x-Centre.x+t_mask/2)];
            image[cptj].Im=image[cptj].Im*masque[t_mask*(y-Centre.y+t_mask/2)+(x-Centre.x+t_mask/2)];
        }
    }
}

///découpe une une fenetre de dimension dim_dest, coin haut gauche coin, dans src de taille dim_src.
void coupeCplx(nbCplx *src, nbCplx *dest, Var2D dim_src, Var2D dim_dest, Var2D coin)
{
//    size_t nbPixSrc=dim_src.x*dim_src.y;
    size_t cpt_destX,cpt_destY, cpt_dest1D,
           cpt_srcX,cpt_srcY, cpt_src1D;

    for(cpt_destX=0; cpt_destX<dim_dest.x; cpt_destX++)
    {
        for(cpt_destY=0; cpt_destY<dim_dest.y; cpt_destY++)
        {

            cpt_dest1D=cpt_destX+cpt_destY*dim_dest.x;///coord 1D destination

            cpt_srcX=coin.x+cpt_destX;///coord X src
            cpt_srcY=coin.y+cpt_destY;///coord Y src
            cpt_src1D=cpt_srcX+cpt_srcY*dim_src.x;///coord 1D source

            dest[cpt_dest1D].Re=src[cpt_src1D].Re;
            dest[cpt_dest1D].Im=src[cpt_src1D].Im;

        }

    }

}

/*
void methodeCarre(int NbPixROI2d, double *holo1,  double *holo2,  double *holo3,  double *holo4)
{
///Calcul du déphasage réel (méthode de Carré), et du taux de modulation

    /// variable pour méthode de carré---------------------------
    double* holo1Re=new double[NbPixROI2d];
    double* holo2Re=new double[NbPixROI2d];
    double* holo3Re=new double[NbPixROI2d];
    double* holo4Re=new double[NbPixROI2d];
    double* txModulFrange=new double[NbPixROI2d];

    double* phaseMod2piIm=new double[NbPixROI2d];
    double* TfPhaseMod2piIm=new double[NbPixROI2d];
    double* TfPhaseMod2pi=new double[NbPixROI2d];

            memset(TfPhaseMod2piIm,0,sizeof(TfPhaseMod2piIm));

            double txModulFrange=0;

            double* holo1Im=new double[NbPixROI2d];
            double* holo2Im=new double[NbPixROI2d];
            double* holo3Im=new double[NbPixROI2d];
            double* holo4Im=new double[NbPixROI2d];
            //partie imaginaire nulle en entrée
            memset(holo1Im,0,sizeof(holo1Im));
            memset(holo2Im,0,sizeof(holo2Im));
            memset(holo3Im,0,sizeof(holo3Im));
            memset(holo4Im,0,sizeof(holo4Im));
    double denominAlpha=0;
    double numeratAlpha=0;
    double argTanAlpha=0;//\delta=2*alpha, bref vaudrait mieux trouver 45° et tan alpha autour de 1

    double* dephasageCarre=new double[NbPixROI2d];

    double *phaseMod2pi=new double[NbPixROI2d];//calcul de la phase [2pi] de l'hologramme, méthode de carré
    memset(phaseMod2pi,0,sizeof(phaseMod2pi));
/// fin déclaration variable carré----------------------------
    double somDephasageCarre=0;
    int nb_tan=0, cptPhiAberrant=0;

    //calcul de la valeur du saut de phase et de la phase elle même
    for(int pixel=0; pixel<NbPixROI2d; pixel++)
    {
        holo1Re[pixel]=(double)holo1[pixel];
        holo2Re[pixel]=(double)holo2[pixel];
        holo3Re[pixel]=(double)holo3[pixel];
        holo4Re[pixel]=(double)holo4[pixel];

        denominAlpha=((holo1Re[pixel]-holo4Re[pixel])+(holo2Re[pixel]-holo3Re[pixel]));
        if(denominAlpha!=0)
        {
            numeratAlpha=(3*(holo2Re[pixel]-holo3Re[pixel])-(holo1Re[pixel]-holo4Re[pixel]));
            argTanAlpha=numeratAlpha/denominAlpha;//formule de Carré
            if(argTanAlpha>0)
            {
                // dephasageCarre[pixel]=2*atan(sqrt(argTanAlpha))*180/3.14;
                dephasageCarre[pixel]=sqrt(argTanAlpha);
            }
            else
            {
                cptPhiAberrant++;
                dephasageCarre[pixel]=1;
            }
        }
        else
        {
            numeratAlpha=1;
            argTanAlpha=numeratAlpha/denominAlpha;//formule de Carré
            if(argTanAlpha>0)
            {
                // dephasageCarre[pixel]=2*atan(sqrt(argTanAlpha))*180/3.14;
                dephasageCarre[pixel]=sqrt(argTanAlpha);

            }
            else
            {
                cptPhiAberrant++;
                dephasageCarre[pixel]=1;
            }

        }//fin calcul Dphi
        ///calcul phase
        double denominPhase=holo2Re[pixel]+holo3Re[pixel]-holo1Re[pixel]-holo4Re[pixel];
        double argNumPhase=(3*(holo2Re[pixel]-holo3Re[pixel])-holo1Re[pixel]+holo4Re[pixel])*(holo1Re[pixel]-holo4Re[pixel]+holo2Re[pixel]-holo3Re[pixel]);
        if(argNumPhase>0)
        {
            if (denominPhase==0)
            {
                denominPhase=1;//solution la plus facile!
                phaseMod2pi[pixel]=atan(sqrt(argNumPhase)/denominPhase);
            }
            else
            {
                phaseMod2pi[pixel]=phaseMod2pi[pixel]=atan(sqrt(argNumPhase)/denominPhase);
                //  phaseMod2pi[pixel]=0;//sqrt((3*(holo2Re[pixel]-holo3Re[pixel])-holo1Re[pixel]+holo4Re[pixel])*(holo1Re[pixel]-holo4Re[pixel]+holo2Re[pixel]-holo3Re[pixel]))/denominPhase;
            }
        }
        else
        {
            phaseMod2pi[pixel]=0;
        }
        if(dephasageCarre[pixel]<1.4&&dephasageCarre[pixel]>0.6)//la valeur de la tangente doit valoir environ tan(0.5*pi/2)=1, on élimine les valeurs aberrantes>10
        {
            somDephasageCarre=dephasageCarre[pixel]+somDephasageCarre;
            txModulFrange[pixel]=2*sqrt((pow(holo4Re[pixel]-holo2Re[pixel],2)+pow(holo1Re[pixel]-holo3Re[pixel],2)))/(holo1Re[pixel]+holo2Re[pixel]+holo3Re[pixel]+holo4Re[pixel]);
            nb_tan++;
        }
    }
    //cout<<"Singularités de Delta phi : "<<cptPhiAberrant<<"%"<<endl;
    ///démodulation
    // phaseMod2pi=circshift(phaseMod2pi, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);
     holo1Re=circshift(holo1Re, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);
      holo2Re=circshift(holo2Re, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);
      holo3Re=circshift(holo3Re, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);
      holo4Re=circshift(holo4Re, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);*/
    //TF2D(phaseMod2pi,phaseMod2piIm,TfPhaseMod2pi,TfPhaseMod2piIm,DIMX_CCD2,DIMY_CCD2);
    //TfPhaseMod2pi=circshift(TfPhaseMod2pi, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);
///SAV tf phase
//SAV(TfPhaseMod2pi, NbPixROI2d, "/home/mat/tomo_test/IDP/TfPhaseMod2pi/TfPhaseMod2pi", t_char,"a+b");
///SAV delta
//SAV(dephasageCarre, NbPixROI2d, "/home/mat/tomo_test/IDP/deltaCarre/deltaCarre", t_char,"a+b");
///sav phase
//SAV(phaseMod2pi, NbPixROI2d, "/home/mat/tomo_test/IDP/phase/phase", t_char,"a+b");
///sav txModulation
//SAV(txModulFrange, NbPixROI2d, "/home/mat/tomo_test/IDP/txModulation/txModulation", t_char,"a+b");
///-------
    // double moyenne_tan_alpha=somme_tan_alpha/N;
    /*  for(int pixel=0;pixel<N;pixel++)
      		{
       			if(tan_alpha[pixel]<1.4&&tan_alpha[pixel]>0.6)//la valeur de la tangente doit valoir environ tan(0.5*pi/2)=1, on élimine les valeurs aberrantes>10
      			{
                 somme_ecart_carre=somme_ecart_carre+(tan_alpha[pixel]-moyenne_tan_alpha)*(tan_alpha[pixel]-moyenne_tan_alpha);//         calcul écart type de la valeur du saut de phase
      			}
      		}
     double ecart_type=sqrt(somme_ecart_carre/N);*/
    //affichage saut de phase
    // printf("la moyenne du phase shifting sur cette image vaut: %f avec nb_tan :%f et ecart type : %f\n",moyenne_tan_alpha,sqrt(nb_tan),ecart_type);
    //allocation résultat du phase shifting attention à ne pas désallouer hors de la boucle for
    //cout<<"Décalage phase estimée : "<<atan(moyenne_tan_alpha)*180/3.14*2<<"°"<<endl;
    //	cout<<"Taux moyen de modulation des Franges"<<txModulFrange<<endl;
    //Libérer la mémoire des 4 images expérimentales+methode carré
    //Libération mémoire carré.
   /* delete[] holo1Re,delete[] holo2Re,delete[] holo3Re,delete[] holo4Re;
    delete[] holo1, delete[] holo2, delete[] holo3, delete[] holo4;
    delete[] phaseMod2piIm;
    delete[] TfPhaseMod2piIm;
    delete[] TfPhaseMod2pi;
    delete[] txModulFrange;
    delete[] holo1Im, holo2Im, holo3Im, holo4Im;
    delete[] TfHolo1Re, TfHolo2Re, TfHolo3Re, TfHolo4Re;
     delete[] TfHolo1Im, TfHolo2Im, TfHolo3Im, TfHolo4Im;

    delete[] dephasageCarre;
    delete[] phaseMod2pi;
///Fin méthode de carré
}
*/
int CreerZoneFresnel(double *FresnelRe,double * FresnelIm, Var2D dim, Var2D centre,float d, float lambda)
{
    for(int cpty=0; cpty<dim.y; cpty++)
    {
        int hauteur=cpty*dim.x;
        for(int cptx=0; cptx<dim.x; cptx++)
        {
            int cpt=hauteur+cptx;
            FresnelRe[cpt]=cos(3.14159*(pow(cptx-centre.x,2)+pow(cpty-centre.y,2))/(lambda*d));
            FresnelIm[cpt]=sin(3.14159*(pow(cptx-centre.x,2)+pow(cpty-centre.y,2))/(lambda*d));
        }
    }
    return 1;
}
///#########lecture  d'un fichire binaire 3D, connaissant sa taille et son type de données
///read a  binary file. (path, 3D table, data format (double=64), Nbpixels to be read). for data format, can use enum PRECISION, cf "projet.h"
int get_bin_file_size(string chemin)
{
    size_t lTaille, nb_elmnt_lu;//size and number of elements
    //unsigned short int dimData=precision/8;//taille en octet d'un element.
    FILE* pFichier = NULL;
    pFichier = fopen(chemin.c_str(), "r");  //ouverture de ce fichier en écriture binaire

    if(pFichier==NULL){
        fputs("Impossible d'ouvrir le fichier\n",stderr);
        cout<<chemin<<endl;
        exit (1);// obtenir la longueur du fichier, comparer avec donnée entrée.
    }
    else{
        fseek(pFichier,0,SEEK_END);//trouver la fin de fichier
        lTaille = ftell (pFichier);//retourne la position courante (en octet) du curseur de fichier : ici, position de la fin du fichier
        rewind(pFichier);
        fclose(pFichier);
    }
return lTaille;
}

///#########lecture  d'un fichire binaire 3D, connaissant sa taille et son type de données
///read a  binary file. (path, 3D table, data format (double=64), Nbpixels to be read). for data format, can use enum PRECISION, cf "projet.h"
void lire_bin(string chemin, double resultat[], short int precision, const size_t NbPix)
{
    size_t lTaille, nb_elmnt_lu;//size and number of elements
    unsigned short int dimData=precision/8;//taille en octet d'un element.
    FILE* pFichier = NULL;
    pFichier = fopen(chemin.c_str(), "r");  //ouverture de ce fichier en écriture binaire

    if(pFichier==NULL){
        fputs("Impossible d'ouvrir le fichier\n",stderr);
        cout<<chemin<<endl;
        exit (1);// obtenir la longueur du fichier, comparer avec donnée entrée.
    }
    else{
        fseek(pFichier,0,SEEK_END);//trouver la fin de fichier
        lTaille = ftell (pFichier);//retourne la position courante (en octet) du curseur de fichier : ici, position de la fin du fichier->taille du fichier
        cout<<"Fichier "<<chemin<<endl;
         printf("taille trouvée en octet par ftell %li, taille estimée : %i\n",lTaille, NbPix*dimData);//
        rewind(pFichier);

        if(NbPix*dimData!=lTaille)
            cout<<"Taille du fichier "<<chemin <<" incompatible avec les dimensions\n"<<endl;

        nb_elmnt_lu = fread (resultat,1,lTaille,pFichier);//lecture
        //nb_elmnt_lu = fread (&resultat_vector[0], 1,lTaille,pFichier);

        if(nb_elmnt_lu!=lTaille){
            cout<<"Problème lors de la lecture du fichier "<<chemin<<endl;
            cout<<"Nombre d'éléments lus="<<nb_elmnt_lu<<endl;
        }
        fclose(pFichier);
    }

}

void SAV_Tiff2D(double *var_sav, string chemin, const size_t dim)
{
    TIFF *tif= TIFFOpen(chemin.c_str(), "a");
    TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, dim);
    TIFFSetField (tif, TIFFTAG_IMAGELENGTH, dim);
    TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField (tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField (tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

    tsize_t strip_size = TIFFStripSize (tif);
    tstrip_t strips_num = TIFFNumberOfStrips (tif);

    float* strip_buf=(float*)_TIFFmalloc(strip_size);
    for (unsigned int s=0; s<strips_num; s++)
    {
        for (unsigned int col=0; col<dim; col++)
        {
            unsigned int cpt=col+dim*s;
            strip_buf[col]=(float)var_sav[cpt];
        }
        TIFFWriteEncodedStrip (tif, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif);
    TIFFClose(tif);
}
/*
void SAVCplx(std::vector<complex<double> > var_sav, string partie, std::string chemin, enum PRECISION2 precision, char options[])
{
    //double* var_sav = &v[0];
    unsigned int cpt;
    unsigned int NbPix=var_sav.size();
    // cout<<"taille volume="<<NbPix<<endl;
    FILE *fichier_ID;
    fichier_ID= fopen(chemin.c_str(), options);
    if(fichier_ID==0)
        cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

    switch(precision)
    {
    case t_double:  //64 bit
    {
        double tampon=0;
        if(partie=="Re"|| partie=="re")
        {
            for(cpt=0; cpt<NbPix; cpt++)
            {
                tampon=var_sav[cpt].real();

                fwrite(&tampon,sizeof(tampon),1,fichier_ID);
            }
        }
        if(partie=="Im"|| partie=="im")
        {
            for(cpt=0; cpt<NbPix; cpt++)
            {
                tampon=var_sav[cpt].imag();
                fwrite(&tampon,sizeof(tampon),1,fichier_ID);
            }
        }
        break;
    }
    case t_float: //32 bits float
    {
        float tampon=0;

        if(partie=="Re"|| partie=="re")
        {

            for(cpt=0; cpt<NbPix; cpt++)
            {
                tampon=var_sav[cpt].real();
                //cout<<"tampon="<<tampon<<endl;
                fwrite(&tampon,sizeof(tampon),1,fichier_ID);
            }
        }
        if(partie=="Im"|| partie=="im")
        {
            for(cpt=0; cpt<NbPix; cpt++)
            {
                tampon=var_sav[cpt].imag();
                fwrite(&tampon,sizeof(tampon),1,fichier_ID);
            }
        }
        break;
    }
    }
    fclose(fichier_ID);
}
*/

void SAV_Tiff2D(std::vector<double> var_sav, string chemin, double taille_pixel)
{
    size_t dim=pow(var_sav.size(),0.5);
    float xres, yres;
    uint16  res_unit;
    TIFF *tif_out= TIFFOpen(chemin.c_str(), "w");//"w"->ecraser, "a"->ajouter
    TIFFSetField (tif_out, TIFFTAG_IMAGEWIDTH, dim);
    TIFFSetField (tif_out, TIFFTAG_IMAGELENGTH, dim);
    TIFFSetField (tif_out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tif_out, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif_out, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField (tif_out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField (tif_out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField (tif_out, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif_out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tif_out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    //fixer la résolution
    xres = yres = 0.01/taille_pixel; //nbpixel par resunit (par centimetre)=0.01m/taille_pixel en metre
    res_unit = RESUNIT_CENTIMETER;
    TIFFSetField(tif_out, TIFFTAG_XRESOLUTION, xres);
    TIFFSetField(tif_out, TIFFTAG_YRESOLUTION, yres);
    TIFFSetField(tif_out, TIFFTAG_RESOLUTIONUNIT, res_unit);

    tsize_t strip_size = TIFFStripSize (tif_out);
    tstrip_t strips_num = TIFFNumberOfStrips (tif_out);

    float* strip_buf=(float*)_TIFFmalloc(strip_size);
    for (unsigned int s=0; s<strips_num; s++)
    {
        for (unsigned int col=0; col<dim; col++)
        {
            unsigned int cpt=col+dim*s;
            strip_buf[col]=(float)var_sav[cpt];
        }
        TIFFWriteEncodedStrip (tif_out, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif_out);
    TIFFClose(tif_out);
}

void SAV_Tiff2DCplx(std::vector<complex<double>> var_sav, string partie, string chemin, double taille_pixel)
{
    const size_t dim=pow(var_sav.size(),0.5);
    float xres, yres;
    uint16  res_unit;
    TIFF *tif_out= TIFFOpen(chemin.c_str(), "w");//"w"->ecraser, "a"->ajouter
    TIFFSetField (tif_out, TIFFTAG_IMAGEWIDTH, dim);
    TIFFSetField (tif_out, TIFFTAG_IMAGELENGTH, dim);
    TIFFSetField (tif_out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tif_out, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif_out, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField (tif_out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField (tif_out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField (tif_out, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif_out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tif_out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    //fixer la résolution
    xres = yres = 0.01/taille_pixel; //nbpixel par resunit (par centimetre)=0.01m/taille_pixel en metre
    res_unit = RESUNIT_CENTIMETER;
    TIFFSetField(tif_out, TIFFTAG_XRESOLUTION, xres);
    TIFFSetField(tif_out, TIFFTAG_YRESOLUTION, yres);
    TIFFSetField(tif_out, TIFFTAG_RESOLUTIONUNIT, res_unit);

    tsize_t strip_size = TIFFStripSize (tif_out);
    tstrip_t strips_num = TIFFNumberOfStrips (tif_out);

    float* strip_buf=(float*)_TIFFmalloc(strip_size);
    for (size_t s=0; s<strips_num; s++)
    {
        for (unsigned int col=0; col<dim; col++)
        {
            unsigned int cpt=col+dim*s;
            if(partie=="Re" || partie=="re")
                strip_buf[col]=(float)var_sav[cpt].real();
            else
            {
                if(partie=="Im"|| partie=="im")
                    strip_buf[col]=(float)var_sav[cpt].imag();
                else
                    cout<<"La chaine fixant la partie doit être Re,re, Im ou im"<<endl;
            }
        }
        TIFFWriteEncodedStrip (tif_out, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif_out);
    TIFFClose(tif_out);
}

void SAV3D_Tiff(vector<double> var_sav, string chemin, double taille_pixel)
{
    const int dim=(int)round(std::pow(var_sav.size(), 1.0/3.0));
    uint32 image_width, image_height, dimz;
    float xres, yres;
    uint16 spp, bpp,  res_unit; //zpage;photo,
    TIFF *out;
    int x, y;
    //  float buffer2D[dim * dim];
    float *buffer2D=new float[dim*dim];
    out = TIFFOpen(chemin.c_str(), "w");
    if (!out)
        fprintf (stderr, "Can't open  for writing\n");

    image_width = dim;
    image_height = dim;
    dimz=dim;
    spp = 1; /* Samples per pixel */
    bpp = 32; /* Bits per sample */
    // photo = PHOTOMETRIC_MINISBLACK;

    for (int num_page = 0; num_page < dim; num_page++)//z=page
    {
        for (y = 0; y < dim; y++)
            for(x = 0; x < dim; x++)
            {
                buffer2D[y * dim + x] = (float)var_sav[y * dim + x+num_page*dim*dim];

            }


        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image_width / spp);
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, image_height);
        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
        TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        // TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
        TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
        TIFFSetField (out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP); //image en Floating point
        /* It is good to set resolutions too (but it is not nesessary) */
        xres = yres = 0.01/taille_pixel; //nbpixel par resunit (par centimetre, on multiplie par 0.01 pour tout passer en mètre)
        res_unit = RESUNIT_CENTIMETER;
        TIFFSetField(out, TIFFTAG_XRESOLUTION, xres);
        TIFFSetField(out, TIFFTAG_YRESOLUTION, yres);
        TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, res_unit);
        /* We are writing single page of the multipage file */
        TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
        /* Set the page number */
        TIFFSetField(out, TIFFTAG_PAGENUMBER, num_page, dimz);

        for (y = 0; y < image_height; y++)
            TIFFWriteScanline(out, &buffer2D[y * image_width], y, 0);

        TIFFWriteDirectory(out);
    }
    delete[] buffer2D;
    TIFFClose(out);
}
void SAV3D_Tiff(vector<complex <double>> var_sav, Var3D const dim, string partie, string chemin, double taille_pixel)
{
//    int dim=round(std::pow(var_sav.size(), 1.0/3.0));
    uint32 image_width, image_height;// dimz;
    float xres, yres;
    uint16 spp, bpp,  res_unit;//photo,zpage;
    TIFF *out;
    size_t x, y;// z;
    float *buffer2D=new float[dim.x * dim.y];
    out = TIFFOpen(chemin.c_str(), "w");
    if (!out)
        fprintf (stderr, "Can't open  for writing\n");
    image_width = dim.x;
    image_height = dim.y;
    //  dimz=dim.z;
    spp = 1; /* Samples per pixel */
    bpp = 32; /* Bits per sample */
    // photo = PHOTOMETRIC_MINISBLACK;
    size_t num_page=0;
    for(num_page = 0; num_page < dim.z; num_page++) //z=page
    {

        int nbPix_plan=num_page*dim.x*dim.y;
        ///reel
        if(partie=="Re" || partie=="re")
        {
            //  #pragma omp parallel for private(y)
            for (y = 0; y < dim.y; y++)
            {
                size_t num_lgn=y*dim.x;
                for(x = 0; x < dim.x; x++)
                {
                    buffer2D[num_lgn + x] = (float)var_sav[num_lgn + x+nbPix_plan].real();
                }
            }
        }
        ///imag
        else
        {
            if(partie=="Im" || partie=="im")
            {
                //  #pragma omp parallel for private(y)
                for (y = 0; y < dim.y; y++)
                {
                    size_t num_lgn=y*dim.x;
                    for(x = 0; x < dim.x; x++)
                    {
                        buffer2D[num_lgn + x] = (float)var_sav[num_lgn + x+nbPix_plan].imag();
                    }
                }
            }
            else
                cout<<"Partie non identifiée : Re || re, Im, || im"<<endl;
        }



//z=page
        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image_width / spp);
        //TIFFSetField(out, TIFFTAG_COMPRESSION, LZW_SUPPORT);
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, image_height);
        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
        TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        // TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
        TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
        TIFFSetField (out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP); //image en Floating point
        /* It is good to set resolutions too (but it is not nesessary) */
        xres = yres = 0.01/taille_pixel; //nbpixel par resunit (par centimetre, on multiplie par 0.01 pour tout passer en mètre)
        res_unit = RESUNIT_CENTIMETER;
        TIFFSetField(out, TIFFTAG_XRESOLUTION, xres);
        TIFFSetField(out, TIFFTAG_YRESOLUTION, yres);
        TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, res_unit);
        /* We are writing single page of the multipage file */
        TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
        /* Set the page number */
        TIFFSetField(out, TIFFTAG_PAGENUMBER, num_page, dim.z);
//auto start_tiff = std::chrono::system_clock::now();



        for (y = 0; y < image_height; y++) //écriture d'une page numérotée num_page, ligne par ligne (y).
        {
            TIFFWriteScanline(out, &buffer2D[y * image_width], y, 0);
        }
        // auto end_tiff = std::chrono::system_clock::now();
        //  auto elapsed_tiff = end_tiff - start_tiff;

        //  cout<<"numpage="<<num_page<<endl;

//   std::cout <<"Temps ecriture 1 plan tiff= "<< elapsed_tiff.count()/(pow(10,9)) << '\n';

        TIFFWriteDirectory(out);

    }

    delete[] buffer2D;
    TIFFClose(out);
}


void SAV3D_Tiff(vector<complex <double>> var_sav, string partie, string chemin, double taille_pixel)
{
    const size_t dim=round(std::pow(var_sav.size(), 1.0/3.0));
    uint32 image_width, image_height, dimz;
    float xres, yres;
    uint16 spp, bpp, res_unit;//photo, zpage;
    TIFF *out;
    size_t x, y;// z;

    //float *buffer2D=new float[dim * dim];
    std::vector<float> buffer2D(dim * dim);
    out = TIFFOpen(chemin.c_str(), "w");
    if (!out){
        fprintf (stderr, "Can't open  for writing\n");
        return;
    }
    image_width = dim;
    image_height = dim;
    dimz=dim;
    spp = 1; /* Samples per pixel */
    bpp = 32; /* Bits per sample */
    // photo = PHOTOMETRIC_MINISBLACK;
  //  size_t num_page=0;
    for(size_t num_page = 0; num_page < dim; num_page++) //z=page
    {
        int nbPix_plan=num_page*dim*dim;
        ///reel
        if(partie=="Re" || partie=="re")
        {
            //  #pragma omp parallel for private(y)
            #pragma omp parallel for
            for (y = 0; y < dim; y++)
            {
                size_t num_lgn=y*dim;
                for(x = 0; x < dim; x++)
                {
                    buffer2D[num_lgn + x] = (float)var_sav[num_lgn + x+nbPix_plan].real();
                }
            }
        }
        ///imag
        else
        {
            if(partie=="Im" || partie=="im")
            {
                //  #pragma omp parallel for private(y)
                for (y = 0; y < dim; y++)
                {
                    size_t num_lgn=y*dim;
                    for(x = 0; x < dim; x++)
                    {
                        buffer2D[num_lgn + x] = (float)var_sav[num_lgn + x+nbPix_plan].imag();
                    }
                }
            }
            else
                cout<<"Partie non identifiée : Re || re, Im, || im"<<endl;
        }
//z=page
        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image_width / spp);
        //TIFFSetField(out, TIFFTAG_COMPRESSION, LZW_SUPPORT);
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, image_height);
        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
        TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        // TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
        TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
        TIFFSetField (out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP); //image en Floating point
        /* It is good to set resolutions too (but it is not nesessary) */
        xres = yres = 0.01/taille_pixel; //nbpixel par resunit (par centimetre, on multiplie par 0.01 pour tout passer en mètre)
        res_unit = RESUNIT_CENTIMETER;
        TIFFSetField(out, TIFFTAG_XRESOLUTION, xres);
        TIFFSetField(out, TIFFTAG_YRESOLUTION, yres);
        TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, res_unit);
        /* We are writing single page of the multipage file */
        TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
        /* Set the page number */
        TIFFSetField(out, TIFFTAG_PAGENUMBER, num_page, dimz);
//auto start_tiff = std::chrono::system_clock::now();
        for (y = 0; y < image_height; y++) //écriture d'une page numérotée num_page, ligne par ligne (y).
        {
            TIFFWriteScanline(out, &buffer2D[y * image_width], y, 0);
        }
        // auto end_tiff = std::chrono::system_clock::now();
        //  auto elapsed_tiff = end_tiff - start_tiff;

        //  cout<<"numpage="<<num_page<<endl;

//   std::cout <<"Temps ecriture 1 plan tiff= "<< elapsed_tiff.count()/(pow(10,9)) << '\n';

        TIFFWriteDirectory(out);

    }

//    delete[] buffer2D;
    TIFFClose(out);
}
void SAV3D_Tiff_Optimized(const std::vector<std::complex<double>>& var_sav, const std::string& partie, const std::string& chemin, double taille_pixel)
{
    const size_t dim = std::round(std::cbrt(var_sav.size()));
    if (dim * dim * dim != var_sav.size()) {
        std::cerr << "Erreur : la taille du volume n'est pas un cube parfait." << std::endl;
        return;
    }

    const uint32 image_width = static_cast<uint32>(dim);
    const uint32 image_height = static_cast<uint32>(dim);
    const uint16 spp = 1;               // Samples per pixel
    const uint16 bpp = 32;              // Bits per sample
    const uint16 res_unit = RESUNIT_CENTIMETER;
    const float resolution = 0.01f / static_cast<float>(taille_pixel); // pixels/mètre -> pixels/cm

    std::vector<float> buffer2D(dim * dim);

    TIFF* out = TIFFOpen(chemin.c_str(), "w");
    if (!out) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier TIFF en écriture." << std::endl;
        return;
    }

    for (uint32 z = 0; z < dim; ++z)
    {
        const size_t offset = z * dim * dim;

        // Préparation du plan en parallèle
        #pragma omp parallel for
        for (size_t y = 0; y < dim; ++y) {
            for (size_t x = 0; x < dim; ++x) {
                size_t index = y * dim + x;
                const std::complex<double>& pixel = var_sav[offset + index];
                buffer2D[index] = (partie == "Re" || partie == "re") ? static_cast<float>(pixel.real())
                                 : (partie == "Im" || partie == "im") ? static_cast<float>(pixel.imag())
                                 : 0.0f;
            }
        }

        if (partie != "Re" && partie != "re" && partie != "Im" && partie != "im") {
            std::cerr << "Erreur : partie non reconnue (" << partie << "), utilisez 'Re' ou 'Im'." << std::endl;
            TIFFClose(out);
            return;
        }

        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image_width);
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, image_height);
        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
        TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
        TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
        TIFFSetField(out, TIFFTAG_XRESOLUTION, resolution);
        TIFFSetField(out, TIFFTAG_YRESOLUTION, resolution);
        TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, res_unit);
        TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
        TIFFSetField(out, TIFFTAG_PAGENUMBER, z, dim);
        //TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
       // TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
        for (uint32 y = 0; y < image_height; ++y) {
            if (TIFFWriteScanline(out, &buffer2D[y * image_width], y, 0) < 0) {
                std::cerr << "Erreur : écriture scanline échouée à z=" << z << ", y=" << y << std::endl;
                TIFFClose(out);
                return;
            }
        }

        TIFFWriteDirectory(out);
    }

    TIFFClose(out);
}

void Import3D_Tiff(vector<double> &imgTiff, string chemin, double taille_pixel)
{
    const size_t dim=round(std::pow(imgTiff.size(), 1.0/3.0));
    uint32 image_width, image_height, dimz;
//    float xres, yres;
    uint16 spp; //autres arguments : photo, res_unit, zpage, bpp
    TIFF *Tiff_id;
    size_t x, y; //z;
    float *buffer2D=new float[dim * dim];
    Tiff_id = TIFFOpen(chemin.c_str(), "r");
    if (!Tiff_id)
        fprintf (stderr, "Can't open  for writing\n");

    image_width = dim;
    image_height = dim;
    dimz=dim;
    spp = 1; /* Samples per pixel */
//    bpp = 32; /* Bits per sample */
    // photo = PHOTOMETRIC_MINISBLACK;
    size_t num_page=0;

    for(num_page = 0; num_page < dim; num_page++) //z=page
    {
//z=page
        TIFFGetField(Tiff_id, TIFFTAG_IMAGEWIDTH, image_width / spp);
        TIFFGetField(Tiff_id, TIFFTAG_IMAGELENGTH, image_height);
        /*  TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
          TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
          TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
         // TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
          TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
          TIFFSetField (out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP); //image en Floating point
          // It is good to set resolutions too (but it is not nesessary) */
        /* xres = yres = 0.01/taille_pixel; //nbpixel par resunit (par centimetre, on multiplie par 0.01 pour tout passer en mètre)
         res_unit = RESUNIT_CENTIMETER;
         TIFFSetField(out, TIFFTAG_XRESOLUTION, xres);
         TIFFSetField(out, TIFFTAG_YRESOLUTION, yres);
         TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, res_unit);
         // We are writing single page of the multipage file
        // TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
         // Set the page number
        //  TIFFSetField(out, TIFFTAG_PAGENUMBER, num_page, dimz);*/

        for (y = 0; y < image_height; y++) //écriture d'une page numérotée num_page, ligne par ligne (y).
        {
            TIFFReadScanline(Tiff_id, &buffer2D[y * image_width], y, 0);
        }
        int nbPix_plan=num_page*dim*dim;
        for (y = 0; y < dim; y++)
        {
            size_t num_lgn=y*dim;
            for(x = 0; x < dim; x++)
            {
                imgTiff[num_lgn + x+nbPix_plan]=buffer2D[num_lgn + x];
            }
        }
    }
    delete[] buffer2D;
    TIFFClose(Tiff_id);
}

void SAV_Tiff3D(nbCplx *var_sav, string chemin_ind, string chemin_abs, int dim)
{
    double *bufferRe=new double[dim*dim];
    double *bufferIm=new double[dim*dim];
    for(int z=0; z < dim; z++)
    {
        for(int y=0; y < dim; y++)
        {
            for(int x=0; x < dim; x++)
            {
                int cpt1= x + y * dim + z * dim * dim;
                int cpt2= x + y * dim;
                bufferRe[cpt2]=var_sav[cpt1].Re;
                bufferIm[cpt2]=var_sav[cpt1].Im;
            }
        }
        SAV_Tiff2D(bufferRe, chemin_ind, dim);
        SAV_Tiff2D(bufferIm, chemin_abs, dim);
    }
    delete[] bufferRe;
    delete[] bufferIm;
}

