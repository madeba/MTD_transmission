#include "scan_functions.h"
#include "math.h"
#include <vector>
#include <numeric>
#include <iostream>
using namespace std;

float2D maj_fleur(float rho, int nbHolo, float *theta)//tension x,y/rayon/nbHolo,/angle
{
    float2D Vout={0,0};
    float t=*theta;//init angle polaire pour la courbe paramétree
    // float LongTot=17.15*rho, delta_curv=LongTot/nbHolo;//longueur totale, delta_abs_curv "régulier"
    //float delta_theta=(1/rho*delta_curv)/(sqrt(15*sin(4*t)*sin(4*t)+1));
    float delta_theta=6.28/nbHolo;
    t=delta_theta+t; ///abscisse curv  point suivant
    Vout.x=rho*cos(4*t)*cos(t);
    Vout.y=rho*cos(4*t)*sin(t);
    *theta=t;///maj theta prog principal
    return Vout;
}

void scan_fleur(float rho, int nbHolo, std::vector<float2D> &Vout_table)//tension x,y/rayon/nbHolo,/angle
{
    double theta=0;//init angle polaire pour la courbe paramétree
    // float LongTot=17.15*rho, delta_curv=LongTot/nbHolo;//longueur totale, delta_abs_curv "régulier"
    //float delta_theta=(1/rho*delta_curv)/(sqrt(15*sin(4*t)*sin(4*t)+1));
    double delta_theta=2*M_PI/nbHolo;
    for(int numHolo=0;numHolo<nbHolo;numHolo++){
        Vout_table[numHolo].x=rho*cos(4*theta)*cos(theta);
        Vout_table[numHolo].y=rho*cos(4*theta)*sin(theta);
        theta=delta_theta+theta; ///abscisse curv  point suivant
    }
}

//concentric circles, export pattern as a table
void scan_annular(float Vmax, float nb_circle, int nbHolo, std::vector<float2D> &Vout_table)//tension x,y/rayon/nbHolo,/angle
{
    double delta_rho=(Vmax/nb_circle);//delta_frequency (voltage) between successive radius
    cout<<"delta_rho="<<delta_rho<<endl;
    double total_length=0;
    double perimetre=0;
    int numHolo=0;
    vector<double> kiz(nbHolo);
    vector<size_t> nbPointCircle(nb_circle);//nb of point in each circle.
    vector<double> rho_tension(nb_circle);//table with the different radius values
    size_t numCircle=0;
    ///créer le tableau des 4 tensions
    for(double rho=Vmax; rho>0; rho=rho-delta_rho){
            rho_tension[numCircle]=rho;
            cout<<"rho_tension="<<rho<<endl;
            numCircle++;
        }
    total_length=2*M_PI*accumulate(rho_tension.begin(),rho_tension.end(),0);

    numCircle=0;
    size_t nbPtRestant=nbHolo;
    //répartir le nombre de point dans le tableau nbPointcircle
    //sans avoir de probleme lié à l'échantillonnage
    for(int numCircle=nb_circle-1; numCircle>=0; numCircle--)
        {
            perimetre=2*M_PI*rho_tension[numCircle];
            if (numCircle!=0)
                {
                    nbPointCircle[numCircle]=round((perimetre/total_length)*nbHolo);
                    nbPtRestant=nbPtRestant-nbPointCircle[numCircle];
                }
            else nbPointCircle[numCircle]=nbPtRestant;
        }
    numCircle=0;
    // for(int numCircle=0;numCircle<nb_circle; numCircle++)
    //balayage en num_circle et nbPoint pour éviterr les pb d'arrondi
    for(int numCircle=0; numCircle<nb_circle; numCircle++)
        {
            int nbPoint=nbPointCircle[numCircle];
            double delta_theta=2*M_PI/nbPoint;
            cout<<"rho_tension="<<rho_tension[numCircle]<<endl;

            for(int cpt=0; cpt<nbPoint; cpt++)
                {
                    Vout_table[numHolo].x=rho_tension[numCircle]*cos(delta_theta*cpt);
                    Vout_table[numHolo].y=rho_tension[numCircle]*sin(delta_theta*cpt);
                    if(pow(Vmax,2)-pow(Vout_table[numHolo].x,2)-pow(Vout_table[numHolo].y,2)>0)
                    kiz[numHolo]=sqrt(pow(Vmax,2)-pow(Vout_table[numHolo].x,2)-pow(Vout_table[numHolo].y,2));
                    else{
                        kiz[numHolo]=0;
                    }
                    //cout<<pow(Vmax,2)-pow(Vout_table[numHolo].x,2)-pow(Vout_table[numHolo].y,2)<<endl;
                    cout<<kiz[numHolo]<<endl;
                    numHolo++;//next scan angle
                }
        }
}
//spirale Fermat avec angle d'or
void scan_fermat(double Vmax, int nbHolo, std::vector<float2D> &Vout_table, Var2D dim)
{
    double theta=0,theta_max=nbHolo;//number of holograms = theta max in radian !!
    double angle_dor=(3*sqrt(5))*M_PI;
    double delta_theta=theta_max/nbHolo;
    cout<<"delta_theta="<<delta_theta<<endl;
    vector<double> centre(dim.x*dim.y);
    double nbPixVolt=(dim.x/(2*Vmax));
    Var2D spec={0,0},specI={0,0};

    for(int numHolo=0;numHolo<nbHolo;numHolo++){
    Vout_table[numHolo].x=sqrt((theta)/theta_max)*Vmax*cos((theta)*angle_dor);
    Vout_table[numHolo].y=sqrt((theta)/theta_max)*Vmax*sin((theta)*angle_dor);

     spec.x=round(Vout_table[numHolo].x*nbPixVolt);
    spec.y=round(Vout_table[numHolo].y*nbPixVolt);
    specI.x=spec.x+dim.x/2;
    specI.y=spec.y+dim.y/2;
    centre[specI.y*dim.x+specI.x]=numHolo;

    theta=theta+delta_theta;
    }
    //redistribuer les points dans un tableau de taille nbHolo
    //en limitant les mouvements brusques du miroir
    int numHolo=0;
    int x=0,y=0;
    for(int cpt=0; cpt<dim.x*dim.y; cpt++){///2) balayage du pattern pour trier à y croissant (et x en zig zag)
        y=cpt/dim.x;
        if(y%2==0) x=dim.x-cpt%dim.x;//si y paire, balayer à x décroissant
        else x=cpt%dim.x;//si y impaire, balayer à x croissant
        ///y'a t-il un point dans le pattern ?
        if(centre[y*dim.x+x]!=0){
            Vout_table[numHolo].x=x/nbPixVolt-Vmax;//tension recalculée à partir de la position du point
            Vout_table[numHolo].y=y/nbPixVolt-Vmax;
            numHolo++;
        }
    }
}
///generate a random cartesian coordinate between [-vmax, +vmax], return the values in Vout_table
void scan_random_cartesian(double Vmax, int nbHolo, std::vector<float2D> &Vout_table, Var2D dim)
{
    double const nbPixVolt=(dim.x/(2*Vmax));
    srand (time(NULL));/* initialize random seed: */
    double nbAleaX=0, nbAleaY=0;
    vector<double> centre(dim.x*dim.y);
    int numHolo=0;
    while(numHolo<nbHolo){
        nbAleaX=Vmax*(2*rand()/(double)RAND_MAX)-1;
        nbAleaY=Vmax*(2*rand()/(double)RAND_MAX)-1;
        if(nbAleaX*nbAleaX+nbAleaY*nbAleaY<Vmax*Vmax){
            Vout_table[numHolo].x=nbAleaX;
            Vout_table[numHolo].y=nbAleaY;
            numHolo++;
        }
    }
    Var2D spec={0,0},specI={0,0};
    //placer les points alatoire sur la grille Vx, Vy
    for(int numHolo=0; numHolo<nbHolo; numHolo++){
        spec.x=round(Vout_table[numHolo].x*nbPixVolt);
        spec.y=round(Vout_table[numHolo].y*nbPixVolt);
        specI.x=spec.x+dim.x/2;
        specI.y=spec.y+dim.y/2;
        centre[specI.y*dim.x+specI.x]=1;
    }
    //redistribuer les points dans un tableau de taille nbHolo
    //en limitant les mouvements brusques du miroir
     numHolo=0;
    //balayer les lignes paires à x croissant//scan
    for( int y=0; y<dim.y; y=y+2){
        for(int x=0; x<dim.x; x++){
            if(centre[y*dim.x+x]==1){
                Vout_table[numHolo].x=x/nbPixVolt-Vmax;
                Vout_table[numHolo].y=y/nbPixVolt-Vmax;
                numHolo++;
            }
        }
    }
    //puis revenir en balayant à x décroissant les lignes impaires
    for( int y=1; y<dim.y; y=y+2){
        for(int x=dim.x; x<dim.x; x--){
            if(centre[y*dim.x+x]==1){
                Vout_table[numHolo].x=x/nbPixVolt-Vmax;
                Vout_table[numHolo].y=y/nbPixVolt-Vmax;
                numHolo++;
            }
        }
    }
}

//generate random voltages in polar coordinates, which modify the points distribution (cos(random)sin(random))
void scan_random_polar3D(double Vmax, int nbHolo, std::vector<float2D> &Vout_table, Var2D dim)
{
    double const nbPixVolt=(dim.x/(2*Vmax));
    srand (time(NULL));/* initialize random seed: */
    double nbAleaPhi=0, nbAleaTheta=0;
    vector<int> centre(dim.x*dim.y);
    int numHolo=0;
    while(numHolo<nbHolo)//generate random (float) voltages [-Vmax;+Vmax], and stock it in a integer table "specular-like" [-dimHolo,+dimHolo]
    {
        nbAleaPhi=acos(2.0*(rand()/(double)RAND_MAX)-1.0);//Phi=[0,pi]
        nbAleaTheta=rand()/(double)RAND_MAX*2*M_PI;//theta=[0;2pi]
        Vout_table[numHolo].x = Vmax*cos(nbAleaTheta)*sin(nbAleaPhi);
        Vout_table[numHolo].y= Vmax*sin(nbAleaTheta)*sin(nbAleaPhi);

        int x=round((Vout_table[numHolo].x+Vmax)*nbPixVolt);
        int y=round((Vout_table[numHolo].y+Vmax)*nbPixVolt);
        if(centre[y*dim.x+x]!=1)//test if the random point is unique
        {
        centre[y*dim.x+x]=1;
        numHolo++;
        }
    }
    //réordonner les points dans un tableau de taille nbHolo
    //afin de limiter les mouvements brusques du miroir
    numHolo=0;
    for( int y=0; y<dim.y; y=y+1){
        if(y%2==0)//y pair : balayage à x croissant
        {
            for( int x=0; x<dim.x; x++){
                if(centre[y*dim.x+x]==1){
                    Vout_table[numHolo].x=x/nbPixVolt-Vmax;
                    Vout_table[numHolo].y=y/nbPixVolt-Vmax;
                    numHolo++;
                 }
            }
        }
        else //y impair: balayage à x décroissant
        {
            for (int x=dim.x-1;x>=0;x--){
                if(centre[y*dim.x+x]==1){
                Vout_table[numHolo].x=x/nbPixVolt-Vmax;
                Vout_table[numHolo].y=y/nbPixVolt-Vmax;
                numHolo++;
                }
            }
        }
    }
}

float2D maj_spiral(float rho, int nbHolo, float *theta)//tension x,y/rayon/nbHolo,/angle
{
    float2D Vout={0,0};
    double numbTurn=10;
    double thetaMax=2*numbTurn*M_PI,thetaMin=0;
    double t=*theta;//init angle polaire pour la courbe paramétree
    //longueur totale, delta_abs_curv "régulier"
    double longTot=floor(0.5*(log(thetaMax+sqrt(thetaMax*thetaMax+1))+thetaMax*sqrt(thetaMax*thetaMax+1)-
                         log(thetaMin+sqrt(thetaMin*thetaMin+1))+thetaMin*sqrt(thetaMin*thetaMin+1)));
    double delta_theta=longTot/nbHolo*1/(sqrt(1+t*t));
    t=delta_theta+t; ///abscisse curv  point suivant
    Vout.x = (rho*t/thetaMax)*sin(t);
    Vout.y =(rho*t/thetaMax)*cos(t);
    *theta=t;///maj theta prog principal
    return Vout;
}
void scan_spiralOS(double Vmax, int nbHolo, vector<float2D> &Vout_table, unsigned short int nbTurn)//tension x,y/rayon/nbHolo,/angle
{
    double const numbTurn=nbTurn;
    int const nbSample=10000;
    vector<float2D> Vout_tableOS(nbSample);
    double thetaMax=2*numbTurn*M_PI,thetaMin=0, theta=0;
    //longueur totale, delta_abs_curv "régulier"
    double longTot=floor(0.5*(log(thetaMax+sqrt(thetaMax*thetaMax+1))+thetaMax*sqrt(thetaMax*thetaMax+1)-
                         log(thetaMin+sqrt(thetaMin*thetaMin+1))+thetaMin*sqrt(thetaMin*thetaMin+1)));

    for(int numHolo=0; numHolo<nbSample; numHolo++){
        double delta_theta=longTot/nbSample*1/(sqrt(1+theta*theta));
        theta=delta_theta+theta; ///abscisse curv  point suivant
        Vout_tableOS[numHolo].x = (Vmax*theta/thetaMax)*sin(theta);
        Vout_tableOS[numHolo].y =(Vmax*theta/thetaMax)*cos(theta);
    }
    int ratioOS=nbSample/nbHolo; ///oversampling ratio
    for(int numHolo=0; numHolo<nbHolo; numHolo++){
        Vout_table[numHolo].x=Vout_tableOS[numHolo*ratioOS].x;
        Vout_table[numHolo].y=Vout_tableOS[numHolo*ratioOS].y;
    }
}

void scan_spiral(double Vmax, int nbHolo, std::vector<float2D> &Vout_table, Var2D dim, unsigned short int nbTurn)//tension x,y/rayon/nbHolo,/angle
{
    double numbTurn=nbTurn, thetaMax=2*numbTurn*M_PI,thetaMin=0, theta=0;
    //longueur totale, delta_abs_curv "régulier"
    double longTot=floor(0.5*(log(thetaMax+sqrt(thetaMax*thetaMax+1))+thetaMax*sqrt(thetaMax*thetaMax+1)-
                         log(thetaMin+sqrt(thetaMin*thetaMin+1))+thetaMin*sqrt(thetaMin*thetaMin+1)));

    for(int numHolo=0; numHolo<nbHolo; numHolo++){
        double delta_theta=longTot/nbHolo*1/(sqrt(1+theta*theta));
        theta=delta_theta+theta; ///abscisse curv  point suivant
        Vout_table[numHolo].x = (Vmax*theta/thetaMax)*sin(theta);
        Vout_table[numHolo].y =(Vmax*theta/thetaMax)*cos(theta);
        //cout<< Vout_table[numHolo].x<<endl;
    }
}

void scan_uniform3D(double Vmax, int nbHolo, std::vector<float2D> &Vout_table)
{
    double surface_moyenne_pt=2*M_PI/nbHolo;
    double distance_moyenne=sqrt(surface_moyenne_pt);
    vector<double> kiz(nbHolo);
    cout<<"distance moyenne="<<distance_moyenne<<endl;

    double theta_max=M_PI/2,theta_min=0;
    int nbCercle=round((theta_max-theta_min)/distance_moyenne);//MTHETA=nbCercle concentrique
    cout<<"nb_total_cercle="<<nbCercle<<endl;
    double delta_theta=(theta_max-theta_min)/nbCercle;
    double theta=0,phi=0,delta_curv_phi=0;
    int numHolo=0;
    double total_circle_length=0;
    for(int num_cercle=0;num_cercle<nbCercle;num_cercle++){//scan polar half diameter
        theta=(num_cercle)*delta_theta;//theta=0, singularité donc on exclut num_cercle==0
        delta_curv_phi=surface_moyenne_pt/delta_theta; //donc d ?
        cout<<"distance delta_curv_phi"<<delta_curv_phi<<endl;
        total_circle_length=total_circle_length+2*M_PI*sin(theta);
    }
    cout<<"total length="<<total_circle_length<<endl;
    double dist_recalc=total_circle_length/nbHolo;
    cout<<"distance moyenne recalculée="<<total_circle_length/nbHolo<<endl;

    for(int num_cercle=0;num_cercle<nbCercle;num_cercle++){//scan polar half diameter
        theta=(num_cercle)*delta_theta;//theta=0, singularité donc on exclut num_cercle==0
        //delta_phi=surface_moyenne_pt/delta_theta; //donc d ?
        delta_curv_phi=dist_recalc;//delta abscisse curviligne sur le cercle concentrique
        cout<<"distance delta_curv_phi"<<delta_curv_phi<<endl;
        int nbPtPhi=round(2*M_PI*sin(theta)/(delta_curv_phi));//M_PHI=nbPoint dans le cercle concentrique actuel=perimetre/delta_phi
        cout<<"---------------------"<<endl;
        cout<<"nbPtcircle="<<nbPtPhi<<endl;
        cout<<"rayon cercle="<<Vmax*sin(theta)<<endl;
        for(int numPt=0;numPt<nbPtPhi;numPt++){//scan  concentric circles
            //phi=numPt*delta_phi;
            phi=2*M_PI*(numPt)/nbPtPhi;
            Vout_table[numHolo].x=Vmax*sin(theta)*cos(phi);
            Vout_table[numHolo].y=Vmax*sin(theta)*sin(phi);
            kiz[numHolo]=Vmax*cos(theta);
            numHolo++;
        }
    }
}

void scan_etoile(double Vmax, int nbHolo, int nbAxe, vector<float2D> &Vout_table)
{
double delta_theta=M_PI/nbAxe;
int nbPtAxe=floor((double)nbHolo/(double)nbAxe);//integer division, rounded to the closest lesser value
vector<double> kiz(nbHolo);
vector<int> nbPointAxe(nbAxe);
int numHolo=0;
int nbPtLeft=nbHolo-(nbPtAxe*nbAxe);//nb pt left if non integer
vector<double> theta(nbAxe);
for(int cpt=0;cpt<nbAxe;cpt++){
theta[cpt]=cpt*delta_theta;
if(cpt>=(nbAxe-nbPtLeft))//if the remainder of the integer division is not zero, put the N remaining points on the  N last  axes
    nbPointAxe[cpt]=nbPtAxe+1;
else nbPointAxe[cpt]=nbPtAxe;
}
for(int cpt=0;cpt<nbAxe;cpt++){
    double delta_rho=2*Vmax/nbPointAxe[cpt];
    double rho=0;
    for(int numPt=0; numPt<nbPointAxe[cpt];numPt++){
        rho=-Vmax+numPt*delta_rho;//parcours rayon à theta fixé
        Vout_table[numHolo].x=rho*cos(theta[cpt]);
        Vout_table[numHolo].y=rho*sin(theta[cpt]);
        kiz[numHolo]=sqrt(Vmax*Vmax-pow(Vout_table[numHolo].x,2)-pow(Vout_table[numHolo].y,2));
        numHolo++;
    }
}
cout<<"nb point placés="<<numHolo<<endl;
}
