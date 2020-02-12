#include "fonctions.h"



///####################fonction################"
void SAV(double *var_sav, int taille, char *chemin, enum PRECISION precision, char options[])
{
        FILE *fichier_ID;
        fichier_ID= fopen(chemin, options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision) {
        case DOUBLE: //64 bit

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        double tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case FLOAT://32 bits float

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        float tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;

        case INT: //32 bit signé

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case UINT://32 bit non signé

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        unsigned int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case CHAR: //8 bits

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        char tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        default:
                break;
        }

        fclose(fichier_ID);
}

void changeDim2D(double* tab, double* tabFinal, Var2D dimInit, Var2D dimFin)
{
        int diffX=round(dimFin.x-dimInit.x)/2;
        int diffY=round(dimFin.y-dimInit.y)/2;

        for(int x=0; x<dimInit.x; x++) {
                for( int y=0; y<dimInit.y; y++) {
                        int cpt_init=y*dimInit.x+x;
                        int cpt_final=(y+diffY)*dimFin.x+x+diffX;
                        tabFinal[cpt_final]=tab[cpt_init];
                }
        }
}



/*void concatener(char chaine1[], char chaine2[], char resultat[])//fonction limitée mais intéret car copie locale (pas besoin de tampon)
{

        //récupérer la longueur du resultat


        size_t LongrChaine1=strlen(chaine1);
        size_t LongrChaine2=strlen(chaine2);
        short unsigned int LongrChaine=LongrChaine1+LongrChaine2;
        char tampon[LongrChaine];//+1 pour zero final?

        strcpy(tampon,chaine1);
        strcat(tampon,chaine2);

        for(short unsigned int cpt=0; cpt<LongrChaine+1; cpt++)
                resultat[cpt]=tampon[cpt];

}*/

//Interpolation3D : attend un volume "troué" et les dimensions en x,y,et z

void interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z)
{
        int z=0;
        // int cpt=x+y*dv0rf+z*taille_y*dv0rf;
        int z_min=taille_z;
        int z_max=0;

        // cout<< "valeur papillon : "<< volume_interp_3D[cpt] << endl;
        for (int x=0; x < taille_x; x++) { //on balaye l'image, référentiel avec centre au milieu
                for (int y=0; y<taille_y; y++) { //on balaye l'image,
                        //printf("y:%i  x: : %i ",y,x);
                        //cpt=x+y*taille_x+z*taille_y*taille_x;
                        // 				for (int z = 0; z < taille_z; z++)//on balaye l'image,
                        // 				{cpt=x+y*taille_x+z*taille_y*taille_x;
                        // // 					if(sup_redon[cpt]==599)
                        // // 					{
                        // // 					cout << "x : " << x <<"y : "<< y<< "z : "<<z<<"cpt : "<<cpt<<"sup_redon : " << sup_redon[cpt] << endl;
                        // // 					}
                        //
                        // 				}
                        z=0;
                        //printf("cpt : %i, sup redon : %d\n",x+y*taille_x+z*taille_y*taille_x,volume_interp_3D[x+y*taille_x+z*taille_y*taille_x]);
                        while(z<taille_z) { // pas de boucle for car l'indice est variable
                                float valeur_test=volume_interp_3D[x+y*taille_x+z*taille_y*taille_x];
                                if(valeur_test!=0 && z<taille_z) { //si valeur !=0, alors borne inférieure mais z_min ne peut valoir 19.
                                        //printf("coucou1 sup_redon : %f,x,y,z : %i,%i,%i\n",volume_interp_3D[x+y*taille_x+z*taille_y*taille_x],x,y,z);
                                        z_min=z; //on a trouvé z_min
                                        //printf("x,y: %i,%i z : %i, sup_redon : %f\n",x,y,z, volume_interp_3D[x+y*taille_x+z*taille_y*taille_x]);
                                        z_max=z_min+1; //initialiser z_max au point suivant (sinon while jamais verifie)

                                        while( z_max<taille_z && volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]==0) { //soit trouver prochain z_max, soit fin de colonne
                                                //printf("coucou2\n,%i,%i,%i,sup_redon : %f\n",x,y,z_max,volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]);
                                                z_max=z_max+1;
                                                z=z_max;

                                                if(z_max==taille_z-1) {
                                                        //printf("fin de colonne\n");
                                                }
                                        }
                                        if(z_max!=taille_z && z_max-z_min>1 && volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]!=0) { //il faut au moins un trou pour interpoler
                                                //printf("interpolation x,y,zmin, z_max: %i,%i,%i, %i\n", x,y,z_min,z_max);
                                                //y=ax+b-> interpolation linéaire
                                                double a=(volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]-volume_interp_3D[x+y*taille_x+z_min*taille_y*taille_x])/(z_max-z_min);
                                                double b=volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]-a*z_max;
                                                //printf("hello: z_min: %i, zmax:%i\n",z_min,z_max);
                                                for(int cpt_z=z_min+1; cpt_z<z_max; cpt_z++)

                                                {
                                                        volume_interp_3D[x+y*taille_x+cpt_z*taille_y*taille_x]=a*cpt_z+b;
                                                        //printf("interpolaion finie\n");
                                                }
                                                // 						x;
                                                // 						volume_interp_3D[x+y*taille_x+z_min*taille_y*taille_x];
                                                // 						z_min;
                                                // 						volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x];
                                                // 						z_max;

                                                z=z_max-1; // nouveau compteur=ancienne borne sup (-1 car z++)
                                        }
                                        // redémarrer le compteur de ligne à z_max
                                }
                                z++;
                        }

                }
        }
}
//##############§FIN INTERP3D#####################################################
///#######masque pour &écraser jumeau############""
void antigaussienne(double *tab, int Tx, int sigma, float A, int Exy)
{

        int corr_paire=0;
        if(Tx%2==0) {
                cout<<"taille Tx="<<Tx<<" paire : le masque ne serait pas centré!"<<endl;
                corr_paire=1;
        }

        if(sigma==0)
                sigma=1;
        short  int x,y, Tinf=-round(Tx/2),Tsup=round(Tx/2);
        short unsigned int cptx,cpty;
        float Ex,Ey,sigmax,sigmay;
        Ex=Exy;
        Ey=Exy;
        sigmax=sigma;
        sigmay=sigma;
        if(Tx==1)
                tab[0]=0;
        else {
                for(x=Tinf; x<Tsup+1-corr_paire; x++) {
                        cptx=x+Tsup;
                        for( y=Tinf; y<Tsup+1-corr_paire; y++) {
                                cpty=y+Tsup;
                                tab[cpty*Tx+cptx]=1-A*exp(-(pow((x-Ex),4)/(2*sigmax*sigmax)+pow((y-Ey),4)/(2*sigmay*sigmay)));
                        }
                }
        }

}

//#############################################################"
double *tukey2D(int dimx,int dimy, float alpha)
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
}

void ecrire_rapport(int NXMAX,float rayon,float Rf, float K, int DIMX_CCD2,int coin_x, int coin_y,short int precision_exportation,string chemin,int nb_proj,float n1,float NA,float Tp, int G)
{
        time_t date;
        time(&date);
        char nom_rapport[12];
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
        for(int pixel=0; pixel<2*t_image*2*t_image; pixel++) {
                masque[pixel]=0;
        }
        for(int x=centreX-t_mask/2; x<centreX+t_mask/2; x++) {
                for(int y=centreY-t_mask/2; y<centreY+t_mask/2; y++) {
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

}

/////////////////////////////////////////////////////////////////////////////////////
//Fonction transferant les pixels de l'image 2D de taille (nx,ny) dans un tableau 1D de taille (nx*ny)
// nécessite la bibliotheque imageMagick
/////////////////////////////////////////////////////////////////////////////////////

/*void rempli_tableau(unsigned char *finalArray, string path, Var2D coin, Var2D taille)
{

}*/

/////////////////////////////////////////////////////////////////////////////////////
//fonction de calcul TF2D
/////////////////////////////////////////////////////////////////////////////////////
void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y)
{
        fftw_plan_with_nthreads(3);
        int N=taille_x*taille_y;
        fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
        fftw_plan p;
        //Réservation memoire
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        //Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<N; cpt++) {
                in[cpt][0]=entree_reelle[cpt];
                in[cpt][1]=entree_imag[cpt];
        }
        //calcul du plan, parametre servant a calculer et optimiser le FFT
        p=fftw_plan_dft_2d( taille_x,  taille_y, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p); /* repeat as needed */

        for(int cpt=0; cpt<(N); cpt++) {

                fft_reel[cpt]=out[cpt][0]; //division par N^2 pour normaliser la fftw qui n'est pas normalisée
                fft_imag[cpt]=out[cpt][1];
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}

void TF2D_INV(double entree_reelle[],double entree_imag[],double sortie_reelle[],double sortie_imag[],int taille_x,int taille_y)
{
        fftw_plan_with_nthreads(3);
        int N=taille_x*taille_y;
        fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
        fftw_plan p;
        //Réservation memoire
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        //Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<N; cpt++) {
                in[cpt][0]=entree_reelle[cpt];
                in[cpt][1]=entree_imag[cpt];
        }
        //calcul du plan, parametre servant a calculer et optimiser le FFT
        p=fftw_plan_dft_2d( taille_x,  taille_y, in, out,FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(p); /* repeat as needed */

        for(int cpt=0; cpt<(N); cpt++) {

                sortie_reelle[cpt]=out[cpt][0];
                sortie_imag[cpt]=out[cpt][1];
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}

/////////////////////////////////////////////////////////////////////////////////////
//fonction de circshift

/////////////////////////////////////////////////////////////////////////////////////
void circshift2(double* entree, double* result, Var2D dim,Var2D decal)
{
        //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo

        decal.y=decal.y%dim.y;
        decal.x=decal.x%dim.x;

        for(int yi=0; yi<dim.y; yi++) {
                int yi_decal=yi+decal.y;
                //si dépassement des dimensions
                if(yi+decal.y>dim.y-1)
                        yi_decal=-dim.y+yi+decal.y;
                //si décage négatif
                if(yi+decal.y<0)
                        yi_decal=dim.y+yi+decal.y;
                int hauteur_decal=yi_decal*dim.x;
                int hauteur=yi*dim.x;
                for(int xi=0; xi<dim.x; xi++) {
                        int xi_decal=xi+decal.x;

                        if(xi+decal.x>dim.x-1)
                                xi_decal=-dim.x+xi+decal.x;
                        if(xi+decal.x<0)
                                xi_decal=dim.x+xi+decal.x;

                        int pos1D_init=hauteur+xi;
                        int pos1D_shift=hauteur_decal+xi_decal;//nouvelle position après décalage
                        result[pos1D_shift]=entree[pos1D_init];
                }
        }


}
void circshift3(double* entree, double* result, Var2D dim,Var2D decal)
{
        //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
        decal.y=decal.y%dim.y;
        decal.x=decal.x%dim.x;

        for(int yi=0; yi<decal.y; yi++) {
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


/*double* circshift(double* entree, int dimx,int dimy,int decal_x,int decal_y)
{
        double * entree_shift=new double[dimx*dimy];//fuite memoire!
        for(int yi=0; yi<decal_y; yi++) {
                for(int xi=0; xi<decal_x; xi++)

                {
                        int pixel=yi*dimx+xi;
                        int pixel_shift=(yi+decal_y)*dimx+xi+decal_x;
                        //1er quadrant vers 4 eme
                        entree_shift[pixel_shift]=entree[pixel];
                        //4 eme quadrant vers 1er
                        entree_shift[pixel]=entree[pixel_shift];
                        //2eme vers 3eme
                        entree_shift[(yi+decal_y)*dimx+xi]=entree[pixel+decal_x];
                        //3eme vers 2eme
                        entree_shift[pixel+decal_x]=entree[(yi+decal_y)*dimx+xi];
                }
        }

        return entree_shift;

}*/
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D)
{
        decal3D.x=decal3D.x%dimFinal3D.x;//élmiiner les "modulos"
        decal3D.y=decal3D.y%dimFinal3D.y;
        decal3D.z=decal3D.z%dimFinal3D.y;

        unsigned short int xi,yi,zi=0;
        short int x2,y2,z2=0; //signé car une fois décalé, peuvent être négatifs!
        const unsigned int taille_plan =dimFinal3D.x*dimFinal3D.y;

        for(zi=0; zi<dimFinal3D.z; zi++) {
                if(zi+decal3D.z>dimFinal3D.z-1) { //dépassement à droite
                        z2=zi+decal3D.z-dimFinal3D.z;
                } else {
                        if(zi+decal3D.z<0) { //dépassement à gauche
                                z2=dimFinal3D.z+(decal3D.z+zi);
                        } else {
                                z2=zi+decal3D.z;
                        }
                }
                int nb_pixelz_decal=z2*taille_plan;
                unsigned int nb_pixelz=zi*taille_plan;
                for(yi=0; yi<dimFinal3D.y; yi++) {
                        if(yi+decal3D.y>dimFinal3D.y-1) { //dépassement à droite
                                y2=yi+decal3D.y-dimFinal3D.y;
                        } else {
                                if(yi+decal3D.y<0) { //dépassement à gauche
                                        y2=dimFinal3D.y+(decal3D.y+yi);
                                } else {
                                        y2=yi+decal3D.y;
                                }
                        }
                        int nb_lignes=yi*dimFinal3D.x;
                        int nb_lignes_decal=y2*dimFinal3D.x;

                        for(xi=0; xi<dimFinal3D.x; xi++) {
                                if(xi+decal3D.x>dimFinal3D.x-1) { //dépassement à droite
                                        x2=xi+decal3D.x-dimFinal3D.x;
                                } else {
                                        if(xi+decal3D.x<0) { //dépassement à gauche
                                                x2=dimFinal3D.x+(decal3D.x+xi);
                                        } else {
                                                x2=xi+decal3D.x;
                                        }
                                }

                                // pixel=zi*dimFinal3D.x*dimFinal3D.y+yi*dimFinal3D.x+xi;
                                //pixel_decal=z2*dimFinal3D.x*dimFinal3D.y+y2*dimFinal3D.x+x2;
                                //volume3D_shift[pixel_decal]=volume3D[pixel];
                                //volume3D_shift[z2*dimFinal3D.x*dimFinal3D.y+y2*dimFinal3D.x+x2]=volume3D[zi*dimFinal3D.x*dimFinal3D.y+yi*dimFinal3D.x+xi];
                                volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2]=volume3D[nb_pixelz+nb_lignes+xi];
                                //memcpy(volume3D_shift+(nb_pixelz_decal+nb_lignes_decal+x2), volume3D+(nb_pixelz+nb_lignes+xi), 8); commande ok mais plus lent !?

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

        for(int x=centreX-t_mask/2; x<centreX+t_mask/2; x++) {
                for(int y=centreY-t_mask/2; y<centreY+t_mask/2; y++) {
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

        for(int x=centreX-t_mask/2; x<centreX+t_mask/2; x++) {
                for(int y=centreY-t_mask/2; y<centreY+t_mask/2; y++) {
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



// soient deux images allouées en row-major
// découpe une fenêtre de img_src et l'écrit dans img_dst
// il faut préciser la taille des deux images et la position du coin supérieur gauche de la fenêtre

void crop_window(double* img_src, size_t src_dimx, size_t src_dimy, double* img_dst, size_t dst_dimx, size_t dst_dimy, size_t ul_corner_x, size_t ul_corner_y)
{
        assert(dst_dimx < src_dimx);
        assert(dst_dimy < src_dimy);

        assert((ul_corner_x + dst_dimx) <= src_dimx);
        assert((ul_corner_y + dst_dimy) <= src_dimy);

        double *start_src = img_src + (ul_corner_y * src_dimx);

        double *curs_src = start_src;
        size_t pos_dst = 0;
        size_t i, j;

        for (i = 0; i < dst_dimy; i++) {
                for (j = ul_corner_x; j < ul_corner_x + dst_dimx; j++) {
                        img_dst[pos_dst] = curs_src[j];
                        pos_dst++;
                }
                curs_src += src_dimx;
        }
}
