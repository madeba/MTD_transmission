#include "fonctions.h"
#include <tiffio.h>
//#include <cxcore.h>
#include <chrono>
//#include <highgui.h>
//#include <cv.h>
#define PI 3.14159
//using namespace cv;



 //extrait la valeur (entier) de la chaine "token" du fichier "chemin_fic"

string extract_string(std::string token,  std::string chemin_fic)
{
    ifstream fichier(chemin_fic.c_str(), ios::in);  // on ouvre en lecture
    string ligne,motcle,valeurMot,separ=" ";
    string valeur;
    vector<std::string> tokens;

    if(fichier)  // si l'ouverture a fonctionné
    {
        while(!fichier.eof()){
            getline(fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire
            tokens.push_back(ligne);
        }
    }
    else
        cerr << "Impossible d'ouvrir le fichier !"<< chemin_fic<< endl;

    int nb_tok=tokens.size();
    for(int cpt=0;cpt<nb_tok;cpt++){
        ligne=tokens[cpt];
        if(ligne!=""){
            int pos_separ=ligne.find(separ);
            int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin)
            if(motcle==token){
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
        while(!fichier.eof()){
            getline(fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire
            tokens.push_back(ligne);
        }
    }
    else
        cerr << "Impossible d'ouvrir le fichier !"<< chemin_fic<< endl;

    int nb_tok=tokens.size();
    for(int cpt=0;cpt<nb_tok;cpt++){
        ligne=tokens[cpt];
        if(ligne!=""){
            int pos_separ=ligne.find(separ);
            int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin)
            if(motcle==token){
            valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));
            cout<<motcle<<"="<<valeurMot<<endl;
            valeur=atof(valeurMot.c_str());
            }
        }
    }
    fichier.close();
    return valeur;
}


float ecrire_val(string token, float valeur_token, string chemin_fic)
{    ///on sauvegarde le fichier, au cas où...
     ifstream src(chemin_fic.c_str() ,ios::binary);
     ofstream dst(chemin_fic+"_SAV" ,ios::binary);
     dst<<src.rdbuf();
     src.close();
     dst.close();

     string ligne, motcle,valeurMot,nlle_ligne,string_val_token,separ=" ";
     std::ostringstream ss;
     ifstream flux_fichier(chemin_fic.c_str(), ios::in);//ouverture en lecture ecriture
     vector<std::string> fichier_tmp;
     if(flux_fichier)  // si l'ouverture a réussi
     {    while(!flux_fichier.eof()){
            getline(flux_fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire (premier caractère ligne!=#)
            {
                int pos_separ=ligne.find(separ);//trouver l'espace
               // cout<<"position de l'espace="<<pos_separ<<endl;
//                int long_separ=separ.length();
                motcle = ligne.substr(0, pos_separ);//copie une portion de lachaine entre 0 et pos_separ=)isole le mot-clé
                if(motcle==token){
                ss << valeur_token;
                string_val_token=ss.str();//convertir float to string //valeur_token
                nlle_ligne=token+" "+string_val_token;
                //cout<<"nouvelle ligne"<<nlle_ligne<<endl;
                fichier_tmp.push_back(nlle_ligne);
                }
                else
                    fichier_tmp.push_back(ligne);//le motclé n'est pas le bon, on ecrit juste la ligne
            }
            else{///pas de mot clé :  ligne de commentaire
            fichier_tmp.push_back(ligne);
            }
        }
     }
     else  // sinon erreur
            cerr << "Erreur d'ouverture du fichier"<<chemin_fic << endl;
    flux_fichier.close();
    ofstream flux_ecriture(chemin_fic.c_str(), ios::out);  // ouverture en écriture
    for(int cpt=0;cpt<fichier_tmp.size();cpt++)
    flux_ecriture<<fichier_tmp[cpt]<<endl;//ecriture dans le fichier
    flux_ecriture.close();  // je referme le fichier
     return 0;
}


