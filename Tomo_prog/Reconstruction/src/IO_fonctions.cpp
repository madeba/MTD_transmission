#include "IO_fonctions.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <fstream>
#include <sstream>
using namespace std;
using namespace cv;



///charger en mémoire le tableau de mot-clés à partir d'un fichier de config situé à chemin_fic.
void init_tab_val(string chemin_fic, vector<string> &tab_val)
{
    //effacer l'ancienne version->il faudrait initialiser dans le corps du programme!
    while (!tab_val.empty())
        tab_val.pop_back();

    string ligne;
    ifstream flux_fichier(chemin_fic.c_str(), ios::in);//ouverture en lecture ecriture

    if(flux_fichier)  // si l'ouverture a réussi
    {
        while(!flux_fichier.eof())
        {
            getline(flux_fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            tab_val.push_back(ligne);
        }
    }
    else  // sinon erreur
        cerr << "Erreur d'ouverture du fichier"<<chemin_fic << endl;
    flux_fichier.close();
}


//modifier en mémoire la valeur des tokens (mots-clés)
void modif_tab_val(string token,string string_valeur_token,vector<string> &tab_val)
{
    std::vector<std::string> fichier_new;
    string ligne,motcle,string_val_token,nlle_ligne,separ=" ";

    for(size_t cpt=0; cpt<tab_val.size(); cpt++)
    {
        ligne=tab_val[cpt];

        if(ligne!="#")//si pas ligne de commentaire (premier caractère ligne!=#)
        {
            int pos_separ=ligne.find(separ);//trouver l'espace (séparateur mot-clé valeur)
            //cout<<"position de l'espace="<<pos_separ<<endl;
            // int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//copie une portion de lachaine entre 0 et pos_separ=isole le mot-clé
            if(motcle==token)
            {
                nlle_ligne=token+" "+string_valeur_token;
                //cout<<"nouvelle ligne"<<nlle_ligne<<endl;
                fichier_new.push_back(nlle_ligne);
            }
            else
                fichier_new.push_back(ligne);//le motclé n'est pas le bon, on ecrit juste la ligne
        }
        else ///caratère=# :  ligne de commentaire
        {
            fichier_new.push_back(ligne);
        }
    }
    //effacer l'ancienne version
    while (!tab_val.empty())
    {
        tab_val.pop_back();
    }
    //insérer la nouvelle version
    for(size_t cpt=0; cpt<fichier_new.size(); cpt++)
    {
        tab_val.push_back(fichier_new[cpt]);
    }
    //effacer le tampon
    for(size_t cpt=0; cpt<fichier_new.size(); cpt++)
    {
        fichier_new.pop_back();
    }
}

void sav_val(string chemin_fic,vector<string> &tab_val)
{
    ///on sauvegarde le fichier, au cas où...
    ifstream src(chemin_fic.c_str(),ios::binary);
    ofstream dst(chemin_fic+"_SAV",ios::binary);
    dst<<src.rdbuf();
    src.close();
    dst.close();
    ofstream flux_ecriture(chemin_fic.c_str(), ios::out);  // ouverture en écriture
    for(size_t cpt=0; cpt<tab_val.size(); cpt++){
        flux_ecriture<<tab_val[cpt]<<endl;//ecriture dans le fichier
    }
    flux_ecriture.close();  //  fermer
}
