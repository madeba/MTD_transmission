/***************************************************************
 * Name:      Gui_TomoMain.h
 * Purpose:   Defines Application Frame
 * Author:    Matthieu Debailleul ()
 * Created:   2017-03-28
 * Copyright: Matthieu Debailleul ()
 * License:
 **************************************************************/
#include <iostream>

#ifndef GUI_TOMOMAIN_H
#define GUI_TOMOMAIN_H

#ifndef WX_PRECOMP
    #include <wx/wx.h>
#endif

#include <fstream>
#include <sstream>
#include <math.h>
#include <wx/frame.h>
#include<vector>

class Gui_TomoFrame: public wxFrame
{   ///La classe de la fenêtre principale est publique
    public:

        Gui_TomoFrame(wxFrame *frame, const wxString& title);
        ~Gui_TomoFrame();
         //float extract_val(std::string token,  std::string chemin_fic);
         std::string extract_string(std::string token,  std::string chemin_fic);

         std::string extract_string(std::string token,  std::string chemin_fic,std::string string_defaut);
         void sav_val(std::string chemin_fic,std::vector<std::string> &tab_val);
         float extract_val(std::string token,  std::string chemin_fic);
         void init_tab_val(std::string chemin_fic, std::vector<std::string> &tab_val);
         void modif_tab_val(std::string token,std::string string_valeur_token,std::vector<std::string> &tab_val);
         void refreshValue(std::string chemin_fic, std::string chemin_recon);
         void sav_all();
         void sav_all_pc_acquis();
         void sav_all_pc_recon();

  protected:
      ///ENUM servant à identifier les boutons et menus pour la table d'évènement
        enum
        {
            idMenuQuit = 1000,
            wxID_EXIT,
            idMenuAbout,
            idBoutonManip,
            idBoutonImage,
            idBoutonTraitement,
            idBoutonReconstruction,
            id_boolBorn,
            idBoolVolkov,
            idBoolDeroul,
            idBoolAber,
            idBoolAmpliRef,
            idBoolExportOTF,
            idBoutonOpenDir,
            idBoutonOpenDirResult,
            idBoutonRecons,
            idBoutonOpenMask,
            idBoutonRecalc,
            idBoutonSAV,
            idMenuSAV,
            idBoutonGPS,
            idComboScan,
        };

        double Vxmax,Vymax,Vxmin,Vymin;
        int fx0,fy0,nbHolo,nbCirclesAnnular;
        unsigned int  DIM_FINAL=512, NXMAX=110;
        bool b_BORN=0, b_VOLKOV=1, b_DEROUL=1, b_ABER=1, b_EXPORT_OTF=0;
       // std::vector<std::string> fichier_tmp;
        std::vector<std::string> tab_val_recon;
        std::vector<std::string> tab_val_manip;
       // std::vector<std::string> tab_val_manip_pc_acquis;///pc acquisition
        std::vector<std::string> tab_val_gui_conf;


    // textFicManip=new wxTextCtrl(zone1,-1,"",wxPoint(160,195),wxSize(100,20));


        ///Déclaration des attributs de type Boutons et menus en protected
        std::string chemin_config_GUI,chemin_recon, chemin_config_manip,chemin_rep_acquis,repertoire_config_defaut,chemin_result,repertoire_config;
        //std::string repertoire_config_pc_acquis,chemin_config_manip_pc_acquis;
        wxButton* BoutonManip,*BoutonImage;
        wxButton* BoutonTraitement;
        wxButton* BoutonRecons;
        wxButton* BoutonOpenDir,*BoutonOpenDirResult,*BoutonOpenMask;

        wxButton* BoutonGPS;
        wxBitmapButton *BoutonRecalc, *BoutonSAV;
        wxCheckBox *m_born,*m_volkov,*m_Aber,*m_Deroul,*m_AmpliRef,*m_export_OTF;

        wxComboBox* ComboBox_Scan;

        wxStaticText* t;
        wxTextCtrl *editX,*editY,*editResult,*editDirAcquis,*editFicMask,*editFicManip,*editDirResultAcquis;//champ de texte editable
        wxTextCtrl *editNbHolo,*editnbCirclesAnnular,*editDimFinal;
        wxTextCtrl *editCX,*editCY,*editNXMAX,*editVxmin,*editVymin,*editVxmax,*editVymax,*editDeltaVx,*editDeltaVy,*editNAcondLim;
        wxTextCtrl *editDnMin,*editDnMax, *editIteration;
        wxTextCtrl *editKappaMin,*editKappaMax, *editKappa;

        wxStaticText *textX,*textY,*textResult,*textDirAcquis,*textFicManip;
        wxStaticText *textDnMin,*textDnMax;
        wxStaticText *textKappaMin,*textKappaMax;
        wxStaticText *titre_Acquis,*titre_Pretraitement, *titre_Recons,*titre_HorsAxe,*titre_Balayage,*titre_Interval_Balayage;//texte des champs
        wxStaticText *textNbHolo,*textNbCirclesAnnular,*textScanPattern, *textDimFinal,*textCX,*textCY,*textNXMAX,*textVxmin,*textVymin,*textVxmax,*textVymax,*textDeltaVy,*textNAcondLim,*textDeltaVx,*textOffx,*textOffy;//texte des champs
        wxStaticText *titre_GPS,*titre_interval_indice,*titre_iterationGPS,*textIteration;
        wxStaticText  *titre_interval_kappa;


  private:
        void OnClose(wxCloseEvent& event);
        void OnQuit(wxCommandEvent& event);
        void OnAbout(wxCommandEvent& event);
        void OnComboScan(wxCommandEvent& event);
        void OnBoutonManip(wxCommandEvent& event);
        void OnBoutonImage(wxCommandEvent& event);
        void OnBoutonTraitement(wxCommandEvent& event);
        void OnBoutonRecons(wxCommandEvent& event);
        void OnBoutonBorn(wxCommandEvent& event);
        void OnBoutonExportOTF(wxCommandEvent& event);
        void OnBoutonVolkov(wxCommandEvent& WXUNUSED(event));
        void OnBoutonAber(wxCommandEvent& event);
        void OnBoutonAmpliRef(wxCommandEvent& event);
        void OnBoutonDeroul(wxCommandEvent& event);
        void OnBoutonOpenDir(wxCommandEvent& event);
        void OnBoutonOpenDirResult(wxCommandEvent& event);
        void OnBoutonOpenMask(wxCommandEvent& event);
        //void OnBoutonRecalc(wxCommandEvent& WXUNUSED(event));
        void OnBoutonSAV(wxCommandEvent& WXUNUSED(event));
        void OnBoutonGPS(wxCommandEvent& event);

        DECLARE_EVENT_TABLE()
};


#endif // GUI_TOMOMAIN_H
