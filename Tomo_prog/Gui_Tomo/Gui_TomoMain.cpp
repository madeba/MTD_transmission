/***************************************************************
 * Name:      Gui_TomoMain.cpp
 * Purpose:   Code for Application Frame
 * Author:    Matthieu Debailleul ()
 * Created:   2017-03-28
 * Copyright: Matthieu Debailleul ()
 * License:
 **************************************************************/
#include "Gui_TomoApp.h"
#include "fonctions.h"
#ifdef WX_PRECOMP
#include "wx_pch.h"
#endif

#ifdef __BORLANDC__
#pragma hdrstop
#endif //__BORLANDC__

#include "Gui_TomoMain.h"
#include<algorithm>//fonction remove() sur string
#include <wx/notebook.h>///gestion des onglets
#include "wx/popupwin.h"
#include "IO_fonctions.h"
using namespace std;

//helper functions
enum wxbuildinfoformat
{
    short_f, long_f
};

wxString wxbuildinfo(wxbuildinfoformat format)
{
    wxString wxbuild(wxVERSION_STRING);

    if (format == long_f )
    {
#if defined(__WXMSW__)
        wxbuild << _T("-Windows");
#elif defined(__WXMAC__)
        wxbuild << _T("-Mac");
#elif defined(__UNIX__)
        wxbuild << _T("-Linux");
#endif

#if wxUSE_UNICODE
        wxbuild << _T("-Unicode build");
#else
        wxbuild << _T("-ANSI build");
#endif // wxUSE_UNICODE
    }

    return wxbuild;
}
using namespace std;
///table d'évènement liant identifiant des boutons et methode
BEGIN_EVENT_TABLE(Gui_TomoFrame, wxFrame)
    EVT_CLOSE(Gui_TomoFrame::OnClose)
    EVT_MENU(wxID_EXIT, Gui_TomoFrame::OnQuit)
    EVT_MENU(idMenuAbout, Gui_TomoFrame::OnAbout)
    EVT_BUTTON(idBoutonManip, Gui_TomoFrame::OnBoutonManip)//Bouton Acquisition
    EVT_BUTTON(idBoutonImage, Gui_TomoFrame::OnBoutonImage)//Bouton Acquisition
    EVT_BUTTON(idBoutonTraitement, Gui_TomoFrame::OnBoutonTraitement)
    EVT_BUTTON(idBoutonRecons, Gui_TomoFrame::OnBoutonRecons)
    // EVT_BUTTON(idBoutonCalcul, Gui_TomoFrame::OnBoutonCalcul)
    EVT_CHECKBOX(id_boolBorn,Gui_TomoFrame::OnBoutonBorn)
    EVT_CHECKBOX(idBoolExportOTF,Gui_TomoFrame::OnBoutonExportOTF)
    EVT_CHECKBOX(idBoolVolkov,Gui_TomoFrame::OnBoutonVolkov)
    EVT_CHECKBOX(idBoolDeroul,Gui_TomoFrame::OnBoutonDeroul)
    EVT_CHECKBOX(idBoolAber,Gui_TomoFrame::OnBoutonAber)
    EVT_CHECKBOX(idBoolAmpliRef,Gui_TomoFrame::OnBoutonAmpliRef)
    EVT_BUTTON(idBoutonOpenDir, Gui_TomoFrame::OnBoutonOpenDir)
    EVT_BUTTON(idBoutonOpenDirResult, Gui_TomoFrame::OnBoutonOpenDirResult)
    EVT_BUTTON(idBoutonOpenMask, Gui_TomoFrame::OnBoutonOpenMask)
    EVT_BUTTON(idBoutonGPS, Gui_TomoFrame::OnBoutonGPS)
    //EVT_BUTTON(idBoutonRecalc, Gui_TomoFrame::OnBoutonRecalc)
    EVT_BUTTON(idBoutonSAV, Gui_TomoFrame::OnBoutonSAV)
    EVT_MENU(idMenuSAV, Gui_TomoFrame::OnBoutonSAV)
    EVT_COMBOBOX(idComboScan, Gui_TomoFrame::OnComboScan)//evenement pour la combo box
END_EVENT_TABLE()


///constructeur de la fenêtre, on y met toutes les init de boutons/menus (faisable aussi dans le App.main).
Gui_TomoFrame::Gui_TomoFrame(wxFrame *frame, const wxString& title)
    : wxFrame(frame, -1, title)
{
#if wxUSE_MENUS
    // create a menu bar
    wxMenuBar* mbar = new wxMenuBar();
    wxMenu* fileMenu = new wxMenu(_T(""));

    fileMenu->Append(idMenuSAV, _("&Sauver\tAlt-S"), wxT("Sauver les paramètres"));
    fileMenu->Append(wxID_EXIT, _("&Quitter\tAlt-F4"), _("Quitter l'application"));
    mbar->Append(fileMenu, _("&File"));

    wxMenu* helpMenu = new wxMenu(_T(""));
    helpMenu->Append(idMenuAbout, _("&A propos\tF1"), _("Information sur l'application"));
    mbar->Append(helpMenu, _("&Aide"));

    SetMenuBar(mbar);

#endif // wxUSE_MENUS

#if wxUSE_STATUSBAR
    // create a status bar with some information about the used wxWidgets version
    CreateStatusBar(2);
    SetStatusText(_("Tomo Mithra 0.1!"),0);
    SetStatusText(wxbuildinfo(short_f), 1);
#endif // wxUSE_STATUSBAR



///Initialiser les deux tableaux stockant les valeurs de recon.txt et config_manip.txt en mémoire

    string home=getenv("HOME");
    string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
    chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
    cout<<"lecture des chemins dans le fichier "<<chemin_config_GUI<<endl;

    repertoire_config_defaut=extract_string("CHEMIN_CONFIG_PC_ACQUIS",home+fin_chemin_gui_tomo);///repertoire config_manip.txt si pc acquisition (defaut)

    //repertoire_config=extract_string("CHEMIN_CONFIG_ACQUIS",home+fin_chemin_gui_tomo);///repertoire config_manip.txt : Inutile ?

    chemin_rep_acquis=extract_string("CHEMIN_ACQUIS",home+fin_chemin_gui_tomo);///repertoire des données
    if(IsPathOk(chemin_rep_acquis)==0) cout<<"\033[1;31m No Acquisition directory "<<chemin_rep_acquis <<". Please create it \033[1;0m"<<endl;//https://en.wikipedia.org/wiki/ANSI_escape_code#Colors

    repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);///repertoire config_manip.txt
    if(IsPathOk(repertoire_config)==0) cout<<"\033[1;31m No config directory "<<repertoire_config <<". Please create it \033[1;0m"<<endl;
    //cout<<"Repertoire_config="<<repertoire_config<<endl;
    chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);///repertoire des résultats
    if(IsPathOk(chemin_result)==0) cout<<"\033[1;31m No results directory "<<chemin_result <<". Please create it \033[1;0m"<<endl;

    std::remove(chemin_result.begin(), chemin_result.end(), ' ');

    chemin_config_manip=repertoire_config+"/config_manip.txt";
    cout<<"fichier config_manip.txt: "<<chemin_config_manip<<endl;
    chemin_recon=repertoire_config+"/recon.txt";
    cout<<"fichier recon.txt : "<<chemin_recon<<endl;
    b_PC_ACQUIS=stof(extract_string("PC_ACQUIS",chemin_config_GUI),0);
    ///tester l'existence du fichiers de config:  config_manip.txt dans le chemin indiqué par l'utilisateur
    ifstream fichier_config_manip(chemin_config_manip.c_str(), ios::in);
    if(!fichier_config_manip)
    {
        cout<<endl<<"fichier ::::::::::::::: "<<chemin_config_manip<<" inexistant"<<endl;
        ///on copie le fichier defaut dans le repertoire de config...
        ifstream src(repertoire_config_defaut+"config_manip.txt",ios::binary);
        ofstream dst(chemin_config_manip,ios::binary);
        dst<<src.rdbuf();
        src.close();
        dst.close();
        // chemin_recon=repertoire_config_defaut+"/recon.txt";
        cout<<"Utilisation du fichier par defaut : "<<repertoire_config_defaut+chemin_config_manip<<endl;
    }

    ///tester l'existence du fichiers de config : recon.txt dans le chemin indiqué par l'utilisateur
    ifstream fichier_recon(chemin_recon.c_str(), ios::in);
    if(!fichier_recon)
    {
        //  cout<<endl<<"fichier "<<chemin_config_manip<<" inexistant"<<endl;
        //  chemin_config_manip=repertoire_config_defaut+"/config_manip.txt";
        // cout<<"Utilisation du fichier par defaut : "<<chemin_config_manip<<endl;
        cout<<endl<<"fichier "<<chemin_recon<<" inexistant"<<endl;
        ///on copie le fichier defaut dans le repertoire de config...
        ifstream src(repertoire_config_defaut+"recon.txt",ios::binary);
        ofstream dst(chemin_recon,ios::binary);
        dst<<src.rdbuf();
        src.close();
        dst.close();
        // chemin_recon=repertoire_config_defaut+"/recon.txt";
        cout<<"Utilisation du fichier par defaut : "<<repertoire_config_defaut+chemin_recon<<endl;
    }

//    chemin_config_manip_pc_acquis=repertoire_config_pc_acquis+"/config_manip.txt";///chemin vers config_manip.txt  si on est sur un pc acquisition

//    cout<<"config_manip_pc_acquis="<<chemin_config_manip_pc_acquis<<endl;
    vector<string> &r_tab_val_gui_conf=tab_val_gui_conf;///tableau pour Gui_tomo.conf
    vector<string> &r_tab_val_recon=tab_val_recon;///tableau pour recon.txt
    vector<string> &r_tab_val_manip=tab_val_manip;///tableau pour config_manip.txt
    // vector<string> &r_tab_val_manip_pc_acquis=tab_val_manip_pc_acquis;///tableau config_manip.txt si on est sur un pc acquisition

    this->init_tab_val(chemin_recon,r_tab_val_recon);///
    this->init_tab_val(chemin_config_manip,r_tab_val_manip);///
    this->init_tab_val(home+fin_chemin_gui_tomo,r_tab_val_gui_conf);
    // this->init_tab_val(chemin_config_manip_pc_acquis,r_tab_val_manip_pc_acquis);///tableau si on est sur un pc acquisition

     bool b_DARK_THEME=stof(extract_string("DARK_THEME",chemin_config_manip),0);
        ///color bright theme
        wxColor *backgroundColor = new wxColor(200,200,200);
       wxColor *panelColor = new wxColor(235,235,220);
        wxColor *staticTxtColor =  new wxColor(23,42,70);
    ///dark theme
    if(b_DARK_THEME==1){
     *backgroundColor =  wxColor(5,10,20);
     *panelColor = wxColor(23,42,70);
     *staticTxtColor =   wxColor(200,210,250);
    }
 // Create a top-level panel to hold all the contents of the frame
    // Créer le widget wxNotebook
    wxNotebook* notebook = new wxNotebook(this, wxID_ANY);//crer un widget notebook dont le aprent est la class een cours

    wxPanel* panel_tab1 = new wxPanel(notebook, wxID_ANY);//céer un panel dont la fenetre parent est le wxNotebook
    notebook->AddPage(panel_tab1,wxT("Général"));
    panel_tab1->SetBackgroundColour(*backgroundColor);
    // panel_tab1->SetBackgroundColour(wxColor(40, 40, 40));///background colour for the 1st tab



    ///##########TAB 1 ################################################
    wxBoxSizer *sizer_horizontal = new wxBoxSizer(wxHORIZONTAL);
    sizer_horizontal->SetMinSize(270, 100);
    panel_tab1->SetSizer(sizer_horizontal);//affecter l'organiseur (sizer) horizontal à l'onglet  1.
    ///#################### Partie Gauche #############################
    ///Le panneau gauche a pour parent le page onglet 1
    wxPanel *panel_gauche = new wxPanel(panel_tab1,wxID_ANY);

    panel_gauche->SetBackgroundColour(*backgroundColor);///background color for  the left panel
    panel_gauche->SetForegroundColour(*staticTxtColor);//font color for this panel
    //le panneau gauche sera géré par le sizer horizontal global
    sizer_horizontal->Add(panel_gauche, 1, wxALL | wxEXPAND, 2);
    //->Add(parent, etirable horizontablement,bordure=2 pixels )
    //sans l'option d'étirement, la pnneau gauche se limite à la partie affichée
    ///sizer vertical partie gauche
    wxBoxSizer *sizer_vertical0 = new wxBoxSizer(wxVERTICAL);//les widgets seront organisés verticalement dans ce sizer
    //  sizer_horizontal->Add(sizer_vertical0, 1, wxALL | wxEXPAND, 10);
    //wxPanel *zone10 = new wxPanel(panel_gauche,wxID_ANY, wxPoint(1,1),wxSize(200,40),wxSIMPLE_BORDER );
    wxPanel *zone10 = new wxPanel(panel_gauche,wxID_ANY );

    //wxPanel *zone10 = new wxPanel(panel_gauche);
    sizer_vertical0->Add(zone10, 1, wxALL | wxEXPAND, 1);
    zone10->SetBackgroundColour(*panelColor);
    zone10->SetSize(180,10);//ne fonctionne pas en hauteur?
    ///Titre
    // titre_Acquis=new wxStaticText(zone10,0,"Acquisition",wxPoint(2,2),wxSize(100,40));
    titre_Acquis=new wxStaticText(zone10,0,"Acquisition",wxPoint(2,2),wxSize(100,40));
    titre_Acquis->SetFont( wxFont(12, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL) );
    ///--------------------NB HOLO+N0
    wxPanel *zone11 = new wxPanel(panel_gauche);
    // sizer->Add(leftPanel, true, wxEXPAND | wxTOP | wxBOTTOM | wxLEFT, 20);
    sizer_vertical0->Add(zone11, 0, wxALL | wxEXPAND | wxALIGN_CENTER, 1);// Ajout de cette zone au sizer vertical    //
    zone11->SetBackgroundColour(*panelColor);//fond gris clair
    zone11->SetSize(180,40);
    ///NB_HOLO-----
    textNbHolo=new wxStaticText(zone11,1,"Nombre d'hologrammes  : ",wxPoint(2,2),wxSize(98,50));
    //édition du champ
    editNbHolo=new wxTextCtrl(zone11,1,"",wxPoint(120,4),wxSize(50,25));
    editNbHolo->SetToolTip(wxT("Nombre d'hologrammes à acquérir"));
    //récupérer la valeur NB_HOLO dans fichier de config, convertir en string, puis l'afficher dans la case
    int Nb_Holo=stof(extract_string("NB_HOLO",chemin_config_manip));//stof=string to float

    wxString string_NbHolo = wxString::Format(wxT("%i"),Nb_Holo);
    editNbHolo->SetValue(string_NbHolo);

    ///indice milieu montage
    textN0=new wxStaticText(zone11,1,"Indice n0 ",wxPoint(2,40),wxSize(98,50));
    editN0=new wxTextCtrl(zone11,1,"",wxPoint(120,40),wxSize(52,25));
    float edit_N0=stof(extract_string("N0",chemin_config_manip));//stof=string to float
    wxString string_N0 = wxString::Format(wxT("%.3f"),edit_N0); //%.3f=3chiffres apres virgule
    editN0->SetValue(string_N0);
    editN0->SetToolTip(wxT("Indice milieu de montage"));

    ///------------BALAYAGE
    ///--------------------------SCHEMA BALAYAGE
    wxPanel *zone11b = new wxPanel(panel_gauche);
    zone11b->SetSize(200,50);
    sizer_vertical0->Add(zone11b, 1, wxALL | wxEXPAND | wxALIGN_BOTTOM, 1);// Ajout de cette zone au sizer vertical    //
    zone11b->SetBackgroundColour(*panelColor);

    textScanPattern=new wxStaticText(zone11b,1,wxT("&Schéma"),wxPoint(2,10),wxSize(98,50));
    // string str_scanPattern=stof(extract_string("SCAN_PATTERN",chemin_config_manip));//stof=string to float

    ComboBox_Scan = new wxComboBox(zone11b, idComboScan, _T(""), wxPoint(120,10), wxSize(140,26),
                                   0, NULL, wxTE_PROCESS_ENTER | wxTE_PROCESS_TAB);//Constructor, creating and showing a combobox.
    ComboBox_Scan->SetToolTip(wxT("Schéma de balayage"));
    ComboBox_Scan->Append ( wxT("ANNULAR") );///attention, respecter l'ordre alpĥabétique
    ComboBox_Scan->Append ( wxT("ARCHIMEDE") );
    ComboBox_Scan->Append ( wxT("FERMAT") );
    ComboBox_Scan->Append ( wxT("RANDOM_POLAR") );
    ComboBox_Scan->Append ( wxT("ROSACE") );
    ComboBox_Scan->Append ( wxT("STAR") );
    ComboBox_Scan->Append ( wxT("UNIFORM3D") );
    wxArrayString  listeScanPattern(6);
    listeScanPattern=ComboBox_Scan->GetStrings();//récupérer les label de pattern dans un tableau


    //édition du champ cercles concentriques
    int Nb_circles_annular=stof(  extract_string("NB_CIRCLES_ANNULAR", chemin_config_manip, "8")  );//stof=string to float
    editnbCirclesAnnular=new wxTextCtrl(zone11b,1,"",wxPoint(120,40),wxSize(36,24));
    editnbCirclesAnnular->SetToolTip(wxT("Nombre de cercles concentriques si ANNULAR"));
    if(ComboBox_Scan->GetString(ComboBox_Scan->GetSelection())==wxT("ANNULAR")) editnbCirclesAnnular->Enable(true);
    else    editnbCirclesAnnular->Enable(false);

    wxString string_NbCirclesAnnular = wxString::Format(wxT("%i"),Nb_circles_annular);
    editnbCirclesAnnular->SetValue(string_NbCirclesAnnular);

    wxString defaultPattern=extract_string("SCAN_PATTERN",chemin_config_manip,"ROSACE");//lire le pattern du fichier de config, défaut sur Rosace
    cout<<"listeScanPattern.Index(defaultPattern)="<<listeScanPattern.Index(defaultPattern)<<endl;
    ComboBox_Scan->SetSelection(listeScanPattern.Index(defaultPattern));//par defaut, mettre le pattern du fichier de config

    ///--------------------INTERVALLE de BALAYAGE------------------
    wxPanel *zone12 = new wxPanel(panel_gauche);
    sizer_vertical0->Add(zone12, 0, wxALL | wxEXPAND | wxALIGN_CENTER, 1);// Ajout de cette zone au sizer vertical    //
    //  zone12->SetBackgroundColour(*panelColor);
    zone12->SetBackgroundColour(*panelColor);

    titre_Interval_Balayage=new wxStaticText(zone12,-1," Intervalle de Balayage",wxPoint(0,0),wxSize(100,40));
    titre_Interval_Balayage->SetFont(wxFont(10, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD) );

    textVxmin=new wxStaticText(zone12,-1,"Vxm ",wxPoint(10,30),wxSize(30,40));
    editVxmin=new wxTextCtrl(zone12,-1,"",wxPoint(50,28),wxSize(50,25));
    editVxmin->SetToolTip(wxT("tension mini en volt selon l'axe x"));
    float Vxmin=extract_val("VXMIN",chemin_config_manip);
    wxString string_Vxmin = wxString::Format(wxT("%.2f"),Vxmin);
    editVxmin->SetValue(string_Vxmin);

    textVymin=new wxStaticText(zone12,-1,"Vym ",wxPoint(110,30),wxSize(30,40));
    editVymin=new wxTextCtrl(zone12,-1,"",wxPoint(140,28),wxSize(50,25));
    editVymin->SetToolTip(wxT("tension mini en volt selon l'axe y"));
    float Vymin=extract_val("VYMIN",chemin_config_manip);
    wxString string_Vymin = wxString::Format(wxT("%.2f"),Vymin);
    editVymin->SetValue(string_Vymin);

    textVxmax=new wxStaticText(zone12,-1,"VxM ",wxPoint(10,70),wxSize(30,40));
    editVxmax=new wxTextCtrl(zone12,-1,"",wxPoint(50,68),wxSize(50,25));
    editVxmax->SetToolTip(wxT("tension maxi en volt selon l'axe x"));
    float Vxmax=extract_val("VXMAX",chemin_config_manip);
    wxString string_Vxmax = wxString::Format(wxT("%.2f"),Vxmax);
    editVxmax->SetValue(string_Vxmax);

    textVymax=new wxStaticText(zone12,-1,"VyM ",wxPoint(110,70),wxSize(30,40));
    editVymax=new wxTextCtrl(zone12,-1,"",wxPoint(140,68),wxSize(50,25));
    editVymax->SetToolTip(wxT("tension maxi en volt selon l'axe y"));
    float Vymax=extract_val("VYMAX",chemin_config_manip);
    wxString string_Vymax = wxString::Format(wxT("%.2f"),Vymax);
    editVymax->SetValue(string_Vymax);

    textNAcondLim=new wxStaticText(zone12,-1,"NAlim",wxPoint(10,110),wxSize(50,40));
    editNAcondLim=new wxTextCtrl(zone12,-1,"",wxPoint(50,110),wxSize(50,25));
    editNAcondLim->SetToolTip(wxT("coefficient limitant NA balayage, valeur=[0,1]"));
    float NACondLim=extract_val("NA_COND_LIM",chemin_config_manip);
    wxString string_LimNACond = wxString::Format(wxT("%.2f"),NACondLim);
    editNAcondLim->SetValue(string_LimNACond);

    ////-------------------------
    ///--------------------LANCER ACQUISITION------------------
    wxPanel *zone13 = new wxPanel(panel_gauche);
    sizer_vertical0->Add(zone13, 1, wxALL | wxEXPAND | wxALIGN_CENTER, 1);// Ajout de cette zone au sizer vertical    //
    zone13->SetBackgroundColour(*panelColor);
    /// Dialog fichier de données;
    BoutonOpenDir=new wxButton(zone13, idBoutonOpenDir, wxT("&Données"), wxPoint(5,10), wxDefaultSize, 0);
    editDirAcquis=new wxTextCtrl(zone13,-1,chemin_rep_acquis,wxPoint(105,14), wxSize(180,26));
    editDirAcquis->SetToolTip(wxT("Répertoire des acquisitions : "+chemin_rep_acquis));
    // textFicManip=new wxTextCtrl(zone10,-1,"",wxPoint(160,195),wxSize(100,20));
    ///bouton Acquisition
    BoutonManip= new wxButton(zone13, idBoutonManip, wxT("&Acquisition"), wxPoint(194,60), wxDefaultSize, 0);
    if(b_PC_ACQUIS==0) BoutonManip->Enable(false);
    //Bouton image
    // BoutonImage= new wxButton(zone13, idBoutonImage, wxT("&Image"), wxPoint(4,60), wxSize(40,26), 0);
    //BoutonImage->SetToolTip(wxT("Acquisition + Pretraitement + Reconstruction"));

    ///-----------------------------
    ///--------------------------
    panel_gauche->SetSizer(sizer_vertical0);//affecter l'organiseur (sizer) horizontal à l'onglet  1.


///------------------------------------------------------------------------


    ///#################### Partie centrale : PRETRAITEMENT#############################
    //*************Création du wxBoxSizer vertical pour la partie centrale**********
    ///Le panneau central a pour parent le page onglet 1
    wxPanel *panel_centre = new wxPanel(panel_tab1,wxID_ANY);
    //le panneau gauche sera géré par le sizer horizontal global
    sizer_horizontal->Add(panel_centre, 1, wxALL | wxEXPAND, 2);
    wxBoxSizer *sizer_vertical1 = new wxBoxSizer(wxVERTICAL);
    panel_centre->SetBackgroundColour(*backgroundColor);
    panel_centre->SetForegroundColour(*staticTxtColor);
    //-------------zone du haut -------------------------------------------------------------------------------------
    wxPanel *zone20 = new wxPanel(panel_centre);
    zone20->SetBackgroundColour(*panelColor);                                  //
    sizer_vertical1->Add(zone20, 0, wxALL | wxEXPAND | wxALIGN_CENTER, 1);// Ajout de cette zone au sizer vertical    //
    titre_Pretraitement=new wxStaticText(zone20,-1,wxT(" Prétraitement"),wxPoint(0,0),wxSize(100,40));               //
    titre_Pretraitement->SetFont(wxFont(12, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL) );      //
    //-------------zone du milieu ------------------------------------------------------------------------------------
    wxPanel *zone21 = new wxPanel(panel_centre, wxID_ANY, wxDefaultPosition, wxSize(-1, 160));
    zone21->SetBackgroundColour(*panelColor);

    // Ajout de cette zone au sizer vertical
    sizer_vertical1->Add(zone21, 0, wxALL | wxEXPAND, 1);

    m_born = new wxCheckBox(zone21, id_boolBorn, wxT("Utiliser Born"), wxPoint(20, 11));
    m_born->SetToolTip(wxT("0=Rytov"));
    m_born->SetValue(stof(extract_string("BORN",chemin_recon,"0")));

    m_volkov = new wxCheckBox(zone21, idBoolVolkov, wxT("Utiliser Volkov"), wxPoint(20, 41));
    m_volkov->SetToolTip(wxT("0=Herraez"));
    m_volkov->SetValue(stof(extract_string("VOLKOV",chemin_recon,"1")));
    //m_born->SetValue(true);
    m_Deroul = new wxCheckBox(zone21, idBoolDeroul, wxT("Déroulement"), wxPoint(20, 71));
    m_Deroul->SetValue(stof(extract_string("DEROUL",chemin_recon,"1")));
    m_Aber = new wxCheckBox(zone21, idBoolAber, wxT("Correction aberration"), wxPoint(20, 101));
    m_Aber->SetValue(stof(extract_string("C_ABER",chemin_recon)));

    m_AmpliRef = new wxCheckBox(zone21, idBoolAmpliRef, wxT("Division ampli_ref"), wxPoint(20, 131));
    m_AmpliRef->SetValue(stof(extract_string("AMPLI_REF",chemin_recon,"1")));

    //-------------zone du bas---------------------------------------------------------------------------------------//
    wxPanel *zone22 = new wxPanel(panel_centre);
    zone22->SetBackgroundColour(*panelColor);
    // Ajout de cette zone au sizer vertical
    sizer_vertical1->Add(zone22, 2, wxALL | wxEXPAND, 1);
    /// Dialog fichier;
    BoutonOpenMask=new wxButton(zone22, idBoutonOpenMask, wxT("&Masque"), wxPoint(5,10), wxSize(70,30), 0);
    editFicMask=new wxTextCtrl(zone22,-1,chemin_rep_acquis+"/Mask.png",wxPoint(85,6), wxSize(200,40));
    editFicMask->SetToolTip(wxT("Masque pour correction d'aberration"));

    BoutonRecons= new wxButton(zone22, idBoutonTraitement, wxT("&Prétraitement"), wxPoint(180,160), wxDefaultSize, 0);


    panel_centre->SetSizer(sizer_vertical1);//affecter l'organiseur (sizer) horizontal à l'onglet  1.


    ///--------------------REconstruction---------------------------------------------
    ///#############################Partie droite : Reconstruction#################################

    ///Le panneau droit a pour parent le page onglet 1
    wxPanel *panel_droit = new wxPanel(panel_tab1,wxID_ANY);
    panel_droit->SetBackgroundColour(*backgroundColor);
    panel_droit->SetForegroundColour(*staticTxtColor);
    sizer_horizontal->Add(panel_droit, 1, wxALL | wxEXPAND, 2); //le panneau gauche sera géré par le sizer horizontal global

    wxBoxSizer *sizer_vertical2 = new wxBoxSizer(wxVERTICAL);

    wxPanel *zone30= new wxPanel(panel_droit);
    zone30->SetBackgroundColour(*panelColor);
    // Ajout de cette zone au sizer vertical
    sizer_vertical2->Add(zone30, 0, wxALL | wxEXPAND | wxALIGN_CENTER, 1);
    titre_Recons=new wxStaticText(zone30,-1,wxT(" Reconstruction"),wxPoint(0,0),wxSize(100,40));
    titre_Recons->SetFont(wxFont(12, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL) );
    //-------------zone Hors axe -------------------------------------------------------------------------------------
    wxPanel *zone31 = new wxPanel(panel_droit);
    zone31->SetBackgroundColour(*panelColor);
    sizer_vertical2->Add(zone31, 0, wxALL | wxEXPAND | wxALIGN_CENTER, 1);
    titre_HorsAxe=new wxStaticText(zone31,-1,wxT(" Hors-axe"),wxPoint(0,0),wxSize(100,40));
    titre_HorsAxe->SetFont(wxFont(10, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD) );

    //CX-----
    textCX=new wxStaticText(zone31,-1,"fx0 ",wxPoint(10,30),wxSize(26,40));
    editCX=new wxTextCtrl(zone31,-1,"",wxPoint(36,25),wxSize(50,25));
    editCX->SetToolTip(wxT("Centre X du spectre"));
    int CX=stof(extract_string("CIRCLE_CX",chemin_config_manip));
    wxString string_CX = wxString::Format(wxT("%i"),CX);
    editCX->SetValue(string_CX);
    //CY------
    textCY=new wxStaticText(zone31,-1,"fy0 ",wxPoint(91,30),wxSize(28,40));
    editCY=new wxTextCtrl(zone31,-1,"",wxPoint(115,25),wxSize(50,25));
    int CY=stof(extract_string("CIRCLE_CY",chemin_config_manip));
    wxString string_CY = wxString::Format(wxT("%i"),CY);
    editCY->SetValue(string_CY);
    editCY->SetToolTip(wxT("Centre Y du spectre"));
    //NXMAX------
    textNXMAX=new wxStaticText(zone31,-1,"NXMAX ",wxPoint(168,30),wxSize(60,40));
    editNXMAX=new wxTextCtrl(zone31,-1,"",wxPoint(225,25),wxSize(50,25));
    int NXMAX=stof(extract_string("NXMAX",chemin_config_manip));
    wxString string_NXMAX = wxString::Format(wxT("%i"),NXMAX);
    editNXMAX->SetValue(string_NXMAX);
    editNXMAX->SetToolTip(wxT("Demi largeur du spectre"));

    //DIMFINALE
    textDimFinal=new wxStaticText(zone31,-1,"DimFinal ",wxPoint(166,65),wxSize(60,40));
    editDimFinal=new wxTextCtrl(zone31,-1,"",wxPoint(225,60),wxSize(50,25));
    int DimFinal=stof(extract_string("DIM_FINAL",chemin_recon,"512"));
    wxString string_DimFinal = wxString::Format(wxT("%i"),DimFinal);
    editDimFinal->SetValue(string_DimFinal);
    editDimFinal->SetToolTip(wxT("Dimension du volume final"));

    //-------------zone du bas -------------------------------------------------------------------------------------
    wxPanel *zone32 = new wxPanel(panel_droit);
    zone32->SetBackgroundColour(*panelColor);
    // Ajout de cette zone au sizer vertical
    sizer_vertical2->Add(zone32, 1, wxALL | wxEXPAND | wxALIGN_CENTER, 1);


    m_export_OTF = new wxCheckBox(zone32, idBoolExportOTF, wxT("Exporter OTF"), wxPoint(5, 50));
    m_export_OTF->SetToolTip(wxT("Exporter OTF pour reconstruction GP"));
    m_export_OTF->SetValue(stof(extract_string("EXPORT_OTF",chemin_recon)));

    BoutonOpenDirResult=new wxButton(zone32, idBoutonOpenDirResult, wxT("& Résultats "), wxPoint(5,10), wxDefaultSize, 0);
    editDirResultAcquis=new wxTextCtrl(zone32,-1,chemin_result,wxPoint(105,12),wxSize(170,30));
    editDirResultAcquis->SetToolTip(wxT("Répertoire résultats"));
    // textFicManip=new wxTextCtrl(zone10,-1,"",wxPoint(160,195),wxSize(100,20));
    BoutonRecons= new wxButton(zone32, idBoutonRecons, wxT("&Reconstruction"), wxPoint(174,214), wxDefaultSize, 0);

    // BoutonOpenDir=new wxButton(zone13, idBoutonOpenDir, wxT("&Données"), wxPoint(5,10), wxDefaultSize, 0);
    // editDirAcquis=new wxTextCtrl(zone13,-1,chemin_rep_acquis,wxPoint(105,14), wxSize(175,26));
    //  editDirAcquis->SetToolTip(wxT("Répertoire des acquisitions"));

    panel_droit->SetSizer(sizer_vertical2);//affecter l'organiseur (sizer) horizontal à l'onglet  1.


    ///##########TAB 2 ################################################
    wxPanel* panel_tab2 = new wxPanel(notebook, wxID_ANY,wxPoint(1,1),wxSize(100,40));//céer un panel dont la fenetre parent est le wxNotebook
    notebook->AddPage(panel_tab2,wxT("Avancé"));
    panel_tab2->SetBackgroundColour(wxColor(10, 25, 50));
    panel_tab2->SetForegroundColour(wxColor(200, 210, 250));
    wxBoxSizer *sizer_horizontalB = new wxBoxSizer(wxHORIZONTAL);
    sizer_horizontalB->SetMinSize(200, 100);
    panel_tab2->SetSizer(sizer_horizontalB);//affecter l'organiseur (sizer) horizontal à l'onglet  1.
    ///#################### Partie Gauche #############################
    ///Le panneau gauche a pour parent le page onglet 2
    wxPanel *panel_gauche_tabB = new wxPanel(panel_tab2,wxID_ANY,wxPoint(1,1),wxSize(200,40));
    panel_gauche_tabB->SetBackgroundColour(wxColor(10, 25, 50));
    //le panneau gauche sera géré par le sizer horizontal global
    sizer_horizontalB->Add(panel_gauche_tabB, 0, wxALL | wxEXPAND, 2);
    //->Add(parent, etirable horizontablement,bordure=2 pixels )
    //sans l'option d'étirement, la pnneau gauche se limite à la partie affichée


    ///sizer vertical partie gauche
    wxBoxSizer *sizer_vertical_tabB_0 = new wxBoxSizer(wxVERTICAL);
    //  sizer_horizontal->Add(sizer_vertical0, 1, wxALL | wxEXPAND, 10);
    wxPanel *zoneB10 = new wxPanel(panel_gauche_tabB,wxID_ANY, wxPoint(1,1),wxSize(200,30) );
    sizer_vertical_tabB_0->Add(zoneB10, 0, wxALL |wxALIGN_LEFT, 1);//zoneB10, non étirable verticalement
    zoneB10->SetBackgroundColour(*panelColor);
    //zoneB10->SetBor
    ///Titre
    titre_GPS=new wxStaticText(zoneB10,0,"Traitement GPS",wxPoint(2,2),wxSize(100,20));
    titre_GPS->SetFont( wxFont(12, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL) );

    ///--------------------Intervalle indice------------------
    wxPanel *zoneB11 = new wxPanel(panel_gauche_tabB, wxID_ANY, wxPoint(1,1),wxSize(200,40));
    sizer_vertical_tabB_0->Add(zoneB11, 1, wxALL |wxALIGN_LEFT, 1);// Ajout de cette zone au sizer vertical    //
    zoneB11->SetBackgroundColour(*panelColor);

    titre_interval_indice=new wxStaticText(zoneB11,-1,"Gamme d'indice",wxPoint(0,0),wxSize(100,40));
    titre_interval_indice->SetFont(wxFont(10, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD) );

    textDnMin=new wxStaticText(zoneB11,-1,wxT("Δn Min "),wxPoint(10,30),wxSize(80,40));
    textDnMin->SetFont(wxFont(10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD) );
    editDnMin=new wxTextCtrl(zoneB11,-1,"",wxPoint(70,26),wxSize(50,25));
    editDnMin->SetToolTip(wxT("Delta n minimum"));
    float DnMin=extract_val("DELTA_NMIN",chemin_recon);
    wxString string_DnMin = wxString::Format(wxT("%.2f"),DnMin);
    editDnMin->SetValue(string_DnMin);


    textDnMax=new wxStaticText(zoneB11,-1,wxT("Δn Max "),wxPoint(10,60),wxSize(80,40));
    textDnMax->SetFont(wxFont(10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD ) );
    editDnMax=new wxTextCtrl(zoneB11,-1,"",wxPoint(70,58),wxSize(50,25));
    editDnMax->SetToolTip(wxT("Delta n maximum"));
    float DnMax=extract_val("DELTA_NMAX",chemin_recon);
    wxString string_DnMax = wxString::Format(wxT("%.2f"),DnMax);
    editDnMax->SetValue(string_DnMax);

    ///--------------------Intervalle coefficient extinction ------------------
    wxPanel *zoneB12 = new wxPanel(panel_gauche_tabB, wxID_ANY, wxPoint(1,1),wxSize(200,40));
    sizer_vertical_tabB_0->Add(zoneB12, 1, wxALL |wxALIGN_LEFT, 1);// Ajout de cette zone au sizer vertical    //
    zoneB12->SetBackgroundColour(*panelColor);

    titre_interval_kappa=new wxStaticText(zoneB12,-1,"Coefficient d'extinction",wxPoint(0,0),wxSize(100,40));
    titre_interval_kappa->SetFont(wxFont(10, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD) );

    textKappaMin=new wxStaticText(zoneB12,-1,wxT("κ Min "),wxPoint(10,30),wxSize(80,40));
    textKappaMin->SetFont(wxFont(10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD) );
    editKappaMin=new wxTextCtrl(zoneB12,-1,"",wxPoint(70,26),wxSize(80,25));
    editKappaMin->SetToolTip(wxT("Kappa minimum"));
    float KappaMin=extract_val("KAPPA_MIN",chemin_recon);
    wxString string_KappaMin = wxString::Format(wxT("%.2f"),KappaMin);
    editKappaMin->SetValue(string_KappaMin);


    textKappaMax=new wxStaticText(zoneB12,-1,wxT("κ Max "),wxPoint(10,60),wxSize(80,40));
    textKappaMax->SetFont(wxFont(10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD ) );
    editKappaMax=new wxTextCtrl(zoneB12,-1,"",wxPoint(70,58),wxSize(80,25));
    editKappaMax->SetToolTip(wxT("Coef d'extinction maximum"));
    float KappaMax=extract_val("KAPPA_MAX",chemin_recon);
    wxString string_KappaMax = wxString::Format(wxT("%.2f"), KappaMax);
    editKappaMax->SetValue(string_KappaMax);

    ///----------NB ITERATION GPS-----------------
    wxPanel *zoneB13 = new wxPanel(panel_gauche_tabB, wxID_ANY, wxPoint(1,1),wxSize(200,40));
    sizer_vertical_tabB_0->Add(zoneB13, 1, wxALL |wxALIGN_LEFT, 1);// Ajout de cette zone au sizer vertical    //
    zoneB13->SetBackgroundColour(*panelColor);

    titre_iterationGPS=new wxStaticText(zoneB13,-1,wxT("Nombre d'itérations"),wxPoint(4,6),wxSize(200,40));
    titre_iterationGPS->SetFont(wxFont(10, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD) );

    textIteration=new wxStaticText(zoneB13,-1,wxT("NB_ITER"),wxPoint(6,39),wxSize(150,40));
    textIteration->SetFont(wxFont(10, wxFONTFAMILY_MODERN, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD) );

    editIteration=new wxTextCtrl(zoneB13,-1,"",wxPoint(70,35),wxSize(80,25));
    editIteration->SetToolTip(wxT("Nombre d'itération GPS"));
    unsigned int NbIteration=extract_val("NB_ITER_GPS",chemin_recon);
    wxString string_Iteration = wxString::Format(wxT("%.i"),NbIteration);
    editIteration->SetValue(string_Iteration);

    BoutonGPS= new wxButton(zoneB13, idBoutonGPS, wxT("&GPS"), wxPoint(140,65), wxSize(50,28), 0);

    panel_gauche_tabB->SetSizer(sizer_vertical_tabB_0);
}
///########################################functions linked to widgets ###############"
Gui_TomoFrame::~Gui_TomoFrame()
{

}

void Gui_TomoFrame::OnClose(wxCloseEvent &event)
{
    Destroy();
}

void Gui_TomoFrame::OnQuit(wxCommandEvent &event)
{
    Destroy();
}
void Gui_TomoFrame::OnAbout(wxCommandEvent &event)
{
    wxString msg = ("Interface graphique \npour microscope tomographique \nV.0.1");
    wxMessageBox(msg, _("Bienvenue..."));

}
void Gui_TomoFrame::OnBoutonManip(wxCommandEvent &event)///Acquisition bouton
{
    // float Rx=abs(Vxmax-Vxmin);
//    float offX=Vxmax-Rx;
    //  float Ry=abs(Vymax-Vymin);
    //float offY=Vymax-Ry;

    // wxString param1 =editDirAcquis->GetValue();
    //  wxString string_NbHolo=editNbHolo->GetValue();
    // wxString cmd_final=wxT("manipTomo5 -ni ")+editNbHolo->GetValue()+wxT(" -voffset ")+wxString::Format(wxT("%.2f"), offX)+" "+
    // wxString::Format(wxT("%.2f"),offY)+wxT(" -vfleur ")+wxString::Format(wxT("%.2f"),Rx)+" "+wxString::Format(wxT("%.2f"),Ry);
    sav_all();
    wxExecute(wxT("tomo_manip"));

//    cout<<cmd_final<<endl;
    // wxExecute(cmd_final);
}
void Gui_TomoFrame::OnBoutonImage(wxCommandEvent &event)///image complète
{
    // float Rx=abs(Vxmax-Vxmin);
//    float offX=Vxmax-Rx;
    //  float Ry=abs(Vymax-Vymin);
    //float offY=Vymax-Ry;

    // wxString param1 =editDirAcquis->GetValue();
    //  wxString string_NbHolo=editNbHolo->GetValue();
    // wxString cmd_final=wxT("manipTomo5 -ni ")+editNbHolo->GetValue()+wxT(" -voffset ")+wxString::Format(wxT("%.2f"), offX)+" "+
    // wxString::Format(wxT("%.2f"),offY)+wxT(" -vfleur ")+wxString::Format(wxT("%.2f"),Rx)+" "+wxString::Format(wxT("%.2f"),Ry);
    sav_all();
    wxExecute(wxT("tomo_manip && tomo_pretraitement && tomo_reconsruction"));

//    cout<<cmd_final<<endl;
    // wxExecute(cmd_final);
}


void Gui_TomoFrame::OnBoutonTraitement(wxCommandEvent &event)
{
    sav_all();
    wxExecute(wxT("tomo_pretraitement"));
}
void Gui_TomoFrame::OnBoutonGPS(wxCommandEvent &event)
{
    sav_all();
    wxExecute(wxT("tomo_GPS"));
}
void Gui_TomoFrame::OnBoutonRecons(wxCommandEvent &event)
{
    // wxExecute(wxT("Tomo_pretraitement -i ")+editDirResultAcquis->GetValue());
    sav_all();
    wxExecute(wxT("tomo_reconstruction"));
}

void Gui_TomoFrame::OnBoutonBorn(wxCommandEvent& WXUNUSED(event))
{
    if (m_born->GetValue()==true)
    {
        modif_tab_val("BORN","1",tab_val_recon);
    }
    else
    {
        modif_tab_val("BORN","0",tab_val_recon);
        modif_tab_val("DEROUL", "1",tab_val_recon);
        m_Deroul->SetValue("1");
    }
}

void Gui_TomoFrame::OnBoutonExportOTF(wxCommandEvent& WXUNUSED(event))
{
    if (m_export_OTF->GetValue()==true)
    {
        modif_tab_val("EXPORT_OTF","1",tab_val_recon);
    }
    else
    {
        modif_tab_val("EXPORT_OTF","0",tab_val_recon);
    }
}

void Gui_TomoFrame::OnBoutonDeroul(wxCommandEvent& WXUNUSED(event))
{
    if(m_Deroul->GetValue()==true)
    {
        modif_tab_val("DEROUL", "1",tab_val_recon);
    }
    else
    {
        modif_tab_val("DEROUL", "0",tab_val_recon);
    }
}

///modifier aussi la sauvegarde, sinon fichier non modifié
void Gui_TomoFrame::OnBoutonVolkov(wxCommandEvent& WXUNUSED(event))
{
    if(m_volkov->GetValue()==true)
    {
        modif_tab_val("VOLKOV", "1",tab_val_recon);
    }
    else
    {
        modif_tab_val("VOLKOV", "0",tab_val_recon);
    }
}

void Gui_TomoFrame::OnBoutonAber(wxCommandEvent& WXUNUSED(event))
{
    if (m_Aber->GetValue()==true)
    {
        modif_tab_val("C_ABER", "1",tab_val_recon);
    }
    else
    {
        modif_tab_val("C_ABER", "0",tab_val_recon);
    }
}

void Gui_TomoFrame::OnBoutonAmpliRef(wxCommandEvent& WXUNUSED(event))
{
    if (m_AmpliRef->GetValue()==true)
    {
        modif_tab_val("AMPLI_REF", "1",tab_val_recon);
    }
    else
    {
        modif_tab_val("AMPLI_REF", "0",tab_val_recon);
    }
}

void Gui_TomoFrame::OnBoutonOpenDir(wxCommandEvent& WXUNUSED(event))
{
    // wxString defaultPath = chemin_rep_acquis;//wxT("/ramdisk/Acquis/");
    wxDirDialog* OpenDirDial = new wxDirDialog(this, _("Choisir un répertoire"),editDirAcquis->GetValue());

    if ( OpenDirDial->ShowModal() == wxID_OK )
    {
        wxString path = OpenDirDial->GetPath();
        editDirAcquis->SetValue(path);
    }
    //récupérer le nouveau chemin
    string dir_config(editDirAcquis->GetValue().ToStdString());//ToStdString =wxwidget 3.0. Non compatible wxWidget 2.8

    chemin_config_manip=dir_config+"/config_manip.txt";
    chemin_recon=dir_config+"/recon.txt";
    refreshValue(chemin_config_manip, chemin_recon);

}


void Gui_TomoFrame::OnBoutonOpenDirResult(wxCommandEvent& WXUNUSED(event))
{
    wxString defaultPath = chemin_result;
    wxDirDialog* OpenDirDial = new wxDirDialog(this, _("Choisissez un répertoire"),defaultPath);
    cout<<"dans OpenDirResult"<<endl;
    if ( OpenDirDial->ShowModal() == wxID_OK )
    {
        wxString path = OpenDirDial->GetPath();
        editDirResultAcquis->SetValue(path);
    }
}

void Gui_TomoFrame::OnBoutonOpenMask(wxCommandEvent& WXUNUSED(event))
{
    wxString defaultPath = chemin_rep_acquis;
    wxFileDialog* OpenFileDial = new wxFileDialog(this, _("Choisissez un fichier"),defaultPath);
    if ( OpenFileDial->ShowModal() == wxID_OK )
    {
        wxString path = OpenFileDial->GetPath();
        editFicMask->SetValue(path);
    }
}

void Gui_TomoFrame::OnBoutonSAV(wxCommandEvent& WXUNUSED(event))
{
    sav_all();
}




///######################################## General functions ###############"
//rafraichit toutes les valeurs de l'interface (avec setvalue) avec les valeurs du fichier "chemin_fic", "chemin_recon"
void Gui_TomoFrame::refreshValue(string chemin_fic, string chemin_recon)
{
    cout<<endl;
    cout<<"####-- Maj des données avec le fichier : "<<chemin_fic<<" --####"<<endl;
    ///-------------MAJ paramètres manip----------------------------------
    editNbHolo->SetValue(wxString::Format(wxT("%i"),stoi(extract_string("NB_HOLO",chemin_fic,"400"))));
    //editN0->SetValue(wxString::Format(wxT("%.2f"),stoi(extract_string("N0",chemin_fic))));
    editN0->SetValue(wxString::Format( wxT("%.3f"),extract_val("N0",chemin_fic) ));
    editVymin->SetValue(wxString::Format( wxT("%.2f"),extract_val("VYMIN",chemin_fic) ));
    editVxmin->SetValue(wxString::Format(wxT("%.2f"),extract_val("VXMIN",chemin_fic) ));
    editVxmax->SetValue(wxString::Format(wxT("%.2f"), extract_val("VXMAX",chemin_fic)));
    editVymax->SetValue(wxString::Format(wxT("%.2f"),extract_val("VYMAX",chemin_fic)));
    editNAcondLim->SetValue(wxString::Format(wxT("%.2f"),  extract_val("NA_COND_LIM",  chemin_fic)));

    wxArrayString  listeScanPattern(6);
    listeScanPattern=ComboBox_Scan->GetStrings();//récupérer les label de pattern dans un tableau de string
    // wxString defaultPattern=extract_string("SCAN_PATTERN",chemin_fic);//lire le pattern du fichier de config
    ComboBox_Scan->SetSelection(listeScanPattern.Index(extract_string("SCAN_PATTERN",chemin_fic,"ROSACE")));//

    editNXMAX->SetValue(wxString::Format(wxT("%i"),stoi(extract_string("NXMAX",chemin_fic,"0"))));
    editCX->SetValue(wxString::Format(wxT("%i"),stoi(extract_string("CIRCLE_CX",chemin_fic,"0"))));
    editCY->SetValue(wxString::Format(wxT("%i"), stoi(extract_string("CIRCLE_CY",chemin_fic,"0"))));
    editnbCirclesAnnular->SetValue(wxString::Format( wxT("%i"), stoi(extract_string("NB_CIRCLES_ANNULAR",chemin_fic,"8")) ));

    ///-------------MAJ paramètres recon----------------------------------
    m_born->SetValue(stof(extract_string("BORN",chemin_recon,"0")));
    m_volkov->SetValue(stof(extract_string("VOLKOV",chemin_recon,"1")));
    m_Deroul->SetValue(stof(extract_string("DEROUL",chemin_recon,"1")));
    m_Aber->SetValue(stof(extract_string("C_ABER",chemin_recon)));
    m_AmpliRef->SetValue(stof(extract_string("AMPLI_REF",chemin_recon,"1")));
    editDimFinal->SetValue(wxString::Format(wxT("%i"),stoi(extract_string("DIM_FINAL",chemin_recon,"512"))));
    m_export_OTF->SetValue(stof(extract_string("EXPORT_OTF",chemin_recon)));
    editDnMin->SetValue(wxString::Format(wxT("%.2f"), extract_val("DELTA_NMIN",chemin_recon)));
    editDnMax->SetValue(wxString::Format(wxT("%.2f"),extract_val("DELTA_NMAX",chemin_recon)));
    editKappaMin->SetValue(wxString::Format(wxT("%.2f"),extract_val("KAPPA_MIN", chemin_recon)));
    editKappaMax->SetValue(wxString::Format(wxT("%.2f"), extract_val("KAPPA_MAX",chemin_recon)));
    editIteration->SetValue(wxString::Format(wxT("%.i"),(int)extract_val("NB_ITER_GPS",chemin_recon)));
    sav_all();
}


//extract a string from the file whose path is  "chemin_fic",
string Gui_TomoFrame::extract_string(string token,  string chemin_fic)
{
    ifstream fichier(chemin_fic.c_str(), ios::in);  // on ouvre en lecture
    string ligne,motcle,valeurMot,separ=" ";
    string valeur;
    vector<std::string> tokens;//tableau contenant les lignes != token !!

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
        cerr << "Impossible d'ouvrir le fichier : "<< chemin_fic<< endl;

    int nb_tok=tokens.size();//length of thetable
    for(int cpt=0; cpt<nb_tok; cpt++) //parcourir le tableau
    {
        ligne=tokens[cpt];//récupérer une ligne
        if(ligne!="") //si elle n'est aps vide
        {
            int pos_separ=ligne.find(separ);//trouver la position du séparateur (délime mot-clé et valeur)
            int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin) : le mot clé est entre la position 0 et la position du séparateur
            if(motcle==token)
            {
                valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));
                ///remove spaces possibly inserted at the beginning (position=0) of the keyword.
                pos_separ=valeurMot.find(separ);//find space inside the keyword
                // char current_character=valeurMot[pos_separ];
                while(valeurMot[pos_separ]==separ) //is the current character equal to "space" ? then remove it
                {
                    valeurMot=valeurMot.substr(pos_separ+1,valeurMot.length()-1);//on raccourci la chaine de caractère de l'espace
                    pos_separ++;
                }
                cout<<motcle<<"="<<valeurMot<<endl;
                valeur=valeurMot.c_str();
            }
        }
    }
    if(valeur.empty())
    {
        cout<<"mot_clé "<<token<<" inexistant dans le fichier "<<chemin_fic<<endl;
        exit(0);
        //valeur="";
    }
    fichier.close();
    return valeur;
}




//extract a string from the file whose path is  "chemin_fic", but in addition can specify a default value
string Gui_TomoFrame::extract_string(string token,  string chemin_fic,string defaut)
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
        cerr << "Impossible d'ouvrir le fichier : "<< chemin_fic<< endl;
    ///ligne trouvée, on passe au mot clé
    int nb_tok=tokens.size();
    for(int cpt=0; cpt<nb_tok; cpt++)
    {
        ligne=tokens[cpt];
        if(ligne!="")
        {
            int pos_separ=ligne.find(separ);//find cherche la string "separ" et retourne sa position
            int long_separ=separ.length();//retourne la longueur du séparateur, ici 1 car séparateur=espace
            //retourne le mot-clé
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,longuieur) : copie une string depuis 'pos_debut' sur une longueur

            if(motcle==token)
            {
                valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));

                ///remove spaces possibly inserted by mistake at the beginning (position=0) of the keyword.
                pos_separ=valeurMot.find(separ);//find space inside the keyword
                while(valeurMot[pos_separ]==separ) //is the current character equal to "space" ? then remove it
                {
                    valeurMot=valeurMot.substr(pos_separ+1,valeurMot.length()-1);//on raccourci la chaine de caractère de l'espace
                    pos_separ++;
                }
                cout<<motcle<<"="<<valeurMot<<endl;
                valeur=valeurMot.c_str();
            }
        }
    }
    if(valeur.empty())
    {
        cout<<"mot_clé  bidule "<<token<<" inexistant dans le fichier "<<chemin_fic<<endl;
        cout<<"Utlisation valeur par défaut: "<<defaut<<endl;
        valeur=defaut;
    }
    fichier.close();
    return valeur;
}
///charger en mémoire le tableau de mot-clés à partir d'un fichier de config situé à chemin_fic.
void Gui_TomoFrame::init_tab_val(string chemin_fic, vector<string> &tab_val)
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
void Gui_TomoFrame::modif_tab_val(string token,string string_valeur_token,vector<string> &tab_val)
{
    std::vector<std::string> fichier_new;
    string ligne,motcle,string_val_token,nlle_ligne,separ=" ";

    for(size_t cpt=0; cpt<tab_val.size(); cpt++)
    {
        ligne=tab_val[cpt];

        if(ligne!='#')//si pas ligne de commentaire (premier caractère ligne!=#)
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

void Gui_TomoFrame::sav_val(string chemin_fic,vector<string> &tab_val)
{
    ///on sauvegarde le fichier, au cas où...
    ifstream src(chemin_fic.c_str(),ios::binary);
    ofstream dst(chemin_fic+"_SAV",ios::binary);
    dst<<src.rdbuf();
    src.close();
    dst.close();
    ofstream flux_ecriture(chemin_fic.c_str(), ios::out);  // ouverture en écriture
    for(size_t cpt=0; cpt<tab_val.size(); cpt++)
    {
        flux_ecriture<<tab_val[cpt]<<endl;//ecriture dans le fichier
    }
    flux_ecriture.close();  //  fermer
}

///-------les champ edit sont modifiés dans tab_val ici, les booleens sont déjà modifiés dans les fonctions de boutons
//sauvegarde des paramètres dans le dossier des données, tous pc (acquis et recon)
/*void Gui_TomoFrame::sav_all_pc_recon()
{
    //tableau de config de l'interface graphique (chemins)
    modif_tab_val("CHEMIN_RESULT",editDirResultAcquis->GetValue().ToStdString(),tab_val_gui_conf);
    modif_tab_val("CHEMIN_CONFIG",editDirAcquis->GetValue().ToStdString(),tab_val_gui_conf);
    modif_tab_val("CHEMIN_ACQUIS",editDirAcquis->GetValue().ToStdString(),tab_val_gui_conf);

    modif_tab_val("CHEMIN_MASK",editFicMask->GetValue().ToStdString(),tab_val_gui_conf);
 //   modif_tab_val("NB_HOLO",editNbHolo->GetValue().ToStdString(),tab_val_manip);
    //tableau de config  reconstruction
    modif_tab_val("FINAL_ANGLE",editNbHolo->GetValue().ToStdString(),tab_val_recon);

    modif_tab_val("DIM_FINAL",editDimFinal->GetValue().ToStdString(),tab_val_recon);

    //tableau de config  manip
   // modif_tab_val("NXMAX",editNXMAX->GetValue().ToStdString(),tab_val_manip);
   // modif_tab_val("VXMIN",editVxmin->GetValue().ToStdString(),tab_val_manip);
  //  modif_tab_val("VYMIN",editVymin->GetValue().ToStdString(),tab_val_manip);
  //  modif_tab_val("VXMAX",editVxmax->GetValue().ToStdString(),tab_val_manip);
 //   modif_tab_val("VYMAX",editVymax->GetValue().ToStdString(),tab_val_manip);
    //tableau de config  recons
    modif_tab_val("DELTA_NMAX",editDnMax->GetValue().ToStdString(),tab_val_recon);
    modif_tab_val("DELTA_NMIN",editDnMin->GetValue().ToStdString(),tab_val_recon);
    modif_tab_val("KAPPA_MAX",editKappaMax->GetValue().ToStdString(),tab_val_recon);
    modif_tab_val("KAPPA_MIN",editKappaMin->GetValue().ToStdString(),tab_val_recon);
    modif_tab_val("NB_ITER_GPS",editIteration->GetValue().ToStdString(),tab_val_recon);

   // modif_tab_val("CIRCLE_CX",editCX->GetValue().ToStdString(),tab_val_manip);
  //  modif_tab_val("CIRCLE_CY",editCY->GetValue().ToStdString(),tab_val_manip);


    sav_val(chemin_recon,tab_val_recon);
   // sav_val(chemin_config_manip,tab_val_manip);
    sav_val(chemin_config_GUI,tab_val_gui_conf);
}*/
void Gui_TomoFrame::sav_all()
{
    //tableau de config de l'interface graphique (chemins)
    modif_tab_val("CHEMIN_RESULT",editDirResultAcquis->GetValue().ToStdString(),tab_val_gui_conf);
    modif_tab_val("CHEMIN_CONFIG",editDirAcquis->GetValue().ToStdString(),tab_val_gui_conf);
    modif_tab_val("CHEMIN_ACQUIS",editDirAcquis->GetValue().ToStdString(),tab_val_gui_conf);

    modif_tab_val("CHEMIN_MASK",editFicMask->GetValue().ToStdString(),tab_val_gui_conf);
    modif_tab_val("N0",editN0->GetValue().ToStdString(),tab_val_manip);
    modif_tab_val("NB_HOLO",editNbHolo->GetValue().ToStdString(),tab_val_manip);
    //tableau de config  reconstruction
    modif_tab_val("FINAL_ANGLE",editNbHolo->GetValue().ToStdString(),tab_val_recon);

    modif_tab_val("DIM_FINAL",editDimFinal->GetValue().ToStdString(),tab_val_recon);

    //tableau de config  manip
    modif_tab_val("SCAN_PATTERN",ComboBox_Scan->GetString(ComboBox_Scan->GetSelection()).ToStdString(),tab_val_manip);
    modif_tab_val("NB_CIRCLES_ANNULAR",editnbCirclesAnnular->GetValue().ToStdString(),tab_val_manip);
    modif_tab_val("NXMAX",editNXMAX->GetValue().ToStdString(),tab_val_manip);
    modif_tab_val("VXMIN",editVxmin->GetValue().ToStdString(),tab_val_manip);
    modif_tab_val("VYMIN",editVymin->GetValue().ToStdString(),tab_val_manip);
    modif_tab_val("VXMAX",editVxmax->GetValue().ToStdString(),tab_val_manip);
    modif_tab_val("VYMAX",editVymax->GetValue().ToStdString(),tab_val_manip);
    modif_tab_val("NA_COND_LIM",editNAcondLim->GetValue().ToStdString(),tab_val_manip);
    //tableau de config  recons
    modif_tab_val("DELTA_NMAX",editDnMax->GetValue().ToStdString(),tab_val_recon);
    modif_tab_val("DELTA_NMIN",editDnMin->GetValue().ToStdString(),tab_val_recon);
    modif_tab_val("KAPPA_MAX",editKappaMax->GetValue().ToStdString(),tab_val_recon);
    modif_tab_val("KAPPA_MIN",editKappaMin->GetValue().ToStdString(),tab_val_recon);
    modif_tab_val("NB_ITER_GPS",editIteration->GetValue().ToStdString(),tab_val_recon);

    modif_tab_val("CIRCLE_CX",editCX->GetValue().ToStdString(),tab_val_manip);
    modif_tab_val("CIRCLE_CY",editCY->GetValue().ToStdString(),tab_val_manip);


    sav_val(chemin_recon,tab_val_recon);
    sav_val(chemin_config_manip,tab_val_manip);//sav PC acquisition
    sav_val(chemin_config_GUI,tab_val_gui_conf);
}



float Gui_TomoFrame::extract_val(std::string token,  std::string chemin_fic)
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
    fichier.close();
    return valeur;
}
//superflu car déjà géré sav_all ?
void Gui_TomoFrame::OnComboScan(wxCommandEvent &even)
{
    modif_tab_val("SCAN_PATTERN", ComboBox_Scan->GetString(ComboBox_Scan->GetSelection()).ToStdString(),tab_val_manip);
    cout<<"Schéma: "<<ComboBox_Scan->GetString(ComboBox_Scan->GetSelection())<<endl;
    if(ComboBox_Scan->GetString(ComboBox_Scan->GetSelection())==wxT("ANNULAR"))
    {

        cout<<"schema ANNULAR!"<<endl;
        editnbCirclesAnnular->Enable(true);
    }
    else
    {
        editnbCirclesAnnular->Enable(false);
    }
//wxArrayString  liste_scan(4);
//liste_scan=ComboBox_Scan->GetStrings();
//cout<<"liste_scan[2]"<<liste_scan[2]<<endl;
//cout<<"Indice  : "<<ComboBox_Scan->GetCurrentSelection ()<<endl;
//cout<<"GetStringSelection() : "<<ComboBox_Scan->GetStringSelection()<<endl;
//cout<<"GetValue() : "<<ComboBox_Scan->GetValue()<<endl;
}
