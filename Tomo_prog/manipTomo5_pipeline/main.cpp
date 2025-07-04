// *****************************************************************************
//
//      Copyright (c) 2013, Pleora Technologies Inc., All rights reserved.
//
// *****************************************************************************

//
// Shows how to use a PvStream object to acquire images from a GigE Vision or
// USB3 Vision device.
//
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>

using namespace cv;
using namespace std;
#include <string>
#include <sstream>
//pleora lib
#include <PvSampleUtils.h>
#include <PvDevice.h>
#include <PvDeviceGEV.h>
#include <PvDeviceU3V.h>
#include <PvStream.h>
#include <PvStreamGEV.h>
#include <PvStreamU3V.h>
#include <PvPipeline.h>
#include <PvBuffer.h>

#include "Ljack.h"
#include "Ljack_DAC.h"
#include "macros.h"
#include "vectra.h"
#include "string.h"
#include "vChronos.h"
#include <boost/thread.hpp>
#include <boost/date_time.hpp>
#include "vPlotterT.h"

//#include "vImg.hpp"

#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <termios.h>
#include <unistd.h>
#include <fcntl.h>

#include "fonctions.h"
#include "scan_functions.h"



#include <PvDeviceFinderWnd.h>


#include <list>
typedef std::list<PvBuffer *> BufferList;
//float tiptilt_factor = 0.25;
//float tiptilt_factor_x;
//float tiptilt_factor_y;
//float2D VfOffset;
float tiptilt_factor_x = 0;//init rayon X en volt
float tiptilt_factor_y = 0;//init rayon  Y en volt
float2D VfOffset= {-0, -0};
unsigned int b_AmpliRef=0;
float NAcondLim=1;//coefficient de limitation de 'louverture numerique de balyage 0<valeur<=1
float flower_x(float t);
float flower_y(float t);
size_t MAX_IMAGES = 0;

PV_INIT_SIGNAL_HANDLER();

#define BUFFER_COUNT ( 2 )

///
/// Function Prototypes
///
#ifdef PV_GUI_NOT_AVAILABLE
const PvDeviceInfo *SelectDevice( PvSystem *aPvSystem );
#else
const PvDeviceInfo *SelectDevice( PvDeviceFinderWnd *aDeviceFinderWnd );
#endif // PV_GUI_NOT_AVAILABLE
PvDevice *ConnectToDevice( const PvDeviceInfo *aDeviceInfo );
PvStream *OpenStream( const PvDeviceInfo *aDeviceInfo );
void ConfigureStream( PvDevice *aDevice, PvStream *aStream );
PvPipeline* CreatePipeline( PvDevice *aDevice, PvStream *aStream );
void AcquireImages( PvDevice *aDevice, PvStream *aStream, PvPipeline *aPipeline, Ljack_DAC *DAC_tiptilt,string chemin_result, string chemin_acquis, int dimROI, vector<float2D> &Vout_table );
void AcquireIntensiteRef( PvDevice *aDevice, PvStream *aStream, PvPipeline *aPipeline, Ljack_DAC *DAC_tiptilt, string chemin_result, string chemin_acquis, int dimROI);
string IntToString ( int number );


// Main function
//
int main(int argc, char *argv[])
{

//    usage(argc, argv);
    /*  while (argc > 0) {
                  if (!strcmp(argv[0], "-ni") && (argc > 1)) {
                          MAX_IMAGES = atoi(argv[1]);
                          argc--;
                          argv++;
                  }
                  if (!strcmp(argv[0], "-vfleur") && (argc > 2)) {
                          tiptilt_factor_x = atof(argv[1])/10;
                          tiptilt_factor_y = atof(argv[2])/10;
                          argc-=2;
                          argv+=2;
                  }
                  if (!strcmp(argv[0], "-voffset") && (argc > 2)) {
                          VfOffset.x = atof(argv[1]);
                          VfOffset.y = atof(argv[2]);

                          argc-=2;
                          argv+=2;
                  }

                  argc--;
                  argv++;
      }*/
    cout<<"########################## ACQUISITION ##################################"<<endl;
    ///Initialiser les deux tableaux stockant les valeurs de recon.txt et config_manip.txt en mémoire
    string chemin_config_GUI,chemin_recon, chemin_config_manip,chemin_acquis,chemin_result,repertoire_config;
    string home=getenv("HOME");
    string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
    chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
    cout<<"lecture des chemins dans le fichier "<<chemin_config_GUI<<endl;

    // repertoire_config=extract_string("CHEMIN_CONFIG_PC_ACQUIS",home+fin_chemin_gui_tomo);///chemin de config pour la 1ere acquisition, qui va servir à copier recon.txt et
    repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);///chemin de config pour la 1ere acquisition, qui va servir à copier recon.txt et
    ///config_manip.txt dans le repertoire d'acquisition.
    chemin_acquis=extract_string("CHEMIN_ACQUIS",home+fin_chemin_gui_tomo);
    chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);
    chemin_recon=repertoire_config+"/recon.txt";

    chemin_config_manip=repertoire_config+"/config_manip.txt";
    cout<<"fichier config : "<<chemin_config_manip<<endl;

    int dimROI=extract_val("DIM_ROI", chemin_config_manip);
    int NXMAX=extract_val("NXMAX", chemin_config_manip);

    string repAcquis=chemin_acquis;
    Var2D dimHolo= {2*NXMAX,2*NXMAX};
    efface_acquis(repAcquis,chemin_config_manip,chemin_recon);
///!\ attention les fonction labjack de john n'utilisent pas des tensions mais des facteurs "d'attenuation de la tension max". Tout est divisé par 10.
    MAX_IMAGES=extract_val("NB_HOLO",chemin_config_manip);
    NAcondLim=extract_val("NA_COND_LIM",chemin_config_manip);
    tiptilt_factor_x=(abs(extract_val("VXMIN", chemin_config_manip))+abs(extract_val("VXMAX", chemin_config_manip)))*NAcondLim/20;///amplitude de tension/10=(|vmax|+|vmin|)/10
    tiptilt_factor_y=(abs(extract_val("VYMIN", chemin_config_manip))+abs(extract_val("VYMAX", chemin_config_manip)))*NAcondLim/20;

    VfOffset.x=(extract_val("VXMIN", chemin_config_manip)*NAcondLim+tiptilt_factor_x*10);///offset en volt=amplitude tension+Vmin
    VfOffset.y=(extract_val("VYMIN", chemin_config_manip)*NAcondLim+tiptilt_factor_y*10)*NAcondLim;
    b_AmpliRef=extract_val("AMPLI_REF", chemin_recon);

    cout<<"  ####################################"<<endl;
    cout<<"  #  tension MAX x="<<tiptilt_factor_x*10<<"                   #"<<endl;
    cout<<"  #  tension MAX y="<<tiptilt_factor_y*10<<"                     #"<<endl;
    cout<<"  #  VfOffset.x="<<VfOffset.x<<" | VfOffset.y="<<VfOffset.y<<"  #"<<endl;
    cout<<"  #  Nombre d'hologrammes="<<MAX_IMAGES<<"        #"<<endl;
    cout<<"  ####################################" <<endl;



    ///Copie des fichiers de config (depuis repertoire PC acquis : projet_tomo) dans le repertoire d'acquisition


    /*  std::ifstream srce_recon(chemin_recon, std::ios::binary ) ;
      string chemin_recon_dest=chemin_acquis+"/recon.txt";
      cout<<"srce="<<chemin_recon<<endl;
      std::ofstream dest_recon(chemin_recon_dest, std::ios::binary ) ;
      dest_recon << srce_recon.rdbuf();
      dest_recon.close();

      std::ifstream srce_config(chemin_config_manip, std::ios::binary ) ;

      string chemin_config_dest=chemin_acquis+"/config_manip.txt";
      std::ofstream dest_config(chemin_config_dest, std::ios::binary ) ;
      dest_config << srce_config.rdbuf();
      ///-----------------------------------------------------------------

      dest_config.close();*/

    ///initialiser la caméra (détection par PvDeviceFinderWnd, sélection d'un périphérique)
    const PvDeviceInfo *lDeviceInfo = NULL;
    PvDevice *lDevice = NULL;
    PvStream *lStream = NULL;

    PV_SAMPLE_INIT();

    PvDeviceFinderWnd *lDeviceFinderWnd = new PvDeviceFinderWnd();
    lDeviceInfo = SelectDevice( lDeviceFinderWnd );

///----déclaration LabJack
    Ljack LU3HV;

    //Ljack_DAC ljDAC_phase(LU3HV);
   // ASSERT(ljDAC_phase.connect(FIO7_6));

    Ljack_DAC ljDAC_flower(LU3HV);
    ASSERT(ljDAC_flower.connect(FIO5_4));

    ljDAC_flower.set_A_output(VfOffset.x);
    ljDAC_flower.set_B_output(VfOffset.y);

    //generate table containing all the voltages
    string scan_pattern_str=extract_string("SCAN_PATTERN",chemin_config_manip);
    int NbCirclesAnnular=extract_val("NB_CIRCLES_ANNULAR", chemin_config_manip);
    cout<<"Scan pattern = "<<scan_pattern_str<<endl;
    int const rho=10;//Vmax admissible par labjack (la tension max du pattern est géré par le facteur "tiptilt_factor").
    int const nbHolo=MAX_IMAGES;
    vector<float2D> Vout_table(nbHolo);
    bool scan_pattern_exist=0;
    if(scan_pattern_str=="ROSACE") {scan_fleur(rho, nbHolo, Vout_table); scan_pattern_exist=1; }
    if(scan_pattern_str=="FERMAT") {scan_fermat(rho, nbHolo, Vout_table,dimHolo);  scan_pattern_exist=1;}
    if(scan_pattern_str=="ARCHIMEDE"){scan_spiralOS(rho, nbHolo, Vout_table,4); scan_pattern_exist=1;};
    if(scan_pattern_str=="ANNULAR") {scan_annular(rho, NbCirclesAnnular, nbHolo, Vout_table); scan_pattern_exist=1;}
    if(scan_pattern_str=="UNIFORM3D"){scan_uniform3D(rho, nbHolo, Vout_table); scan_pattern_exist=1;};
    if(scan_pattern_str=="RANDOM_POLAR"){ scan_random_polar3D(rho, nbHolo, Vout_table,dimHolo);scan_pattern_exist=1;}
    if(scan_pattern_str=="STAR"){scan_etoile(rho, nbHolo,  3, Vout_table);scan_pattern_exist=1;}
    //-------------------------
    //unidentified scan pattern : warn user
    cout<<"scan_pattern_exist="<<scan_pattern_exist<<endl;
    if(scan_pattern_exist==0)
        cout<<"Cannot handle this scan pattern"<<endl;
///-----
    cout << "PvStreamSample:" << endl << endl;

    if ( NULL != lDeviceInfo )
        {
            lDevice = ConnectToDevice( lDeviceInfo );
            if ( lDevice != NULL )
                {
                    ///-------------------------------------------------------changer la largeur de la ROI---------------------------------------

                    PvGenParameterArray* lParameters = lDevice->GetParameters();///stocker tous les paramètres caméra
                    // Get width parameter - mandatory GigE Vision parameter, it should be there.
                    PvGenParameter *lParameter = lParameters->Get( "Width" );///stocker un seul paramètre caméra
                    PvGenInteger *lWidthParameter = dynamic_cast<PvGenInteger *>( lParameter );
                    //lStream = OpenStream( lDeviceInfo );
                    if ( lWidthParameter == NULL )
                        {
                            cout << "Unable to get the width parameter." << endl;
                        }
                    int64_t lWidth=dimROI;
                    lWidthParameter->SetValue( lWidth );

                    ///changer la hauteur de la ROI---------
                    lParameter = lParameters->Get( "Height" );
                    PvGenInteger *lHeightParameter = dynamic_cast<PvGenInteger *>( lParameter );
                    if ( lHeightParameter == NULL )
                        {
                            cout << "Unable to get the width parameter." << endl;
                        }
                    int64_t lHeight=dimROI;
                    lHeightParameter->SetValue( lHeight );

                    ///changer OffsetX de la ROI---------
                    lParameter = lParameters->Get( "OffsetX" );
                    PvGenInteger *lOffsetXParameter = dynamic_cast<PvGenInteger *>( lParameter );

                    if ( lOffsetXParameter == NULL )
                        {
                            cout << "Unable to get the OffsetX parameter." << endl;
                        }
                    int64_t lOffsetX=round(dimROI/2);///OFFSET de l'image=512
                    lOffsetXParameter->SetValue( lOffsetX );

                    ///changer OffsetY de la ROI---------
                    lParameter = lParameters->Get( "OffsetY" );
                    PvGenInteger *lOffsetYParameter = dynamic_cast<PvGenInteger *>( lParameter );
                    if ( lOffsetYParameter == NULL )
                        {
                            cout << "Unable to get the OffsetY parameter." << endl;
                        }
                    int64_t lOffsetY=round(dimROI/2);///OFFSET de l'image=512
                    lOffsetYParameter->SetValue( lOffsetY );

                    ///--------------------------------------------------------fin changement ROI------------------------------------------------

                    lStream = OpenStream( lDeviceInfo );
                    if ( lStream != NULL )
                        {
                            PvPipeline *lPipeline = NULL;

                            ConfigureStream( lDevice, lStream );
                            lPipeline = CreatePipeline( lDevice, lStream );
                            lPipeline->SetBufferCount(32);
                            if( lPipeline )
                                {
                                    AcquireImages(lDevice, lStream, lPipeline, &ljDAC_flower, chemin_result,chemin_acquis, dimROI,Vout_table);

                                    if(b_AmpliRef==1)
                                        AcquireIntensiteRef( lDevice, lStream, lPipeline, &ljDAC_flower, chemin_result,chemin_acquis, dimROI );
                                    delete lPipeline;
                                }
                            // Close the stream
                            cout << "Closing stream" << endl;
                            lStream->Close();
                            PvStream::Free( lStream );
                        }

                    ///-----Rétablir la valeur Max Largeur et hauteur de la caméra
                    int64_t lMaxWidth = 0;
                    if ( !(lWidthParameter->GetMax( lMaxWidth ).IsOK()) )//met à jour et controle si tout s'est bien passé.
                        {
                            cout << "Error retrieving width max from device" << endl;
                        }
                    lWidthParameter->SetValue( lMaxWidth );
                    int64_t lMaxHeight = 0;
                    if ( !(lHeightParameter->GetMax( lMaxHeight ).IsOK()) )//met à jour et controle si tout s'est bien passé.
                        {
                            cout << "Error retrieving width max from device" << endl;
                        }
                    lHeightParameter->SetValue( lMaxHeight );

                    // Disconnect the device
                    cout << "Disconnecting device" << endl;
                    // cout << "<press a key to exit>" << endl;
                    // PvWaitForKeyPress();
                    lDevice->Disconnect();
                    PvDevice::Free( lDevice );
                }
        }

    // cout << "<press a key to exit>" << endl;
    // PvWaitForKeyPress();

    delete lDeviceFinderWnd;
    lDeviceFinderWnd = NULL;

    PV_SAMPLE_TERMINATE();

    ljDAC_flower.set_A_output(0.0);
    ljDAC_flower.set_B_output(0.0);
  ofstream file;
  file.open ("codebind.txt");
  file << "Please writr this text to a file.\n this text is written using C++\n";
  file.close();
    return 0;
}

///Fonctions

const PvDeviceInfo *SelectDevice( PvDeviceFinderWnd *aDeviceFinderWnd )
{
    const PvDeviceInfo *lDeviceInfo = NULL;
    PvResult lResult;

    if (NULL != aDeviceFinderWnd)
        {
            // Display list of GigE Vision and USB3 Vision devices
            lResult = aDeviceFinderWnd->ShowModal();
            if ( !lResult.IsOK() )
                {
                    // User hit cancel
                    cout << "No device selected." << endl;
                    return NULL;
                }

            // Get the selected device information.
            lDeviceInfo = aDeviceFinderWnd->GetSelected();
        }
    return lDeviceInfo;
}

PvDevice *ConnectToDevice(const PvDeviceInfo *aDeviceInfo )
{
    PvDevice *lDevice;
    PvResult lResult;

    // Connect to the GigE Vision or USB3 Vision device
    cout << "Connecting to " << aDeviceInfo->GetDisplayID().GetAscii() << "." << endl;
    lDevice = PvDevice::CreateAndConnect( aDeviceInfo, &lResult );
    if ( lDevice == NULL )
        {
            cout << "Unable to connect to " << aDeviceInfo->GetDisplayID().GetAscii() << "." << endl;
        }

    return lDevice;
}

PvStream *OpenStream( const PvDeviceInfo *aDeviceInfo )
{
    PvStream *lStream;
    PvResult lResult;

    // Open stream to the GigE Vision or USB3 Vision device
    cout << "Opening stream to device." << endl;
    lStream = PvStream::CreateAndOpen( aDeviceInfo->GetConnectionID(), &lResult );
    if ( lStream == NULL )
        {
            cout << "Unable to stream from " << aDeviceInfo->GetDisplayID().GetAscii() << "." << endl;
        }

    return lStream;
}

void ConfigureStream( PvDevice *aDevice, PvStream *aStream )
{
    // If this is a GigE Vision device, configure GigE Vision specific streaming parameters
    PvDeviceGEV* lDeviceGEV = dynamic_cast<PvDeviceGEV *>( aDevice );
    if ( lDeviceGEV != NULL )
        {
            PvStreamGEV *lStreamGEV = static_cast<PvStreamGEV *>( aStream );

            // Negotiate packet size
            lDeviceGEV->NegotiatePacketSize();

            // Configure device streaming destination
            lDeviceGEV->SetStreamDestination( lStreamGEV->GetLocalIPAddress(), lStreamGEV->GetLocalPort() );
        }
}

PvPipeline *CreatePipeline( PvDevice *aDevice, PvStream *aStream )
{
    // Create the PvPipeline object
    PvPipeline* lPipeline = new PvPipeline( aStream );

    if ( lPipeline != NULL )
        {
            // Reading payload size from device
            uint32_t lSize = aDevice->GetPayloadSize();

            // Set the Buffer count and the Buffer size
            lPipeline->SetBufferCount( BUFFER_COUNT );
            lPipeline->SetBufferSize( lSize );
        }

    return lPipeline;
}

void AcquireImages( PvDevice *aDevice, PvStream *aStream, PvPipeline *aPipeline, Ljack_DAC *DAC_tiptilt, string chemin_result, string chemin_acquis, int dimROI, vector<float2D> &Vout_table)
{
    Mat image;
    Mat ImFleur(1000,1000, CV_16U, Scalar(0));
    float theta=0.393;//init angle courbes parametré, rayon max en volt
    int nbHolo=MAX_IMAGES;
    // Get device parameters need to control streaming
    PvGenParameterArray *lDeviceParams = aDevice->GetParameters();


    // Map the GenICam AcquisitionStart and AcquisitionStop commands
    PvGenCommand *lStart = dynamic_cast<PvGenCommand *>( lDeviceParams->Get( "AcquisitionStart" ) );
    PvGenCommand *lStop = dynamic_cast<PvGenCommand *>( lDeviceParams->Get( "AcquisitionStop" ) );

    // Note: the pipeline must be initialized before we start acquisition
    cout << "Starting pipeline" << endl;
    aPipeline->Start();

    // Get stream parameters
    PvGenParameterArray *lStreamParams = aStream->GetParameters();

    // Map a few GenICam stream stats counters
    PvGenFloat *lFrameRate = dynamic_cast<PvGenFloat *>( lStreamParams->Get( "AcquisitionRate" ) );
    PvGenFloat *lBandwidth = dynamic_cast<PvGenFloat *>( lStreamParams->Get( "Bandwidth" ) );

    ///changer la ROI systématiquement à 1024-----------------------------------

    // Get width parameter - mandatory GigE Vision parameter, it should be there.
    PvGenParameter *lParameter = lDeviceParams->Get( "Width" );
    PvGenInteger *lWidthParameter = dynamic_cast<PvGenInteger *>( lParameter );
//        lStream = OpenStream( lDeviceInfo );
    if ( lWidthParameter == NULL )
        {
            cout << "Unable to get the width parameter." << endl;
        }

    int64_t lWidth=dimROI;///1024
    lWidthParameter->SetValue( lWidth );

    //changer la hauteur de la ROI---------

    // Get ROI width parameter -
    lParameter =lDeviceParams->Get( "Height" );
    PvGenInteger *lHeightParameter = dynamic_cast<PvGenInteger *>( lParameter );

    if ( lHeightParameter == NULL )
        {
            cout << "Unable to get the width parameter." << endl;
        }
    int64_t lHeight=dimROI;///1024;
    lHeightParameter->SetValue( lHeight );

    ///-------fin changement ROI------------------------------------------


    // Enable streaming and send the AcquisitionStart command
    cout << "Enabling streaming and sending AcquisitionStart command." << endl;
    aDevice->StreamEnable();
    lStart->Execute();

    char lDoodle[] = "|\\-|-/";
    int lDoodleIndex = 0;
    double lFrameRateVal = 0.0;
    double lBandwidthVal = 0.0;
    int cpt=0;

    /// Acquire images until the user instructs us to stop.

    cout << endl << "##-- Acquisition Hologrammes --##" << endl;
    ////////////////////////////////////////////////////////////////////////////////////////
    /// BOUCLE ACQUISITION

    // synchro instrumentale:
    size_t cpt_img=0;

    boost::posix_time::milliseconds delay_jack( 0 ); //1.2 + 1.2 +1 );
    boost::posix_time::milliseconds delay_tilt( 0);//attention bug : boost 1.67 n'accepte pas de flottant !
    boost::posix_time::milliseconds delay_cam( 0 );

    vChronos vTime("prise d'images--");
    vTime.clear();

    vTime.start();    //while ( !PvKbHit() )//tant que touche pas pressée
    //FlowerPlotter -> compute_next();
    // float Vfleur.x = FlowerPlotter -> get_x() * tiptilt_factor;

    // float Vfleur.y = FlowerPlotter -> get_y() * tiptilt_factor;

    float2D Vfleur= {0*tiptilt_factor_x,0*tiptilt_factor_y};

//    float2D Vfleur={6,6};//init pour intensite ref
    //boost::this_thread::sleep( delay_jack );
    // cout<<"Valeur init Vx, Vy="<<Vfleur.x<<","<<Vfleur.y<<endl;
    DAC_tiptilt -> set_A_output(Vfleur.x*tiptilt_factor_x+VfOffset.x); // + offset_vx);
    DAC_tiptilt -> set_B_output(Vfleur.y*tiptilt_factor_y+VfOffset.y); // + offset_vy);

    //scan_fleur(rho,MAX_IMAGES,Vout_table);
    //while ( (! _kbhit())  && (cpt_img < MAX_IMAGES) )
    while ( (cpt_img < MAX_IMAGES) )
        {

            PvBuffer *lBuffer = NULL;
            PvResult lOperationResult;


            // Retrieve next buffer
            PvResult lResult = aPipeline->RetrieveNextBuffer( &lBuffer, 1000, &lOperationResult );
            //boost::this_thread::sleep( delay_cam );


            if ( lResult.IsOK() )///Buffer acquis?
                {
                    if ( lOperationResult.IsOK() ) ///opération finie sans erreur? (time out..)
                        {
                            PvPayloadType lType;
                            // We now have a valid buffer. This is where you would typically process the buffer.
                            lFrameRate->GetValue( lFrameRateVal );

                            lBandwidth->GetValue( lBandwidthVal );

                            // If the buffer contains an image, display width and height.
                            uint32_t lWidth = 0, lHeight = 0;
                            lType = lBuffer->GetPayloadType();

                            cout << fixed << setprecision( 1 );
                            cout << lDoodle[ lDoodleIndex ];
                            // cout << " BlockID: " << uppercase << hex << setfill( '0' ) << setw( 16 ) << lBuffer->GetBlockID();

                            if ( lType == PvPayloadTypeImage )///le buffer contient bien une image? Tout est ok, on traite le flux
                                {
                                    PvImage *lImage = lBuffer->GetImage();// Get image specific buffer interface.

                                    unsigned char *img = lImage->GetDataPointer();// Get image data pointer so we can pass it to CV::MAT container

                                    // Read width, height.
                                    lWidth = lImage->GetWidth();
                                    lHeight = lImage->GetHeight();
                                    cout << "  W: " << dec << lWidth << " H: " << lHeight;

                                    // Copy/convert Pleora Vision image pointer to cv::Mat container
                                    cv::Mat lframe(lHeight,lWidth,CV_8UC1,img, cv::Mat::AUTO_STEP);
                                    // namedWindow( "Display window", WINDOW_AUTOSIZE );// Create a window for display.
                                    // imshow("Display window",lframe);
                                    // cv::waitKey(1);
                                    imwrite(chemin_acquis+format("/i%03d.pgm", cpt_img), lframe);

                                    //Vfleur=maj_fleur(Vfleur, rho, nbHolo, &theta);///param Vfleur inutile, à supprimer ?
                                    Vfleur.x=Vout_table[cpt_img].x;
                                    Vfleur.y=Vout_table[cpt_img].y;
                                   // cout<<"vout.x"<<Vout_table[cpt_img].x<<endl;
                                   // cout<<"vout.y"<<Vout_table[cpt_img].y<<endl;
                                  //  cout<<"vfleur.x"<<Vfleur.x<<endl;
                                  //  cout<<"vfleur.y"<<Vfleur.y<<endl;
                                   // cout<<endl;
                                    DAC_tiptilt -> set_A_output(Vfleur.x*tiptilt_factor_x+VfOffset.x);
                                    DAC_tiptilt -> set_B_output(Vfleur.y*tiptilt_factor_y+VfOffset.y);

                                   //  cout<<Vfleur.x*tiptilt_factor_x+VfOffset.x<<endl;
                                   //  cout<<Vfleur.y*tiptilt_factor_y+VfOffset.y<<endl;
                                    cpt_img++; //incrémenter le numero de l'image
                                    //   ImFleur.at<unsigned short int>(100*(Vfleur.x*tiptilt_factor_x+VfOffset.x)+325,100*(Vfleur.y*tiptilt_factor_y+VfOffset.y)+650)=cpt_img;
                                    //  boost::this_thread::sleep( delay_tilt );

                                }
                            else
                                {
                                    cout << " (buffer does not contain image)";
                                }
                            cout << "  " << lFrameRateVal << " FPS  " << ( lBandwidthVal / 1000000.0 ) << " Mb/s   \r";
                        }
                    else
                        {
                            // Non OK operational result
                            cout << lDoodle[ lDoodleIndex ] << " " << lOperationResult.GetCodeString().GetAscii() << "\r";
                        }

                    // Release the buffer back to the pipeline
                    aPipeline->ReleaseBuffer( lBuffer );
                }
            else
                {
                    // Retrieve buffer failure
                    cout << lDoodle[ lDoodleIndex ] << " " << lResult.GetCodeString().GetAscii() << "\r";
                }

            ++lDoodleIndex %= 6;
        }



    // imwrite(chemin_result+"/balayage_fleur.png", ImFleur);
    //  PvGetChar(); // Flush key buffer for next stop.
    //  cout << endl << endl;


    // Tell the device to stop sending images.
    cout << "Sending AcquisitionStop command to the device" << endl;
    lStop->Execute();

    // Disable streaming on the device
    cout << "Disable streaming on the controller." << endl;
    aDevice->StreamDisable();


///-----Rétablir la valeur Max Largeur et hauteur de la caméra
    int64_t lMaxWidth = 0;
    if ( !(lWidthParameter->GetMax( lMaxWidth ).IsOK()) )//met à jour et controle si tout s'est bien passé.
        {
            cout << "Error retrieving width max from device" << endl;
        }
    lWidthParameter->SetValue( lMaxWidth );
    int64_t lMaxHeight = 0;
    if ( !(lHeightParameter->GetMax( lMaxHeight ).IsOK()) )//met à jour et controle si tout s'est bien passé.
        {
            cout << "Error retrieving width max from device" << endl;
        }
    lHeightParameter->SetValue( lMaxHeight );
    //-----Fin Rétablir la valeur Max Largeur et hauteur de la caméra



    // Stop the pipeline
    cout << "Stop pipeline" << endl;
    aPipeline->Stop();
}

std::string IntToString ( int number )
{
    std::ostringstream oss;

    // Works just like cout
    oss<< number;

    // Return the underlying string
    return oss.str();
}

// *****************************************************************************

float rayon_balayage = 10.00f;

float
flower_x(float t)
{
    float c4t10 = rayon_balayage * cos(4 * t);
    return c4t10 * cos(t);
}


float
flower_y(float t)
{
    float c4t10 = rayon_balayage * cos(4 * t);
    return c4t10 * sin(t);
}

void AcquireIntensiteRef( PvDevice *aDevice, PvStream *aStream, PvPipeline *aPipeline, Ljack_DAC *DAC_tiptilt, string chemin_result, string chemin_acquis, int dimROI)
{
    DAC_tiptilt -> set_A_output(0); // + offset_vx);
    DAC_tiptilt -> set_B_output(8); // + offset_vy);
    Mat image;
    // Get device parameters need to control streaming
    PvGenParameterArray *lDeviceParams = aDevice->GetParameters();
    // Map the GenICam AcquisitionStart and AcquisitionStop commands
    PvGenCommand *lStart = dynamic_cast<PvGenCommand *>( lDeviceParams->Get( "AcquisitionStart" ) );
    PvGenCommand *lStop = dynamic_cast<PvGenCommand *>( lDeviceParams->Get( "AcquisitionStop" ) );

    // Note: the pipeline must be initialized before we start acquisition
    cout << "Starting pipeline" << endl;
    aPipeline->Start();

    // Get stream parameters
    PvGenParameterArray *lStreamParams = aStream->GetParameters();

    // Map a few GenICam stream stats counters
    PvGenFloat *lFrameRate = dynamic_cast<PvGenFloat *>( lStreamParams->Get( "AcquisitionRate" ) );
    PvGenFloat *lBandwidth = dynamic_cast<PvGenFloat *>( lStreamParams->Get( "Bandwidth" ) );

    ///changer la ROI systématiquement à 1024-----------------------------------

    // Get width parameter - mandatory GigE Vision parameter, it should be there.
    PvGenParameter *lParameter = lDeviceParams->Get( "Width" );
    PvGenInteger *lWidthParameter = dynamic_cast<PvGenInteger *>( lParameter );
    //lStream = OpenStream( lDeviceInfo );
    if ( lWidthParameter == NULL )
        {
            cout << "Unable to get the width parameter." << endl;
        }

    int64_t lWidth=dimROI;///1024
    lWidthParameter->SetValue( lWidth );

    //changer la hauteur de la ROI---------

    // Get ROI width parameter -
    lParameter =lDeviceParams->Get( "Height" );
    PvGenInteger *lHeightParameter = dynamic_cast<PvGenInteger *>( lParameter );

    if ( lHeightParameter == NULL )
        {
            cout << "Unable to get the width parameter." << endl;
        }
    int64_t lHeight=dimROI;///1024;
    lHeightParameter->SetValue( lHeight );

    ///-------fin changement ROI------------------------------------------

    ///changer OffsetX de la ROI---------
    lParameter = lDeviceParams->Get( "OffsetX" );
    PvGenInteger *lOffsetXParameter = dynamic_cast<PvGenInteger *>( lParameter );

    if ( lOffsetXParameter == NULL )
        {
            cout << "Unable to get the OffsetX parameter." << endl;
        }
    int64_t lOffsetX=round(dimROI/2);///OFFSET de l'image=512
    lOffsetXParameter->SetValue( lOffsetX );

    ///changer OffsetY de la ROI---------
    lParameter = lDeviceParams->Get( "OffsetY" );
    PvGenInteger *lOffsetYParameter = dynamic_cast<PvGenInteger *>( lParameter );
    if ( lOffsetYParameter == NULL )
        {
            cout << "Unable to get the OffsetY parameter." << endl;
        }
    int64_t lOffsetY=round(dimROI/2);///OFFSET de l'image=512
    lOffsetYParameter->SetValue( lOffsetY );
///----------------------------------fin changement offset
    // Enable streaming and send the AcquisitionStart command
    aDevice->StreamEnable();
    lStart->Execute();

    char lDoodle[] = "|\\-|-/";
    int lDoodleIndex = 0;
    double lFrameRateVal = 0.0;
    double lBandwidthVal = 0.0;
    int cpt=0;

    /// Acquire images until the user instructs us to stop.
    cout << endl << "##-- Acquisition Référence --##" << endl;

    ////////////////////////////////////////////////////////////////////////////////////////
    /// BOUCLE ACQUISITION

    // synchro instrumentale:
    size_t cpt_img=0;

    boost::posix_time::milliseconds delay_jack( 1 ); //1.2 + 1.2 +1 );
    boost::posix_time::milliseconds delay_tilt( 1);
    boost::posix_time::milliseconds delay_cam( 100 );


    PvBuffer *lBuffer = NULL;
    PvResult lOperationResult;

    // Retrieve next buffer
    PvResult lResult = aPipeline->RetrieveNextBuffer( &lBuffer, 1000, &lOperationResult );
    boost::this_thread::sleep( delay_cam );


    if ( lResult.IsOK() )///Buffer acquis?
        {
            if ( lOperationResult.IsOK() ) ///opération finie sans erreur? (time out..)
                {
                    PvPayloadType lType;
                    lFrameRate->GetValue( lFrameRateVal );
                    lBandwidth->GetValue( lBandwidthVal );

                    // If the buffer contains an image, display width and height.
                    uint32_t lWidth = 0, lHeight = 0;
                    lType = lBuffer->GetPayloadType();

                    cout << fixed << setprecision( 1 );
                    cout << lDoodle[ lDoodleIndex ];
                    // cout << " BlockID: " << uppercase << hex << setfill( '0' ) << setw( 16 ) << lBuffer->GetBlockID();

                    if ( lType == PvPayloadTypeImage )///le buffer contient bien une image? Tout est ok, on traite le flux
                        {
                            DAC_tiptilt -> set_A_output(9.5);
                            DAC_tiptilt -> set_B_output(9.5);

                            // Get image specific buffer interface.
                            PvImage *lImage = lBuffer->GetImage();

                            // Get image data pointer so we can pass it to CV::MAT container
                            unsigned char *img = lImage->GetDataPointer();

                            // Read width, height.
                            lWidth = lImage->GetWidth();
                            lHeight = lImage->GetHeight();
                            cout << "  W: " << dec << lWidth << " H: " << lHeight;

                            // Copy/convert Pleora Vision image pointer to cv::Mat container
                            cv::Mat lframe(lHeight,lWidth,CV_8UC1,img, cv::Mat::AUTO_STEP);

                            imwrite(chemin_acquis+"/Intensite_ref.pgm", lframe);
                            cout<<endl;
                            cout<<"Acquisition intensite reference dans "<<chemin_acquis<<endl;
                        }
                    else
                        {
                            cout << " (buffer does not contain image)";
                        }
                    cout << "  " << lFrameRateVal << " FPS  " << ( lBandwidthVal / 1000000.0 ) << " Mb/s   \r";
                }
            else
                {
                    // Non OK operational result
                    cout << lDoodle[ lDoodleIndex ] << " " << lOperationResult.GetCodeString().GetAscii() << "\r";
                }
            // Release the buffer back to the pipeline
            aPipeline->ReleaseBuffer( lBuffer );
        }
    else
        {
            // Retrieve buffer failure
            cout << lDoodle[ lDoodleIndex ] << " " << lResult.GetCodeString().GetAscii() << "\r";
        }

    ++lDoodleIndex %= 6;

    // imwrite(chemin_result+"/balayage_fleur.png", ImFleur);
    //  PvGetChar(); // Flush key buffer for next stop.
    //  cout << endl << endl;

    // Tell the device to stop sending images.
    cout << "Sending AcquisitionStop command to the device" << endl;
    lStop->Execute();

    // Disable streaming on the device
    cout << "Disable streaming on the controller." << endl;
    aDevice->StreamDisable();


///-----Rétablir la valeur Max Largeur et hauteur de la caméra
    int64_t lMaxWidth = 0;
    if ( !(lWidthParameter->GetMax( lMaxWidth ).IsOK()) )//met à jour et controle si tout s'est bien passé.
        {
            cout << "Error retrieving width max from device" << endl;
        }
    lWidthParameter->SetValue( lMaxWidth );
    int64_t lMaxHeight = 0;
    if ( !(lHeightParameter->GetMax( lMaxHeight ).IsOK()) )//met à jour et controle si tout s'est bien passé.
        {
            cout << "Error retrieving width max from device" << endl;
        }
    lHeightParameter->SetValue( lMaxHeight );
    //-----Fin Rétablir la valeur Max Largeur et hauteur de la caméra



    // Stop the pipeline
    cout << "Stop pipeline" << endl;
    aPipeline->Stop();
}

