// *****************************************************************************
//
//      Copyright (c) 2013, Pleora Technologies Inc., All rights reserved.
//
// *****************************************************************************

//
// Shows how to use a PvStream object to acquire images from a GigE Vision or
// USB3 Vision device.
//

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <iostream>
using namespace cv;
using namespace std;
#include <string>
#include <sstream>
///Genicam/pleora
#include <PvSampleUtils.h>
#include <PvDevice.h>
#include <PvDeviceGEV.h>
#include <PvDeviceU3V.h>
#include <PvStream.h>
#include <PvStreamGEV.h>
#include <PvStreamU3V.h>
#include <PvPipeline.h>
#include <PvBuffer.h>
///labjack
//#include "Ljack.h"
//#include "Ljack_DAC.h"
///-----------------
//#include "macros.h"
//#include "vectra.h"
#include "string.h"
#include "vChronos.h"
#include <boost/thread.hpp>
#include <boost/date_time.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <termios.h>
#include <unistd.h>
#include <fcntl.h>
#include <fftw3.h>
#include <tiffio.h>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "fonctions.h"
#include "manip.h"
#include <complex>

#include "Point2D.h";

using namespace std;
using namespace cv;
typedef struct{
int a,b,c;
}coord3D;

typedef struct {
float x,y;
}float2D;

int _kbhit(void)//fction attente clavier->openCV?
{
  struct termios oldt, newt;
  int ch;
  int oldf;

  tcgetattr(STDIN_FILENO, &oldt);
  newt = oldt;
  newt.c_lflag &= ~(ICANON | ECHO);
  tcsetattr(STDIN_FILENO, TCSANOW, &newt);
  oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
  fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);

  ch = getchar();

  tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
  fcntl(STDIN_FILENO, F_SETFL, oldf);

  if(ch != EOF)
  {
    ungetc(ch, stdin);
    return 1;
  }

  return 0;
}

#include <PvDeviceFinderWnd.h>


#include <list>
typedef std::list<PvBuffer *> BufferList;


PV_INIT_SIGNAL_HANDLER();

#define BUFFER_COUNT ( 1 )

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
void AcquireImages( PvDevice *aDevice, PvStream *aStream, PvPipeline *aPipeline, coord3D cercle, int64_t dim, int64_t offset);
string IntToString ( int number );
//void TF2D(double *entree, double *fft, int taille);
//void circshift2D(double* entree, double* sortie, int taille, int decal);
void double2mat(vector<complex<double>> entree, Mat &sortie);
//void mat2double(Mat entree, double* sortie, int taille);
void mat2double(Mat entree,vector<double> &sortie);
//void SAV_Tiff(double *var_sav, char *chemin, int dim);
void normlog(double* entree, double* sortie, int taille);
void calcul_histo(Mat src, Mat dst, int hist_h, int hist_w);
//void decal2DCplxGen(nbCplx* entree, nbCplx* result, Var2D dim,Var2D decal);
//void TF2Dcplx(nbCplx *entree, nbCplx *fft, Var2D dim);
//void TF2Dcplx(double *entree, nbCplx *fft, double dim);
//float extract_val(string token,  string chemin_fic);
/* -------------------------------------------------------------------------- */
// Usage
/* -------------------------------------------------------------------------- */


static void
usage(int argc, char **argv)
{
        if ((argc - 1) == 0) {
                printf("Affichage de la TF pour hors-axe\n");
                printf("Usage: %s <paramètres obligatoires> <paramètres optionnels>\n", argv[0]);
                printf("Paramètres obligatoires: \n");
                printf("-r <rayon du cercle>:  \n");
                printf("-c <cx> <cy>: centre du cercle sur l'image 1024x1024 (orig°: top-left corner) \n");

                exit(EXIT_FAILURE);
        }

}
//
// Main function
//

int main(int argc, char *argv[])
{
//    usage(argc, argv);



       ///Initialiser les deux tableaux stockant les valeurs de recon.txt et config_manip.txt en mémoire

    manip m1;
    int cx,cy,rayon;

    cx = m1.circle_cx;
    cy = m1.circle_cy;
    rayon =m1.NXMAX; //attentino, il s'agit du rayon du disque dans l'hologramme.
//    while (argc > 0) {
//                if (!strcmp(argv[0], "-r") && (argc > 1)) {
//                        rayon = atoi(argv[1]);
//                        argc--;
//                        argv++;
//                }
//                if (!strcmp(argv[0], "-c") && (argc > 2)) {
//                        cx = atoi(argv[1]);
//                        cy = atoi(argv[2]);
//
//                        argc-=2;
//                        argv+=2;
//                }
//
//                argc--;
//                argv++;
//        }
    ///Initialiser les deux tableaux stockant les valeurs de recon.txt et config_manip.txt en mémoire
    string chemin_config_GUI,chemin_recon, chemin_config_manip,chemin_result,repertoire_config;
    string home=getenv("HOME");
    string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
    chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
    cout<<"lecture des chemins dans le fichier "<<chemin_config_GUI<<endl;

    repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);

    chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);
    chemin_recon=repertoire_config+"recon.txt";
    chemin_config_manip=repertoire_config+"config_manip.txt";
float tiptilt_factor_x = 0;//init rayon X en volt
float tiptilt_factor_y = 0;//init rayon  Y en volt
float2D VfOffset={-0, -0};
tiptilt_factor_x=(abs(extract_val("VXMIN", chemin_config_manip))+abs(extract_val("VXMAX", chemin_config_manip)))/20;///amplitude de tension/10=(|vmax|+|vmin|)/10
tiptilt_factor_y=(abs(extract_val("VYMIN", chemin_config_manip))+abs(extract_val("VYMAX", chemin_config_manip)))/20;
VfOffset.x=extract_val("VXMIN", chemin_config_manip)+tiptilt_factor_x*10;///offset en volt=amplitude tension+Vmin
VfOffset.y=extract_val("VYMIN", chemin_config_manip)+tiptilt_factor_y*10;
///----déclaration LabJack
  /*  Ljack LU3HV;
    Ljack_DAC ljDAC_phase(LU3HV);
    ASSERT(ljDAC_phase.connect(FIO7_6));

    Ljack_DAC ljDAC_flower(LU3HV);
    ASSERT(ljDAC_flower.connect(FIO5_4));


    ljDAC_flower.set_A_output(VfOffset.x);
    ljDAC_flower.set_B_output(VfOffset.y);
*/
///------------------------------------------
    coord3D Param_cercle={cx,cy,rayon};
    cout<<"###############################"<<endl;
    cout<<"# paramètres hors axe :       #"<<endl<<"#"<<" cx="<<cx<<", cy="<<cy<<", rayon : "<<rayon<<" #"<<endl;
    cout<<"# ou                          #"<<endl<<"#"<<"cx="<<cx<<", cy="<<-cy+m1.dimROI<<"(bas droit) #"<<endl;
    cout<<"###############################"<<endl;
    cout<<"##################################################"<<endl;
    cout<<"#Pressez 's' pour screenshot, 'Esc' pour quitter.#"<<endl;
    cout<<"##################################################"<<endl;
    const PvDeviceInfo *lDeviceInfo = NULL;
    PvDevice *lDevice = NULL;
    PvStream *lStream = NULL;

    PV_SAMPLE_INIT();

    PvDeviceFinderWnd *lDeviceFinderWnd = new PvDeviceFinderWnd();
    lDeviceInfo = SelectDevice( lDeviceFinderWnd );

///-----
    cout << "PvStreamSample:" << endl << endl;

    if ( NULL != lDeviceInfo )
    {
        lDevice = ConnectToDevice( lDeviceInfo );

        if ( lDevice != NULL )
        {
            ///-------------------------------------------------------changer la largeur de la ROI---------------------------------------

            PvGenParameterArray* lParameters = lDevice->GetParameters();
            // Get width parameter - mandatory GigE Vision parameter, it should be there.
            PvGenParameter *lParameter = lParameters->Get( "Width" );
            PvGenInteger *lWidthParameter = dynamic_cast<PvGenInteger *>( lParameter );
            //lStream = OpenStream( lDeviceInfo );
            if ( lWidthParameter == NULL )
            {
                cout << "Unable to get the width parameter." << endl;
            }
            int64_t lWidth=m1.dimROI;
            lWidthParameter->SetValue( lWidth );

            ///changer la hauteur de la ROI---------
            lParameter = lParameters->Get( "Height" );
            PvGenInteger *lHeightParameter = dynamic_cast<PvGenInteger *>( lParameter );
            if ( lHeightParameter == NULL )
            {
                cout << "Unable to get the width parameter." << endl;
            }
            int64_t lHeight=m1.dimROI;///;
            lHeightParameter->SetValue( lHeight );

            ///changer OffsetX de la ROI---------
            lParameter = lParameters->Get( "OffsetX" );
            PvGenInteger *lOffsetXParameter = dynamic_cast<PvGenInteger *>( lParameter );

            if ( lOffsetXParameter == NULL )
            {
                cout << "Unable to get the OffsetX parameter." << endl;
            }
            int64_t lOffsetX=m1.dimROI/2;
            lOffsetXParameter->SetValue( lOffsetX );

            ///changer OffsetY de la ROI---------
            lParameter = lParameters->Get( "OffsetY" );
            PvGenInteger *lOffsetYParameter = dynamic_cast<PvGenInteger *>( lParameter );
            if ( lOffsetYParameter == NULL )
            {
                cout << "Unable to get the OffsetY parameter." << endl;
            }
            int64_t lOffsetY=m1.dimROI/2;
            lOffsetYParameter->SetValue( lOffsetY );

///--------------------------------------------------------fin changement ROI------------------------------------------------

            lStream = OpenStream( lDeviceInfo );

            if ( lStream != NULL )
            {
                PvPipeline *lPipeline = NULL;

                ConfigureStream( lDevice, lStream );
                lPipeline = CreatePipeline( lDevice, lStream );
                if( lPipeline )
                {

                    AcquireImages( lDevice, lStream, lPipeline, Param_cercle,lHeight,lOffsetX );
                    delete lPipeline;
                }
                // Close the stream
                cout << "Closing stream" << endl;
                lStream->Close();
                PvStream::Free( lStream );
            }
         //-----Rétablir la valeur Max Largeur et hauteur de la caméra
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
        lDevice->Disconnect();
        PvDevice::Free( lDevice );
        }
    }


    cout << endl;
    cout << "<press a key to exit>" << endl;
    //PvWaitForKeyPress();

    delete lDeviceFinderWnd;
    lDeviceFinderWnd = NULL;

    PV_SAMPLE_TERMINATE();

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

PvDevice *ConnectToDevice( const PvDeviceInfo *aDeviceInfo )
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

void AcquireImages( PvDevice *aDevice, PvStream *aStream, PvPipeline *aPipeline, coord3D cercle, int64_t dim, int64_t offset)
{
    Mat image;
    Mat image_double(4*cercle.c,4*cercle.c,CV_64F);
        cv::Mat MatUborn(2*cercle.c,2*cercle.c,CV_64F);
    Point2D dim2DUBorn(2*cercle.c,2*cercle.c,2*cercle.c);
    Point2D dim2DROI(dim,dim,dim);
    Var2D dimROI={dim,dim}, decalROI={dim/2,dim/2};
    FFT_encaps tf2DROi(dim2DROI);
    FFT_encaps tf2DUborn(dim2DUBorn);
    Var2D coinCropHA={cercle.a-cercle.c,cercle.b-cercle.c};

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

    // Enable streaming and send the AcquisitionStart command
    cout << "Enabling streaming and sending AcquisitionStart command." << endl;
    aDevice->StreamEnable();
    lStart->Execute();

    char lDoodle[] = "|\\-|-/";
    int lDoodleIndex = 0;
    double lFrameRateVal = 0.0;
    double lBandwidthVal = 0.0;

    // Acquire images until the user instructs us to stop.
    cout << endl << "<press a key to stop streaming>" << endl;

    vChronos vTime("prise d'images"); vTime.clear();

    vTime.start();    //while ( !PvKbHit() )//tant que touche pas pressée

    unsigned char c;
    bool SAVE_IMAGE = false;

   // double* tmp= new double [dim*dim];
    vector<double> tmp(dim*dim);
   // nbCplx* fft= new nbCplx [dim*dim];
    vector <complex<double>> fft2D_UBorn(4*cercle.c*cercle.c);
    vector <complex<double>> fft2D_ROI(dim*dim);
    vector <complex <double>> UBorn(4*cercle.c*cercle.c);
//    cv::Mat MatUborn(2*cercle.c,2*cercle.c,CV_8UC1,Scalar( 0,0,0), cv::Mat::AUTO_STEP);

    //double* fft_shift= new double [dim*dim];
    //double* fft_log= new double [dim*dim];
//           int bool_wisdom3D=fftw_import_wisdom_from_filename("/home/mat/tomo_test/wisdom/test1024.wisdom");//charger ou calculer le fichier wisdom
//                if(bool_wisdom3D==0){
//                        cout<<"Calcul wisdom 1024 (~10 min)"<<endl;
//                        //prepare_wisdom(dimROI,"test1024.wisdom");
//                }

    //while(short int numImg=0;numImg<400;numImg++ )
    size_t continuer=1;
    /// Create Window for image
    //namedWindow("Image", CV_WINDOW_AUTOSIZE);
    /// Create Window for spectrum
    //namedWindow("TF", CV_WINDOW_AUTOSIZE);

    //while ( ! _kbhit() )
    while(continuer)
    {
        PvBuffer *lBuffer = NULL;
        PvResult lOperationResult;

        // Retrieve next buffer
        PvResult lResult = aPipeline->RetrieveNextBuffer( &lBuffer, 1000, &lOperationResult );

        if ( lResult.IsOK() )///Buffer acquis?
        {
            if ( lOperationResult.IsOK() ) ///l'opération finie sans erreur? (time out..)
            {
                PvPayloadType lType;
                //
                // We now have a valid buffer. This is where you would typically process the buffer.
                // -----------------------------------------------------------------------------------------
                // ...
                lFrameRate->GetValue( lFrameRateVal );
                lBandwidth->GetValue( lBandwidthVal );

                // If the buffer contains an image, display width and height.
                uint32_t lWidth = 0, lHeight = 0;
                lType = lBuffer->GetPayloadType();

                cout << fixed << setprecision( 1 );
                cout << lDoodle[ lDoodleIndex ];
                cout << " BlockID: " << uppercase << hex << setfill( '0' ) << setw( 16 ) << lBuffer->GetBlockID();

                if ( lType == PvPayloadTypeImage )///le buffer contient bien une image? Tout est ok, on traite le flux
                {
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
                    mat2double(lframe,tmp);
                    //------- write images to disk
                    //string racine="/ramdisk/ACQUIS/";

                    imwrite("/home/tomo/lframe.pgm", lframe);
                    //cvStartWindowThread();
                    /// Create Window
                    //char* source_window_img = "Image";
                    //namedWindow( source_window_img, CV_WINDOW_AUTOSIZE );
                    imshow("Image", lframe);

                    int hist_w = 512; int hist_h = 400;
                    Mat histImage( hist_h, hist_w, CV_8UC3, Scalar( 0,0,0) );
                    calcul_histo(lframe, histImage, hist_h, hist_w);
                    /// Display
                    namedWindow("Histogramme_plan_image", WINDOW_NORMAL);
                    imshow("Histogramme_plan_image", histImage);

                    holo2TF_UBorn(tmp,fft2D_UBorn,coinCropHA,tf2DROi);
                    TF2Dcplx_INV(fftshift2D(fft2D_UBorn),UBorn,tf2DUborn, 1);

                     double2mat(fftshift2D(UBorn),MatUborn);
                    // resize(MatUborn, image_double, cv::Size(8*cercle.c,8*cercle.c));
                   // normalize(image_double, image_double, 0, 1, cv::NORM_MINMAX);
                    normalize(MatUborn, MatUborn, 0, 1, cv::NORM_MINMAX);
                    namedWindow("Champ complexe", WINDOW_NORMAL);
                     imshow("Champ complexe",MatUborn);

                   // SAV_Tiff2D(UBorn,"Im", "/home/tomo/Uborn.tif", 1);
                   // TF2D(tmp, fft,tf2DROi,1);

                   // decal2DCplxGen(fft, fft_shift, decalROI);
                  //  normlog(fft_shift, fft_modlog, lHeight);
                   /* circshift2D(fft, fft_shift, lHeight, lHeight/2);
                    normlog(fft_shift, fft_log, lHeight);
                    //SAV_Tiff(fft_log, "/home/bailleul/image_fft.tif", lHeight);

                    Mat fft_mat(lHeight,lWidth,CV_64F);
                    double2mat(fft_log, fft_mat, lHeight);
                    //int r=cercle.rayon;
                    int y2=-cercle.b+dim;
                    circle(fft_mat, Point(cercle.a,cercle.b), cercle.c,(0,0,255),1,8,0);//tracer cercle hors-axe sur (image, coordonnée, rayon, couleur,epaisseur, typede ligne, ?))
                    circle(fft_mat, Point(dim/2,dim/2),2*cercle.c,(0,0,255),1,8,0);

                    circle(fft_mat, Point(cercle.a,y2), cercle.c,(0,0,255),1,8,0);//tracer cercle hors-axe sur (image, coordonnée, rayon, couleur,epaisseur, typede ligne, ?))

                    //ligne en x et en y, ordre 1 haut droit
                    line(fft_mat, Point(cercle.a-1.15*cercle.c,cercle.b),Point(cercle.a+1.15*cercle.c,cercle.b),(0,0,255),1,8,0);
                    line(fft_mat, Point(cercle.a,cercle.b-1.15*cercle.c),Point(cercle.a,cercle.b+1.15*cercle.c),(0,0,255),1,8,0);

                    //ligne en x et en y, ordre 1 bas droit
                    line(fft_mat, Point(cercle.a-1.15*cercle.c,y2),Point(cercle.a+1.15*cercle.c,y2),(0,0,255),1,8,0);
                    line(fft_mat, Point(cercle.a,y2-1.15*cercle.c),Point(cercle.a,y2+1.15*cercle.c),(0,0,255),1,8,0);
*/
//                    int dstWidth=lframe.cols * 2;
//                    int dstHeight=lframe.rows;
//                    Mat dst(dstHeight,dstWidth, CV_64F);
//                    Rect roi(Rect(0,0,lframe.cols, lframe.rows));
//                    Mat targetROI = dst(roi);
//                    lframe.copyTo(targetROI);
//                    targetROI = dst(Rect(lframe.cols,0,lframe.cols, lframe.rows));
//                    fft_mat=fft_mat/100;
//                    fft_mat.copyTo(targetROI);
//                    imshow("Image et sa TF", dst);

               //     imshow("TF", fft_mat/100 );

                   //circle(Mat& img, Point center, int radius, const Scalar& color, int thickness=1, int lineType=8, int shift=0)

                    c = cvWaitKey(10);
                    if (c == 's')
	                   SAVE_IMAGE = true;
                    else
	                   SAVE_IMAGE = false;

                    if (SAVE_IMAGE)
                    {
	                   // imwrite("last_fmod.pgm", fft_mat);
	                    imwrite("image_show_fourier.pgm",lframe);
	                }
	                if(c==27)//27=code ascii Escape
                    {
                        continuer=0;
                    }


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
    //destroyAllWindows ();
    destroyWindow("Image");
    waitKey(1);
    destroyWindow("TF");
    waitKey(1);
    destroyWindow("Histogramme_plan_image");
    waitKey(1);
    cout<<endl<<"sortie boucle d'aquisition"<<endl;




    PvGetChar(); // Flush key buffer for next stop.
    cout << endl << endl;

    // Tell the device to stop sending images.
    cout << "Sending AcquisitionStop command to the device" << endl;
    lStop->Execute();

    // Disable streaming on the device
    cout << "Disable streaming on the controller." << endl;
    aDevice->StreamDisable();

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

void mat2double(Mat entree,vector<double> &sortie)
{
    int taille=sqrt(sortie.size());
    for(int y=0; y<taille; y++)
    {
        for(int x=0; x<taille; x++)
            {
                int cpt=y*taille+x;
                sortie[cpt]=entree.at<unsigned char>(y,x);
            }
    }
}

void double2mat(vector<complex<double>> entree, Mat& sortie)
{
    int taille=sqrt(entree.size());
    for(int y=0; y<taille; y++)
    {
        for(int x=0; x<taille; x++)
            {
                int cpt=y*taille+x;
                sortie.at<double>(y,x)=(entree[cpt].imag());
            }
    }
}

void SAV_Tiff(double *var_sav, char *chemin, int dim)
{
    TIFF *tif= TIFFOpen(chemin, "wb");
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

void circshift2D(double* entree, double* sortie, int taille, int decal)
{
    //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
    decal=decal%taille;

    for(int yi=0; yi<decal; yi++) {
        for(int xi=0; xi<decal; xi++)
        {
            int pixel=yi*taille+xi;
            int pixel_shift=(yi+decal)*taille+xi+decal;

            //1er quadrant vers 4 eme
            sortie[pixel_shift]=entree[pixel];
            //4 eme quadrant vers 1er
            sortie[pixel]=entree[pixel_shift];
            //2eme vers 3eme
            sortie[(yi+decal)*taille+xi]=entree[pixel+decal];
            //3eme vers 2eme
            sortie[pixel+decal]=entree[(yi+decal)*taille+xi];
        }
    }
}

void normlog(double* entree, double* sortie, int taille)
{
    for (size_t j = 0; j < taille; j++)
    {
        for (size_t i = 0; i < taille; i++)
        {
            int cpt = j* taille + i;
            sortie[cpt] = log(1+sqrt(entree[cpt]*entree[cpt]));
        }
    }

    double max_fft_log = 0;
    int cpt_max = 0;

    for (size_t j = 0; j < taille; j++)
    {
        for (size_t i = 0; i < taille; i++)
        {
            int cpt = j* taille + i;
            if(sortie[cpt] > sortie[cpt_max])
            {
                cpt_max = cpt;
            }
        }
    }
    max_fft_log = sortie[cpt_max];

    for (size_t j = 0; j < taille; j++)
    {
        for (size_t i = 0; i < taille; i++)
        {
            int cpt = j* taille + i;
            sortie[cpt] = sortie[cpt] / max_fft_log * 100 * 255;
        }
    }

}
void calcul_histo(Mat src, Mat dst, int hist_h, int hist_w)
{

    /// Establish the number of bins
    int histSize = 256;

    float range[] = { 0, 256 } ;
    const float* histRange = { range };

    bool uniform = true; bool accumulate = false;

    Mat r_hist;

    /// Compute the histograms:
    calcHist( &src, 1, 0, Mat(), r_hist, 1, &histSize, &histRange, uniform, accumulate );

    int bin_w = cvRound( (double) hist_w/histSize );

    /// Normalize the result to [ 0, dst.rows ]
    normalize(r_hist, r_hist, 0, dst.rows, NORM_MINMAX, -1, Mat() );

    /// Draw for each channel
    for( int i = 1; i < histSize; i++ )
    {
        line( dst, Point( bin_w*(i-1), hist_h - cvRound(r_hist.at<float>(i-1)) ) ,
                         Point( bin_w*(i), hist_h - cvRound(r_hist.at<float>(i)) ),
                         Scalar( 0, 0, 255), 2, 8, 0  );
    }

}
/*
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
*/
/*
void decal2DCplxGen(nbCplx* entree, nbCplx* result, Var2D dim,Var2D decal)
{
        //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
        decal.y=fmod(decal.y,dim.y);
        decal.x=fmod(decal.x,dim.x);
        if(decal.x<0)
        decal.x=dim.x+decal.x;
        if(decal.y<0)
        decal.y=dim.y+decal.y;
       for(size_t yi=0; yi<dim.y-decal.y; yi++) {
                for(size_t xi=0; xi<dim.x-decal.x; xi++)
                { //cout<<"xi,yi="<<xi<<","<<yi<<endl;
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
                        result[pixel_shift]=entree[pixel];
                }
                for(size_t xi=dim.x-decal.x; xi<dim.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
                        result[pixel_shift]=entree[pixel];
                }
        }
              for(size_t yi=dim.y-decal.y; yi<dim.y; yi++) {
                for(size_t xi=0; xi<dim.x-decal.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(-dim.y+yi+decal.y)*dim.x+xi+decal.x;
                        result[pixel_shift]=entree[pixel];
                }
                for(size_t xi=dim.x-decal.x; xi<dim.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(-dim.y+yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
                        result[pixel_shift]=entree[pixel];
                }
        }
}
*/
/*
void TF2Dcplx(double *entree, nbCplx *fft, Var2D dim)
{
        int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        //int nthreads=3;
        fftw_plan_with_nthreads(4);
        int N=dim.x*dim.y;
        fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
        fftw_plan p;
        //Réservation memoire
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        //Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<N; cpt++) {
                in[cpt][0]=entree[cpt];
                in[cpt][1]=0;
        }
        //calcul du plan, parametre servant a calculer et optimiser le FFT
        p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
          // p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_FORWARD,FFTW_USE_WISDOM);
        fftw_execute(p);

        for(int cpt=0; cpt<(N); cpt++) {
                fft[cpt].Re=out[cpt][0]/N; //division par N (dim*dim) pour normaliser la fftw qui n'est pas normalisée
                fft[cpt].Im=out[cpt][1]/N;
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}

void TF2Dcplx(nbCplx *entree, nbCplx *fft, Var2D dim)
{
        int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        //int nthreads=3;
        fftw_plan_with_nthreads(4);
        int N=dim.x*dim.y;
        fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
        fftw_plan p;
        //Réservation memoire
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        //Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<N; cpt++) {
                in[cpt][0]=entree[cpt].Re;
                in[cpt][1]=entree[cpt].Im;
        }
        //calcul du plan, parametre servant a calculer et optimiser le FFT
        p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
          // p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_FORWARD,FFTW_USE_WISDOM);
        fftw_execute(p);

        for(int cpt=0; cpt<(N); cpt++) {
                fft[cpt].Re=out[cpt][0]/N; //division par N (dim*dim) pour normaliser la fftw qui n'est pas normalisée
                fft[cpt].Im=out[cpt][1]/N;
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}
*/
