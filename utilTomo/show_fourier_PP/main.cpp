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
#include <PvSampleUtils.h>
#include <PvDevice.h>
#include <PvDeviceGEV.h>
#include <PvDeviceU3V.h>
#include <PvStream.h>
#include <PvStreamGEV.h>
#include <PvStreamU3V.h>
#include <PvPipeline.h>
#include <PvBuffer.h>


//#include "macros.h"
//#include "vectra.h"
#include "string.h"
//#include "vChronos.h"
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
void TF2D(double *entree, double *fft, int taille);
void circshift2D(double* entree, double* sortie, int taille, int decal);
void double2mat(double* entree, Mat sortie, int taille);
void mat2double(Mat entree, double* sortie, int taille);
void SAV_Tiff(double *var_sav, char *chemin, int dim);
void normlog(double* entree, double* sortie, int taille);
void calcul_histo(Mat src, Mat dst, int hist_h, int hist_w);
bool SetExposureTime( PvDevice *pDevice, float iExposureTime );
int SetExposureMode( PvDevice *pDevice, int iExposureMode );
//float GetExposureTime( PvDevice *pDevice,float iExposureTime );
double GetExposureTime( PvDevice *pDevice);
bool AddExposureTime( PvDevice *pDevice,float DeltaExposureTime );
bool SubExposureTime( PvDevice *pDevice,float DeltaExposureTime );
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
    int f1=m1.dimROI/2-m1.coord_porteuse;//coordonnée porteuse
    int f2=m1.dimROI/2+m1.coord_porteuse;

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



    coord3D Param_cercle={cx,cy,rayon};

     cout<<"###############################"<<endl;
    cout<<"# paramètres hors axe théoriques : "<<endl<<"#"<<"rayon : "<<rayon<<" #"<<endl;
    cout<<"# paramètres hors axe théoriques : (cx1,cy1)=("<<f1<<","<<f1<<")"<<endl;
    cout<<"# paramètres hors axe théoriques : (cx2,cy2)=("<<f1<<","<<f2<<")"<<endl;
    cout<<"# paramètres hors axe théoriques : (cx3,cy4)=("<<f2<<","<<f1<<")"<<endl;
    cout<<"# paramètres hors axe théoriques : (cx4,cy4)=("<<f2<<","<<f2<<")"<<endl;
    cout<<"###############################"<<endl;
    cout<<"# paramètres hors axe utilisés :       #"<<endl<<"#"<<" cx="<<cx<<", cy="<<cy<<", rayon : "<<rayon<<" #"<<endl;
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

//    vChronos vTime("prise d'images"); vTime.clear();

  //  vTime.start();    //while ( !PvKbHit() )//tant que touche pas pressée

    unsigned char c;
    bool SAVE_IMAGE = false;

    double* tmp= new double [dim*dim];
    double* fft= new double [dim*dim];
    double* fft_shift= new double [dim*dim];
    double* fft_log= new double [dim*dim];
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
    string strTpsExpo;
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

            //    cout << fixed << setprecision( 1 );
            //    cout << lDoodle[ lDoodleIndex ];
           //     cout << " BlockID: " << uppercase << hex << setfill( '0' ) << setw( 16 ) << lBuffer->GetBlockID();

                if ( lType == PvPayloadTypeImage )///le buffer contient bien une image? Tout est ok, on traite le flux
                {
                    // Get image specific buffer interface.
                    PvImage *lImage = lBuffer->GetImage();

                    // Get image data pointer so we can pass it to CV::MAT container
                    unsigned char *img = lImage->GetDataPointer();

                    // Read width, height.
                    lWidth = lImage->GetWidth();
                    lHeight = lImage->GetHeight();
               //     cout << "  W: " << dec << lWidth << " H: " << lHeight;

                    // Copy/convert Pleora Vision image pointer to cv::Mat container
                    cv::Mat lframe(lHeight,lWidth,CV_8UC1,img, cv::Mat::AUTO_STEP);
                    //------- write images to disk
                    //string racine="/ramdisk/ACQUIS/";

                    //imwrite("/ramdisk/ACQUIS/.pgm", lframe);
                    //cvStartWindowThread();
                    /// Create Window
                    char* source_window_img = "Image";
                   // namedWindow( source_window_img, CV_WINDOW_AUTOSIZE );

                    imshow("Image", lframe);

                    int hist_w = 512; int hist_h = 400;
                    Mat histImage( hist_h, hist_w, CV_8UC3, Scalar( 0,0,0) );
                    calcul_histo(lframe, histImage, hist_h, hist_w);
                    /// Display
                    namedWindow("Histogramme_plan_image", WINDOW_NORMAL);
                    imshow("Histogramme_plan_image", histImage);

                    mat2double(lframe, tmp, lHeight);
                    //SAV_Tiff(tmp, "/home/bailleul/image.tif", lHeight);
                    TF2D(tmp, fft, lHeight);
                    circshift2D(fft, fft_shift, lHeight, lHeight/2);
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
    strTpsExpo=to_string(GetExposureTime(aDevice));

        putText(fft_mat, strTpsExpo,Point(50,50),FONT_HERSHEY_COMPLEX,0.5,Scalar(255),1);
                    imshow("TF", fft_mat/100 );

                   //circle(Mat& img, Point center, int radius, const Scalar& color, int thickness=1, int lineType=8, int shift=0)

                    c = cv::waitKey(10);
                    if (c == 's')
	                   SAVE_IMAGE = true;
                    else
	                   SAVE_IMAGE = false;

                    if (SAVE_IMAGE)
                    {
	                    imwrite("last_fmod.pgm", fft_mat);
	                    imwrite("image_show_fourier.pgm",lframe);
	                }
	                if(c==27)//27=code ascii Escape
                    {
                        continuer=0;
                    }
                     if(c==100)//112=code ascii P
                    {
                        AddExposureTime( aDevice,20 );
                        cout<<GetExposureTime(aDevice)<<endl;

                    }
                         if(c==113)//112=code ascii P
                    {
                        SubExposureTime( aDevice,20 );
                    }



                }
                else
                {
                   // cout << " (buffer does not contain image)";
                }
             //   cout << "  " << lFrameRateVal << " FPS  " << ( lBandwidthVal / 1000000.0 ) << " Mb/s   \r";
            }
            else
            {
                // Non OK operational result
            //    cout << lDoodle[ lDoodleIndex ] << " " << lOperationResult.GetCodeString().GetAscii() << "\r";
            }

            // Release the buffer back to the pipeline
            aPipeline->ReleaseBuffer( lBuffer );
        }
        else
        {
            // Retrieve buffer failure
          //  cout << lDoodle[ lDoodleIndex ] << " " << lResult.GetCodeString().GetAscii() << "\r";
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



    delete [] tmp;
    delete [] fft;
    delete [] fft_shift;
    delete [] fft_log;

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

void TF2D(double *entree, double *fft, int taille)
{
        int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        //int nthreads=3;
        fftw_plan_with_nthreads(4);
        int N=taille*taille;
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
        p=fftw_plan_dft_2d(taille, taille, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p); /* repeat as needed */

        for(int cpt=0; cpt<(N); cpt++) {
                fft[cpt]=out[cpt][0]/N; //division par N^2 pour normaliser la fftw qui n'est pas normalisée
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}

void mat2double(Mat entree, double* sortie, int taille)
{
    for(int y=0; y<taille; y++)
    {
        for(int x=0; x<taille; x++)
            {
                int cpt=y*taille+x;
                sortie[cpt]=entree.at<unsigned char>(y,x);
            }
    }
}

void double2mat(double* entree, Mat sortie, int taille)
{
    for(int y=0; y<taille; y++)
    {
        for(int x=0; x<taille; x++)
            {
                int cpt=y*taille+x;
                sortie.at<double>(y,x)=entree[cpt];
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

// Sets the exposure mode
//iExposureMode: 0:Timed 1:TriggerWidth
int SetExposureMode( PvDevice *pDevice, int iExposureMode )
{
PvGenEnum* lPvGenEnum=dynamic_cast<PvGenEnum*>(pDevice->GetCommunicationParameters()->Get( "ExposureMode" ) );
if( lPvGenEnum==NULL ) return false;
PvResult lResult = lPvGenEnum->SetValue(iExposureMode);
return (bool)lResult.IsOK();
}


//iExposureTime: Min:0 Max:16777215
bool SetExposureTime( PvDevice *pDevice,float iExposureTime )
{
double value;
//pDevice->StreamDisable();
PvGenFloat *ic_exposure= dynamic_cast<PvGenFloat*>(pDevice ->GetParameters()->Get("ExposureTime"));
PvString lValue, lName, lCategory;
ic_exposure->GetName( lName );
ic_exposure->ToString( lValue );
//cout << lName.GetAscii() << ": " << lValue.GetAscii() << endl;
double TempsExpo=0;
ic_exposure->GetValue(TempsExpo);


cout<< "--------------------------------"<<endl;
if( ic_exposure==NULL ) return true;

    PvResult lResult =ic_exposure->SetValue(iExposureTime);

return lResult.IsOK();
}

//iExposureTime: Min:0 Max:16777215
double  GetExposureTime( PvDevice *pDevice)
{
double value;
//pDevice->StreamDisable();
PvGenFloat *ic_exposure= dynamic_cast<PvGenFloat*>(pDevice ->GetParameters()->Get("ExposureTime"));
PvString lValue, lName, lCategory;
ic_exposure->GetName( lName );
ic_exposure->ToString( lValue );
//cout << lName.GetAscii() << ": " << lValue.GetAscii() << endl;
double TempsExpo=0;
ic_exposure->GetValue(TempsExpo);
//cout<<"tempsExpo="<<TempsExpo<<endl;

//return lResult.IsOK();
return TempsExpo;
}

//iExposureTime: Min:0 Max:16777215
bool AddExposureTime( PvDevice *pDevice,float DeltaExposureTime )
{
double value;
//pDevice->StreamDisable();
PvGenFloat *ic_exposure= dynamic_cast<PvGenFloat*>(pDevice ->GetParameters()->Get("ExposureTime"));
PvString lValue, lName, lCategory;
ic_exposure->GetName( lName );
ic_exposure->ToString( lValue );
cout << lName.GetAscii() << ": " << lValue.GetAscii() << endl;
double TempsExpo=0;
ic_exposure->GetValue(TempsExpo);


cout<< "--------------------------------"<<endl;
if( ic_exposure==NULL ) return true;

    PvResult lResult =ic_exposure->SetValue(TempsExpo+DeltaExposureTime);

return lResult.IsOK();
}
//iExposureTime: Min:0 Max:16777215
bool SubExposureTime( PvDevice *pDevice,float DeltaExposureTime )
{
double value;
//pDevice->StreamDisable();
PvGenFloat *ic_exposure= dynamic_cast<PvGenFloat*>(pDevice ->GetParameters()->Get("ExposureTime"));
PvString lValue, lName, lCategory;
ic_exposure->GetName( lName );
ic_exposure->ToString( lValue );
cout << lName.GetAscii() << ": " << lValue.GetAscii() << endl;
double TempsExpo=0;
ic_exposure->GetValue(TempsExpo);


cout<< "--------------------------------"<<endl;
if( ic_exposure==NULL ) return true;

    PvResult lResult =ic_exposure->SetValue(TempsExpo-DeltaExposureTime);

return lResult.IsOK();
}
//iExposureTime: Min:0 Max:16777215
float GetExposureTime( PvDevice *pDevice,float iExposureTime )
{
double value;
//pDevice->StreamDisable();
PvGenFloat *ic_exposure= dynamic_cast<PvGenFloat*>(pDevice ->GetParameters()->Get("ExposureTime"));
PvString lValue, lName, lCategory;
ic_exposure->GetName( lName );
ic_exposure->ToString( lValue );

cout << lName.GetAscii() << ": " << lValue.GetAscii() << endl;

// PvResult lResult =ic_exposure->GetValue();
return 0;
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
