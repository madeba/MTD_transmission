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

#include "macros.h"
#include "vectra.h"
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
typedef enum PRECISION {CHAR, UINT, INT, FLOAT, DOUBLE};
using namespace std;
using namespace cv;
typedef struct{
int a,b,c;
}coord3D;
typedef struct {
  double Re,Im;
}nbCplx;
typedef struct {
float x,y;
}float2D;
typedef struct {
float x,y;
}Var2D;


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
void AcquireImages( PvDevice *aDevice, PvStream *aStream, PvPipeline *aPipeline, coord3D cercle );
string IntToString ( int number );
void TF2D(double *entree, double *fft, int taille);
void TF2Dcplx(nbCplx *entree, nbCplx *fft, Var2D dim);
void circshift2D(double* entree, double* sortie, int taille, int decal);
void double2mat(double* entree, Mat sortie, int taille);
void mat2double(Mat entree, double* sortie, int taille);
void SAV_Tiff(double *var_sav, char *chemin, int dim);
void normlog(nbCplx* entree, double* sortie, int taille);
void coupe_double(double *src, double *dest, Var2D dim_src, Var2D dim_dest, Var2D coin);
void coupe_cplx(nbCplx *src, nbCplx *dest, Var2D dim_src, Var2D dim_dest, Var2D coin);
void mat2cplx(Mat entree, nbCplx* sortie, int taille);
void cplxMod2mat(nbCplx* entree, Mat sortie, int taille);
void SAV_Tiff_Re(nbCplx *var_sav, char *chemin, int dim);
void decal2DCplxGen(nbCplx* entree, nbCplx* result, Var2D dim,Var2D decal);
void SAV_Re(nbCplx *var_sav, size_t NbPix2D, char *chemin, enum PRECISION precision, char options[]);
void SAV(double *var_sav, int NbPix2D, char *chemin, PRECISION precision, char options[]);
/* -------------------------------------------------------------------------- */
// Usage
/* -------------------------------------------------------------------------- */


static void
usage(int argc, char **argv)
{
        if ((argc - 1) == 0) {
                printf("Afficahe de la TF pour hors-axe\n");
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
    usage(argc, argv);
    int cx=770,cy=220,rayon=110;
    while (argc > 0) {
                if (!strcmp(argv[0], "-r") && (argc > 1)) {
                        rayon = atoi(argv[1]);
                        argc--;
                        argv++;
                }
                if (!strcmp(argv[0], "-c") && (argc > 2)) {
                        cx = atoi(argv[1]);
                        cy = atoi(argv[2]);

                        argc-=2;
                        argv+=2;
                }

                argc--;
                argv++;
        }



    coord3D Param_cercle={cx,cy,rayon};
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
            //changer la largeur de la ROI---------
            PvGenParameterArray* lParameters = lDevice->GetParameters();
            // Get width parameter - mandatory GigE Vision parameter, it should be there.
            PvGenParameter *lParameter = lParameters->Get( "Width" );
            PvGenInteger *lWidthParameter = dynamic_cast<PvGenInteger *>( lParameter );
            //lStream = OpenStream( lDeviceInfo );
            if ( lWidthParameter == NULL )
            {
                cout << "Unable to get the width parameter." << endl;
            }

            int64_t lWidth=1024;
            lWidthParameter->SetValue( lWidth );

              //changer la hauteur de la ROI---------

            // Get ROI width parameter -
            lParameter = lParameters->Get( "Height" );
            PvGenInteger *lHeightParameter = dynamic_cast<PvGenInteger *>( lParameter );

            if ( lHeightParameter == NULL )
            {
                cout << "Unable to get the width parameter." << endl;
            }
            int64_t lHeight=1024;
            lHeightParameter->SetValue( lHeight );
            cout<<"coucou height"<<endl;
            //-------fin changement ROI

            lStream = OpenStream( lDeviceInfo );

            if ( lStream != NULL )
            {
                PvPipeline *lPipeline = NULL;

                ConfigureStream( lDevice, lStream );
                lPipeline = CreatePipeline( lDevice, lStream );
                if( lPipeline )
                {

                    AcquireImages( lDevice, lStream, lPipeline, Param_cercle );
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

void AcquireImages( PvDevice *aDevice, PvStream *aStream, PvPipeline *aPipeline, coord3D cercle)
{
    Mat image;
    Var2D
    dim_holo={2*cercle.c,2*cercle.c},
    coin={512,512};//{round(cercle.a-cercle.c/2),round(cercle.b-cercle.c/2)};
    Var2D dim_src={1024,1024},dim_src_decal={round(dim_src.x/2),round(dim_src.y/2)};
    size_t NbPix2D=dim_src.x*dim_src.y, NbPixMini=2*cercle.c*2*cercle.c;
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

    nbCplx* tmp= new nbCplx [NbPix2D];
    nbCplx* tmp_mini=new nbCplx [NbPixMini];//image 2D recalée
    nbCplx* fft_mini=new nbCplx [NbPixMini];//TF image 2D recalée
    nbCplx* fft= new nbCplx [NbPix2D];
    nbCplx* fft_shift= new nbCplx [NbPix2D];
    double* fft_modlog= new double [NbPix2D];
    double* fft_modlog_mini= new double [NbPix2D];
    double* tmp_mini_log=new double [NbPix2D];


    int bool_wisdom3D=fftw_import_wisdom_from_filename("/home/bailleul/Documents/show_holo_PP/test1024.wisdom");//charger ou calculer le fichier wisdom
                if(bool_wisdom3D==0){
                        cout<<"Calcul wisdom 1024 (~10 min)"<<endl;
                        //prepare_wisdom(dimROI,"test1024.wisdom");
                }

    //while(short int numImg=0;numImg<400;numImg++ )
    size_t continuer=1;
    /// Create Window for image
    namedWindow("Image", CV_WINDOW_AUTOSIZE);
    /// Create Window for spectrum
    namedWindow("TF", CV_WINDOW_AUTOSIZE);

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
                    //------- write images to disk
                    //string racine="/ramdisk/ACQUIS/";

                    //imwrite("/ramdisk/ACQUIS/.pgm", lframe);
                    //cvStartWindowThread();
                    /// Create Window
                    //char* source_window_img = "Image";
                    //namedWindow( source_window_img, CV_WINDOW_AUTOSIZE );
                    imshow("Image", lframe );

                    //mat2double(lframe, tmp, lHeight);
                    mat2cplx(lframe, tmp, lHeight);


                    //SAV_Tiff_Re(tmp, "home/mat/tomo_test/truc.pgm", lHeight);
                    //SAV_Tiff(tmp, "/home/bailleul/image.tif", lHeight);
                    TF2Dcplx(tmp, fft, dim_src);

                    decal2DCplxGen(fft, fft_shift, dim_src, dim_src_decal);
                    normlog(fft_shift, fft_modlog, lHeight);
                    //SAV(fft_modlog,NbPix2D,"/home/mat/tomo_test/fft_mod_log.raw",FLOAT,"wb");
                   /* //SAV_Tiff(fft_log, "/home/bailleul/image_fft.tif", lHeight);

                    Mat fft_mat(lHeight,lWidth,CV_64F);
                    double2mat(fft_modlog, fft_mat, lHeight);
                    // int r=110;
                    circle(fft_mat, Point(cercle.a,cercle.b), cercle.c,(0,0,255),1,8,0);

                    imshow("TF", fft_mat/100 );*/
                    coupe_cplx(fft_shift, fft_mini, dim_src, dim_holo, coin);

                    normlog(fft_mini, fft_modlog_mini, 2*cercle.c);
                    //coupe_double(fft_modlog,fft_modlog_mini,dim_src, dim_holo, coin);

                    //SAV(fft_modlog_mini,NbPixMini,"/home/mat/tomo_test/fft_modlog_mini.raw",FLOAT,"wb");
                    //SAV_Tiff(fft_modlog, "/home/mat/tomo_test/image_fft_crop.tif", 1024);
                    TF2Dcplx(fft_mini, tmp_mini, dim_holo);
                    normlog(tmp_mini, tmp_mini_log, 2*cercle.c);
                    //SAV(tmp_mini_log,NbPix2D,"/home/mat/tomo_test/img_mod.raw",FLOAT,"wb");

                    int taille=2*cercle.c;
                    Mat demodule(taille, taille, CV_64F);
                    double2mat(tmp_mini_log, demodule, taille);

                    /// Create Window
                    char* source_window_img = "Demodule";
                    namedWindow( source_window_img, CV_WINDOW_AUTOSIZE );
                    imshow("Demodule", demodule);

                    c = cvWaitKey(10);
                    if (c == 's')
	                   SAVE_IMAGE = true;

                    else
	                   SAVE_IMAGE = false;

                    if (SAVE_IMAGE)
                    {
	                    //imwrite("last_fmod.pgm", fft_mat);
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
            // Retrieve buffer failuredim
            cout << lDoodle[ lDoodleIndex ] << " " << lResult.GetCodeString().GetAscii() << "\r";
        }

        ++lDoodleIndex %= 6;
    }
    //destroyAllWindows ();
    destroyWindow("Image");
    waitKey(1);
    destroyWindow("TF");
    waitKey(1);


    cout<<endl<<"sortie boucle d'aquisition"<<endl;

    delete [] tmp;
    delete [] fft;
    delete [] tmp_mini;

    delete [] fft_mini;
    delete [] fft_shift;
    delete [] fft_modlog;
     delete [] fft_modlog_mini;
       delete [] tmp_mini_log;

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

/*void SAV(double *var_sav, int NbPix2D, char *chemin, enum PRECISION precision, char options[])
{
        FILE *fichier_ID;
        fichier_ID= fopen(chemin, options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision) {
        case DOUBLE: //64 bit

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        double tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case FLOAT://32 bits float

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        float tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;

        case INT: //32 bit signé

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case UINT://32 bit non signé

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        unsigned int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case CHAR: //8 bits

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        char tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        default:
                break;
        }

        fclose(fichier_ID);
}*/

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

void mat2cplx(Mat entree, nbCplx* sortie, int taille)
{
    for(int y=0; y<taille; y++)
    {
        for(int x=0; x<taille; x++)
            {
                int cpt=y*taille+x;
                sortie[cpt].Re=entree.at<unsigned char>(y,x);

                sortie[cpt].Im=0;
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
void cplxMod2mat(nbCplx * entree, Mat sortie, int taille)
{
    for(int y=0; y<taille; y++)
    {
        for(int x=0; x<taille; x++)
            {
                int cpt=y*taille+x;
                sortie.at<double>(y,x)=entree[cpt].Re*entree[cpt].Re+entree[cpt].Im*entree[cpt].Im;
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

void normlog(nbCplx* entree, double* sortie, int taille)
{
    for (size_t j = 0; j < taille; j++)
    {
        for (size_t i = 0; i < taille; i++)
        {
            int cpt = j* taille + i;
            sortie[cpt] = log(1+sqrt(entree[cpt].Re*entree[cpt].Re+entree[cpt].Im*entree[cpt].Im));
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


///découpe une une fenetre de dimension dim_dest, coin haut gauche coin, dans src de taille dim_src.
void coupe_double(double *src, double *dest, Var2D dim_src, Var2D dim_dest, Var2D coin)
{
    size_t nbPixSrc=dim_src.x*dim_src.y;
    size_t cpt_destX,cpt_destY, cpt_dest1D,
            cpt_srcX,cpt_srcY, cpt_src1D;

    for(cpt_destX=0; cpt_destX<dim_dest.x; cpt_destX++){
        for(cpt_destY=0; cpt_destY<dim_dest.y; cpt_destY++){

            cpt_dest1D=cpt_destX+cpt_destY*dim_dest.x;///coord 1D destination

            cpt_srcX=coin.x+cpt_destX;///coord X src
            cpt_srcY=coin.y+cpt_destY;///coord Y src
            cpt_src1D=cpt_srcX+cpt_srcY*dim_src.x;///coord 1D source

            dest[cpt_dest1D]=src[cpt_src1D];

        }

    }

}

///découpe une une fenetre de dimension dim_dest, coin haut gauche coin, dans src de taille dim_src.
void coupe_cplx(nbCplx *src, nbCplx *dest, Var2D dim_src, Var2D dim_dest, Var2D coin)
{
    size_t nbPixSrc=dim_src.x*dim_src.y;
    size_t cpt_destX,cpt_destY, cpt_dest1D,
            cpt_srcX,cpt_srcY, cpt_src1D;

    for(cpt_destX=0; cpt_destX<dim_dest.x; cpt_destX++){
        for(cpt_destY=0; cpt_destY<dim_dest.y; cpt_destY++){

            cpt_dest1D=cpt_destX+cpt_destY*dim_dest.x;///coord 1D destination


            cpt_srcX=coin.x+cpt_destX;///coord X src
            cpt_srcY=coin.y+cpt_destY;///coord Y src
            cpt_src1D=cpt_srcX+cpt_srcY*dim_src.x;///coord 1D source

            dest[cpt_dest1D].Re=src[cpt_src1D].Re;
            dest[cpt_dest1D].Im=src[cpt_src1D].Im;

        }

    }

}

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
        fftw_execute(p); /* repeat as needed */

        for(int cpt=0; cpt<(N); cpt++) {
                fft[cpt].Re=out[cpt][0]/N; //division par N (dim*dim) pour normaliser la fftw qui n'est pas normalisée
                fft[cpt].Im=out[cpt][1]/N;
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}
void SAV_Re(nbCplx *var_sav, size_t NbPix2D, char *chemin, enum PRECISION precision, char options[])
{       size_t taille=NbPix2D;
        FILE *fichier_ID;
        fichier_ID= fopen(chemin, options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision) {
        case DOUBLE: //64 bit

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        double tampon=var_sav[cpt].Re;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case FLOAT://32 bits float

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        float tampon=var_sav[cpt].Re;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;

        case INT: //32 bit signé

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        int tampon=var_sav[cpt].Re;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case UINT://32 bit non signé

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        unsigned int tampon=var_sav[cpt].Re;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case CHAR: //8 bits

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        char tampon=var_sav[cpt].Re;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        default:
                break;
        }

        fclose(fichier_ID);
}
void SAV_Tiff_Re(nbCplx *var_sav, char *chemin, int dim)
{
    TIFF *tif= TIFFOpen(chemin, "a");
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
            strip_buf[col]=(float)var_sav[cpt].Re;
        }
        TIFFWriteEncodedStrip (tif, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif);
    TIFFClose(tif);
}

///####################fonction################"
void SAV(double *var_sav, int NbPix2D, char *chemin, PRECISION precision, char options[])
{
        FILE *fichier_ID;
        fichier_ID= fopen(chemin, options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision) {
        case DOUBLE: //64 bit

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        double tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case FLOAT://32 bits float

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        float tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;

        case INT: //32 bit signé

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case UINT://32 bit non signé

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        unsigned int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case CHAR: //8 bits

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        char tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        default:
                break;
        }

        fclose(fichier_ID);
}
