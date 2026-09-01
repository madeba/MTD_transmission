#include "zernike.h"
#include <cmath>
#include <iostream>

using namespace std;
///constructeur Zernike
Zernike::Zernike(int width,
                 int height,
                 double radius,
                 double cx,
                 double cy)
{
    W=width;
    H=height;
    R=radius;
    CX=cx;
    CY=cy;
}

double Zernike::radial(int n,
                       int m,
                       double rho)
{
    double result=0;

    int mabs=abs(m);

    for(int k=0;
        k<=(n-mabs)/2;
        k++)
    {
        double num=
            pow(-1,k)*tgamma(n-k+1);

        double den=
            tgamma(k+1)*
            tgamma((n+mabs)/2-k+1)*
            tgamma((n-mabs)/2-k+1);

        result +=
            (num/den)*
            pow(rho,n-2*k);
    }

    return result;
}
///évaluer un polynome de Zernike
double Zernike::evaluateZernike(int n,
                        int m,
                        double rho,
                        double theta)
{
    if(rho>1.0)
        return 0.0;

    double Rnm=radial(n,m,rho);

    if(m>0)
        return Rnm*cos(m*theta);

    else if(m<0)
        return Rnm*sin(-m*theta);

    else
        return Rnm;
}

vector<pair<int,int>>
Zernike::generateModes(int maxMode)
{
    vector<pair<int,int>> modes;//les modes sont indicés par une paire d'entiers

    int count=0;

    for(int n=0;n<20;n++)
    {
        for(int m=-n;m<=n;m+=2)
        {
            modes.push_back({n,m});

            count++;

            if(count>=maxMode)
                return modes;
        }
    }

    return modes;
}
///project on the zerkine Basis (find a better name ? )
Eigen::VectorXd
Zernike::fitPhase(
    const Eigen::MatrixXd& phase,
    int maxMode)
{
    auto modes=generateModes(maxMode);

    vector<pair<int,int>> pixels;

    for(int y=0;y<H;y++)
    {
        for(int x=0;x<W;x++)
        {
            double xn=(x-CX)/R;
            double yn=(y-CY)/R;

            double rho=sqrt(xn*xn+yn*yn);

            if(rho<=1.0)
                pixels.push_back({x,y});
        }
    }

    int Npix=pixels.size();

    Eigen::MatrixXd A(Npix,maxMode);
    Eigen::VectorXd b(Npix);

    for(int p=0;p<Npix;p++)
    {
        int x=pixels[p].first;
        int y=pixels[p].second;

        double xn=(x-CX)/R;
        double yn=(y-CY)/R;

        double rho=sqrt(xn*xn+yn*yn);
        double theta=atan2(yn,xn);

        b(p)=phase(y,x);

        for(int k=0;k<maxMode;k++)
        {
            int n=modes[k].first;
            int m=modes[k].second;

            A(p,k)=evaluateZernike(n,m,
                           rho,
                           theta);
        }
    }

Eigen::VectorXd coeffs =
    A.colPivHouseholderQr().solve(b);//résoudre les moindres carrés b par décomposition QR: phase mesurée, A matrice contenant les polynome de zrnike, x=coef de zernike
//\phi_i=\sum_k(a_kZ_{k,i}) ou matriciellement Aa=b
    return coeffs;
}

Eigen::MatrixXd
Zernike::reconstruct(///reconstruire la phase aberrante->sommer les coef de zernike détectés
    const Eigen::VectorXd& coeffs)
{
    auto modes=
        generateModes(coeffs.size());

    Eigen::MatrixXd phase(H,W);

    phase.setZero();

    for(int y=0;y<H;y++)
    {
        for(int x=0;x<W;x++)
        {
            double xn=(x-CX)/R;
            double yn=(y-CY)/R;

            double rho=sqrt(xn*xn+yn*yn);

            if(rho<=1.0)
            {
                double theta=
                    atan2(yn,xn);

                double value=0;

                for(int k=0;
                    k<coeffs.size();
                    k++)
                {
                    int n=modes[k].first;
                    int m=modes[k].second;

                    value +=
                        coeffs(k)*
                        evaluateZernike(n,m,
                                 rho,
                                 theta);
                }

                phase(y,x)=value;
            }
        }
    }

    return phase;
}

Eigen::MatrixXcd
Zernike::correctWavefront(
    const Eigen::MatrixXcd& field,
    const Eigen::VectorXd& coeffs)
{
    Eigen::MatrixXd phase=
        reconstruct(coeffs);

    Eigen::MatrixXcd out(H,W);

    for(int y=0;y<H;y++)
    {
        for(int x=0;x<W;x++)
        {
            complex<double> corr=
                exp(complex<double>(
                    0,
                    -phase(y,x)));

            out(y,x)=field(y,x)*corr;
        }
    }

    return out;
}
