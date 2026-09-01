#include "Point2D.h"


Point2D::Point2D()
{
    //ctor

}

/*Point2D::Point2D(double x, double y,int m_dim2D)
{
    this->x=x;
    this->y=y;
    this->dim2D=m_dim2D;
}*/


void Point2D::set_xy(double x, double y)
{
    this->x=x;
    this->y=y;

}



Point2D::~Point2D()
{
    //dtor
}

size_t Point2D::cpt2D()
{
    size_t cpt=0;
    return cpt=round(y)*dim2D+round(x);
}
//cpt repère humain centré vers repère informatique
Point2D Point2D::coordI()
{
    Point2D coordI(this->x+dim2D/2,-this->y+dim2D/2, dim2D);
    return coordI;
}

/*Point2D Point2D::coordI2coordH()
{
    Point2D coordI(this->x+dim2D/2,-this->y+dim2D/2,dim2D);
    return coordI;
}*/

void Point2D::decale(int decalx, int decaly)
{
 this->x=x+decalx;
 this->y=y+decaly;
}



