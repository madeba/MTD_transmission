#include "Point3D.h"

using namespace std;

Point3D::Point3D()
{
    //ctor

}

Point3D::Point3D(double x,double y,double z)
{
    this->x=x;
    this->y=y;
    this->z=z;
}

Point3D::Point3D(double x,double y,double z,int neo_dim3D)
{
    this->x=x;
    this->y=y;
    this->z=z;
    this->dim3D=neo_dim3D;
}

Point3D::Point3D(Point2D Pt2D,double z,int dim3D)
{
    this->x=Pt2D.x;
    this->y=Pt2D.y;
    this->z=z;
    this->dim3D=dim3D;
}

Point3D::~Point3D()
{
    //dtor
}

int Point3D::get_dim3D()
{
return this->dim3D;
}

size_t Point3D::cpt3D()
{
    size_t cpt=0;
    return cpt=round(z)*dim3D*dim3D+round(y)*dim3D+round(x);
}

//cpt repère humain centré vers repère informatique
Point3D Point3D::coordI()
{
    Point3D coordI(this->x+dim3D/2,-this->y+dim3D/2,this->z+dim3D/2,dim3D);
    return coordI;
}

Point3D Point3D::coordH()
{
    Point3D coordH((int)(this->x)%dim3D,(int)(this->y)/dim3D,(int)(this->z)/(dim3D*dim3D),dim3D);
    return coordH;
}

void Point3D::decale(int decalx,int decaly,int decalz)
{
 this->x=x+decalx;
 this->y=y+decaly;
 this->z=z+decalz;
}

Point3D Point3D::operator-(Point3D const &Pt2)
{
Point3D result(0,0,0,0);
result.x=this->x-Pt2.x;
result.y=this->y-Pt2.y;
result.z=this->z-Pt2.z;
result.dim3D=this->dim3D;/// /!\ :normalement, on connait pas l'espace du résultat

return result;
}

Point3D Point3D::sym_xoz(){
Point3D sym;
sym.x=this->x;
sym.y=this->y;
sym.z=this->z;
return sym;
}
