#ifndef DEF_STRUCT// Si la constante n'a pas été définie` le fichier n'a jamais été inclus
#define DEF_STRUCT


//Nombre complexe


typedef struct {
  double Re,Im;
}nbCplx;


typedef struct {
  int x,y;
 // nbCplx Val;
}Var2D;

typedef struct {
  int x,y,z;
}Var3D;
#endif // DEF_STRUCT
