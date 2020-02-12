// this code belongs to:
// http://www.chasanc.com/index.php/Coding/PGM-Image-Format.html


#ifndef __VECTRA_PGM__
#define __VECTRA_PGM__


#define VECTRA_UCHAR unsigned char


int PGM_Read_Header(const char* FileName,
                    long *width,
                    long *height,
                    long *max_color,
                    long *f_start);

void
PGM_write_Header(FILE* file,
		 long width,
                 long height,
                 long max_color);


/* lit les données un fichier PGM binaire 8 bits dans une matrice allouée */
int PGM_8b_Get_Data(const char* FileName,
                 long width,
                 long height,
                 long f_start,
                 VECTRA_UCHAR **img_data);


/* lit les données un fichier PGM binaire 8 bits dans un tableau de données alloué en row-major */
int PGM_8b_Get_Data_Linear(const char* FileName,
			long width,
			long height,
			long f_start,
			VECTRA_UCHAR *img_data);


/* écrit un fichier PGM binaire 8 bits de nom spécifié à partir d'une matrice de données */
int PGM_8b_Put_Data(const char* FileName,
                 long width,
                 long height,
                 long max_color,
                 VECTRA_UCHAR **img_data);


/* écrit un fichier PGM binaire 8 bits de nom spécifié à partir du tableau de données en row-major */
int PGM_8b_Put_Data_Linear(const char *FileName,
			long width,
			long height,
			long max_color,
			VECTRA_UCHAR *img_data);


#endif /* __VECTRA_PGM__ */ 
