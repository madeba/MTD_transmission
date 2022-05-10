CFLAGS=-Iinclude -I. -I/usr/include/opencv4
LDFLAGS=-lpthread -lm `pkg-config --libs opencv4` -ltiff  -lfftw3 -lfftw3_threads -laf
OBJDIR   = obj/Release
CC=g++

pretraitement_hors_axe : main.o Correction_aberration.o deroulement_herraez.o fonctions.o FFT_fonctions.o  GPU_fonctions.o manip.o  deroulement_volkov3_GPU.o deroulement_volkov4.o FFTW_init.o IO_fonctions.o vecteur.o Point2D.o

	g++ -o bin/Release/tomo_pretraitement $(OBJDIR)/Correction_aberration.o $(OBJDIR)/src/deroulement_volkov3_GPU.o $(OBJDIR)/deroulement_herraez.o obj/Release/FFT_fonctions.o obj/Release/fonctions.o obj/Release/GPU_fonctions.o $(OBJDIR)/main.o $(OBJDIR)/manip.o  $(OBJDIR)/src/deroulement_volkov4.o $(OBJDIR)/src/FFTW_init.o $(OBJDIR)/src/IO_fonctions.o $(OBJDIR)/src/Point2D.o $(OBJDIR)/src/vecteur.o -O3  $(LDFLAGS)
	

main.o : main.cpp
	g++ -c $< $(CFLAGS) -o $(OBJDIR)/$@ 
		 
Correction_aberration.o : Correction_aberration.cpp
	g++ -c $< $(CFLAGS) -o $(OBJDIR)/$@ 



deroulement_herraez.o : deroulement_herraez.cpp
	g++ -c $< $(CFLAGS)   -o $(OBJDIR)/$@ 

FFT_fonctions.o : FFT_fonctions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@ 
	
fonctions.o : fonctions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@ 
	
GPU_fonctions.o : src/GPU_fonctions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@ 
	
manip.o : manip.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@  

deroulement_volkov4.o : src/deroulement_volkov4.cpp	
	g++  -c $< $(CFLAGS) -o $(OBJDIR)/src/$@ 

deroulement_volkov3_GPU.o : src/deroulement_volkov3_GPU.cpp	
	g++  -c $< $(CFLAGS) -o $(OBJDIR)/src/$@ 

FFTW_init.o : src/FFTW_init.cpp
	g++  -c $<  $(CFLAGS) -o $(OBJDIR)/src/$@ 
	
IO_fonctions.o : src/IO_fonctions.cpp
	g++  -c $<  $(CFLAGS) -o $(OBJDIR)/src/$@ 
	
vecteur.o : src/vecteur.cpp
	g++  -c $<  $(CFLAGS) -o $(OBJDIR)/src/$@ 

Point2D.o : src/Point2D.cpp
	g++  -c $<  $(CFLAGS) -o $(OBJDIR)/src/$@ 
	
	
clean :
	rm -f obj/Release/*.o obj/Release/src/*.o bin/Release/tomo_pretraitement_GPU
mrproper: clean
	rm -f bin/Release/pretraitement_hors_axe_GPU
