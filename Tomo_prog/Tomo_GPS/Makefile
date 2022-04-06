CFLAGS=-std=c++11 -mveclibabi=svml -march=corei7-avx -fexpensive-optimizations -O3 -I/usr/include/opencv -I./include -I/usr/include/opencv4


LDFLAGS=-lpthread -lm -lgomp  -lopencv_core -lopencv_imgproc -lopencv_imgcodecs -lopencv_highgui -ltiff  -lfftw3 -lfftw3_threads 

OBJDIR   = obj/Release
CC=g++

tomo_GPS : main.o  FFT_fonctions.o FFTW_init.o fonctions.o IO_fonctions.o manip.o Point2D.o Point3D.o
	g++ -o bin/Release/tomo_GPS    $(OBJDIR)/main.o  $(OBJDIR)/src/FFT_fonctions.o $(OBJDIR)/src/FFTW_init.o $(OBJDIR)/src/fonctions.o $(OBJDIR)/src/IO_fonctions.o $(OBJDIR)/src/manip.o $(OBJDIR)/src/Point2D.o $(OBJDIR)/src/Point3D.o -O3  $(LDFLAGS)
	
main.o : main.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@ 
  

FFT_fonctions.o : src/FFT_fonctions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/src/$@  
	
FFTW_init.o : src/FFTW_init.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/src/$@  
	
fonctions.o : src/fonctions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/src/$@  

IO_fonctions.o : src/IO_fonctions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/src/$@
	
manip.o : src/manip.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/src/$@	

Point2D.o : src/Point2D.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/src/$@
		
Point3D.o : src/Point3D.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/src/$@
	
clean :
	rm -f obj/Release/*.o obj/Release/src/*.o bin/Release/tomo_GPS
mrproper: clean
	rm -f bin/Release/tomo_GPS
