CFLAGS=-Iinclude -I. -I/usr/include/opencv4
LDFLAGS=-lpthread -lm `pkg-config --libs opencv4` -ltiff  -lfftw3 -lfftw3_threads
OBJDIR   = obj/Release
CC=g++

reconstruction : main.o  fonctions.o FFT_fonctions.o  manip.o   FFT_encaps.o IO_fonctions.o  Point2D.o Point3D.o

	g++ -o bin/Release/tomo_reconstruction   obj/Release/FFT_fonctions.o obj/Release/fonctions.o $(OBJDIR)/main.o $(OBJDIR)/manip.o   $(OBJDIR)/src/FFT_encaps.o $(OBJDIR)/src/IO_fonctions.o $(OBJDIR)/src/Point2D.o $(OBJDIR)/src/Point3D.o -O3  $(LDFLAGS)
	

main.o : main.cpp
	g++ -c $< $(CFLAGS) -o $(OBJDIR)/$@ 
		 

FFT_fonctions.o : src/FFT_fonctions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@ 

fonctions.o : src/fonctions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@ 

  

manip.o : src/manip.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@  
 

FFT_encaps.o : src/FFT_encaps.cpp
	g++  -c $<  $(CFLAGS) -o $(OBJDIR)/src/$@ 
	
IO_fonctions.o : src/IO_fonctions.cpp
	g++  -c $<  $(CFLAGS) -o $(OBJDIR)/src/$@ 
	

Point2D.o : src/Point2D.cpp
	g++  -c $<  $(CFLAGS) -o $(OBJDIR)/src/$@ 

Point3D.o : src/Point3D.cpp
	g++  -c $<  $(CFLAGS) -o $(OBJDIR)/src/$@ 	
	
clean :
	rm -f obj/Release/*.o obj/Release/src/*.o bin/Release/tomo_reconstruction
mrproper: clean
	rm -f bin/Release/tomo_reconstruction
