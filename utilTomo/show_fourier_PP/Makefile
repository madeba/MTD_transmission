CFLAGS=-D_UNIX_ -D_LINUX_ -Iinclude  -I. -I/usr/include/opencv4 -I/opt/pleora/ebus_sdk/Ubuntu-x86_64/include -I/opt/pleora/ebus_sdk/Ubuntu-x86_64/lib


LDFLAGS=-L/opt/pleora/ebus_sdk/Ubuntu-x86_64/lib `pkg-config --libs opencv4` -ltiff   -lboost_system -lboost_filesystem -lboost_chrono  -lPvBase -lPvDevice -lPvBuffer -lPvGenICam -lPvStream -lPvTransmitter -lPvVirtualDevice -lPvCameraBridge -lPvAppUtils   -lPvPersistence -lPvSystem  -lPvSerial -lPvGUI -lSimpleImagingLib -lboost_thread -ltiff  -lfftw3 -lfftw3_threads 

OBJDIR   = obj/Release
CC=g++

tomo_show_fourier : main.o  fonctions.o manip.o 
	g++ -o bin/Release/tomo_show_fourier    $(OBJDIR)/main.o  $(OBJDIR)/src/fonctions.o $(OBJDIR)/src/manip.o  -O3  $(LDFLAGS)
	
main.o : main.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@ 
  

	
fonctions.o : src/fonctions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/src/$@  

	
manip.o : src/manip.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/src/$@	


	
clean :
	rm -f obj/Release/*.o obj/Release/src/*.o bin/Release/tomo_show_fourier
mrproper: clean
	rm -f bin/Release/tomo_show_fourier
