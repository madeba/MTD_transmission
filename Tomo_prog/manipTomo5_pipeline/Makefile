CFLAGS=-D_UNIX_ -D_LINUX_ -Iinclude -Iinclude/LabJack -Iinclude/john -I. -I/usr/include/opencv4 -I/opt/pleora/ebus_sdk/Ubuntu-x86_64/include -I/opt/pleora/ebus_sdk/Ubuntu-x86_64/lib


LDFLAGS=-L/opt/pleora/ebus_sdk/Ubuntu-x86_64/lib -L/opt/lib_tomo/labjack -lpthread -lm `pkg-config --libs opencv4` -ltiff   -lboost_system -lboost_filesystem -lboost_chrono  -lPvBase -lPvDevice -lPvBuffer -lPvGenICam -lPvStream -lPvTransmitter -lPvVirtualDevice -lPvCameraBridge -lPvAppUtils   -lPvPersistence -lPvSystem  -lPvSerial -lPvGUI -lSimpleImagingLib -lboost_thread      -lLjack -llabjackusb -llabjackusbcpp -lmsleep

OBJDIR   = obj/Release
CC=g++

tomo_manip : main.o  fonctions.o   scan_functions.o 
	g++ -o bin/Release/tomo_manip   $(OBJDIR)/main.o $(OBJDIR)/fonctions.o   $(OBJDIR)/scan_functions.o -O3  $(LDFLAGS)
	
fonctions.o : src/fonctions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@ 
  

scan_functions.o : src/scan_functions.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@  
 
main.o : main.cpp
	g++ -c $< $(CFLAGS) -o $(OBJDIR)/$@ 
	
clean :
	rm -f obj/Release/*.o obj/Release/src/*.o bin/Release/tomo_manip
mrproper: clean
	rm -f bin/Release/tomo_manip
