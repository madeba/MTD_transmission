CFLAGS=-I/usr/lib/x86_64-linux-gnu/wx/include/gtk3-unicode-3.0 -I/usr/include/wx-3.0 -D_FILE_OFFSET_BITS=64 -DWXUSINGDLL -D__WXGTK__ -O2 -I/usr/include/wx-3.0 -I/usr/lib/x86_64-linux-gnu/wx/include/base-unicode-3.0 


LDFLAGS=-L/usr/lib/x86_64-linux-gnu -pthread   -lwx_gtk3u_xrc-3.0 -lwx_gtk3u_webview-3.0 -lwx_gtk3u_html-3.0 -lwx_gtk3u_qa-3.0 -lwx_gtk3u_adv-3.0 -lwx_gtk3u_core-3.0 -lwx_baseu_xml-3.0 -lwx_baseu_net-3.0 -lwx_baseu-3.0  -L/usr/lib/x86_64-linux-gnu/

OBJDIR   = obj/Release
CC=g++

tomo_gui : Gui_TomoApp.o Gui_TomoMain.o
	g++ -o bin/Release/tomo_gui    $(OBJDIR)/Gui_TomoApp.o   $(OBJDIR)/Gui_TomoMain.o -O3  $(LDFLAGS)
	
Gui_TomoApp.o : Gui_TomoApp.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@ 
  

Gui_TomoMain.o : Gui_TomoMain.cpp
	g++ -c $< $(CFLAGS)  -o $(OBJDIR)/$@  
 

	
clean :
	rm -f obj/Release/*.o obj/Release/src/*.o bin/Release/tomo_manip
mrproper: clean
	rm -f bin/Release/tomo_gui
