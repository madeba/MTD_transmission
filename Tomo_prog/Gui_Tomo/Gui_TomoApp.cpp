/***************************************************************
 * Name:      Gui_TomoApp.cpp
 * Purpose:   Code for Application Class
 * Author:    Matthieu Debailleul ()
 * Created:   2017-03-28
 * Copyright: Matthieu Debailleul ()
 * License:
 **************************************************************/

#ifdef WX_PRECOMP
#include "wx_pch.h"
#endif

#ifdef __BORLANDC__
#pragma hdrstop
#endif //__BORLANDC__

#include "Gui_TomoApp.h"
#include "Gui_TomoMain.h"
#include <iostream>

using namespace std;

IMPLEMENT_APP(Gui_TomoApp);

bool Gui_TomoApp::OnInit()
{
    Gui_TomoFrame* frame = new Gui_TomoFrame(0L, _("Interface Tomo")) ;

    wxSize frame_size(900,460);
    frame->SetSize(frame_size);
    frame->SetMaxSize(wxSize(960,520));
    frame->Show(true);
    return true;
}
