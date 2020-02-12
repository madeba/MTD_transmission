/***************************************************************
 * Name:      Gui_TomoApp.h
 * Purpose:   Defines Application Class
 * Author:    Matthieu Debailleul ()
 * Created:   2017-03-28
 * Copyright: Matthieu Debailleul ()
 * License:
 **************************************************************/

#ifndef GUI_TOMOAPP_H
#define GUI_TOMOAPP_H

#include <wx/app.h>

class Gui_TomoApp : public wxApp
{
    public:
        virtual bool OnInit();
};

#endif // GUI_TOMOAPP_H
