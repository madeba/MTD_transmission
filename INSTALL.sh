#!/bin/bash



rm ${HOME}/.config/gui_tomo.conf
touch ${HOME}/.config/gui_tomo.conf

DIR_TOMO=`pwd`
DIR_ACQUIS=${DIR_TOMO}/Tomo_acquis/
DIR_RESULT=${DIR_TOMO}/Tomo_result/
DIR_CONFIG=${DIR_TOMO}/Tomo_config/

echo "CHEMIN_RESULT "$DIR_RESULT>>$HOME/.config/gui_tomo.conf
echo "CHEMIN_ACQUIS "$DIR_ACQUIS>>$HOME/.config/gui_tomo.conf
echo "CHEMIN_CONFIG "$DIR_CONFIG>>$HOME/.config/gui_tomo.conf
echo "#uniquement sur pc acquisition">>$HOME/.config/gui_tomo.conf
echo "CHEMIN_CONFIG_PC_ACQUIS "$DIR_CONFIG>>$HOME/.config/gui_tomo.conf
echo "PC_ACQUIS 0">>$HOME/.config/gui_tomo.conf

sudo cp ${DIR_TOMO}/Tomo_prog/Gui_Tomo/bin/Release/tomo_gui /usr/bin/tomo_gui
sudo cp ${DIR_TOMO}/Tomo_prog/pretraitement_hors_axe/bin/Release/tomo_pretraitement /usr/bin/tomo_pretraitement
sudo cp ${DIR_TOMO}/Tomo_prog/Reconstruction/bin/Release/tomo_reconstruction /usr/bin/tomo_reconstruction
sudo cp ${DIR_TOMO}/Tomo_prog/Tomo_GPS/bin/Release/tomo_GPS /usr/bin/tomo_GPS
sudo cp ${DIR_TOMO}/utilTomo/show_fourier_PP/bin/Release/tomo_show_fourier /usr/bin/tomo_show_fourier
sudo cp ${DIR_TOMO}/utilTomo/set_Flower /usr/bin/tomo_set_flower




