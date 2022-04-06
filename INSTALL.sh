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

sudo cp ./bin/tomo_gui /usr/bin/tomo_gui
sudo cp ./bin/tomo_pretraitement /usr/bin/tomo_pretraitement
sudo cp ./bin/tomo_reconstruction /usr/bin/tomo_reconstruction



