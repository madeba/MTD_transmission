#!/bin/bash


rm ${HOME}/.config/gui_tomo.conf

DIR_TOMO=`pwd`
DIR_ACQUIS=${DIR_TOMO}/Tomo_acquis/
DIR_RESULT=${DIR_TOMO}/Tomo_result/
DIR_CONFIG=${DIR_TOMO}/Tomo_config/

sudo rm ./bin/tomo_gui /usr/bin/tomo_gui
sudo rm ./bin/tomo_pretraitement /usr/bin/tomo_pretraitement
sudo rm ./bin/tomo_reconstruction /usr/bin/tomo_reconstruction



