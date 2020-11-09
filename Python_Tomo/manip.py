#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 13:01:32 2020

@author: nicolas
"""
import os
import FileTools as ft

class Manip:
    """
    Class for data processing parametrization
    """
    def __init__(self, dossier_acquis, isdata):
        """
        Definition of the data acquisition parameters

        Parameters
        ----------
        dossier_acquis : str
            Path to the data.
        isdata : bool
            If isdata is True, data are processed, if isdata is False, white images are processed.

        Returns
        -------
        None.

        """
        self.dossier_acquis = dossier_acquis
        self.dossier_data = f"{dossier_acquis}data/"
        self.dossier_blanc = f"{dossier_acquis}blanc/"
        self.isdata = isdata
        self.fichier_config = f"{dossier_acquis}config/config_manip.txt"
        procfolder = self.createresultfolder(dossier_acquis, "Pretraitement", isdata)
        self.dossier_pretraitement = procfolder
        recfolder = self.createresultfolder(dossier_acquis, "Reconstruction", isdata)
        self.dossier_reconstruction = recfolder
        gerchfolder = self.createresultfolder(dossier_acquis, "Gerchberg", isdata)
        self.dossier_gerchberg = gerchfolder
        self.NA = float(ft.readvalue(self.fichier_config, 'NA'))
        self.NIMM = float(ft.readvalue(self.fichier_config, 'N0'))
        self.LAMBDA = float(ft.readvalue(self.fichier_config, 'LAMBDA'))
        self.F_TUBE = float(ft.readvalue(self.fichier_config, 'F_TUBE')) # Tube lens focal length
        self.F_OBJ = float(ft.readvalue(self.fichier_config, 'F_OBJ')) # Microscope objective focal length
        self.PIX = float(ft.readvalue(self.fichier_config, 'TPCAM')) # Physical pixel pitch
        self.RAPFOC = float(ft.readvalue(self.fichier_config, 'RF')) # focal length ratio of the resampling lens dublet
        self.CHEMINMASQUE = f"{self.dossier_data}Image_maask.pgm"
        self.CENTREX = int(ft.readvalue(self.fichier_config, 'CIRCLE_CX')) # Pupil center in Fourier space
        self.CENTREY = int(ft.readvalue(self.fichier_config, 'CIRCLE_CY'))
        self.NB_HOLOTOT = int(ft.readvalue(self.fichier_config, 'NB_HOLO'))

    def createresultfolder(self, dossier_acquis, foldername, isdata):
        """
        Creation of the result folders

        Parameters
        ----------
        dossier_acquis : str
            Path to the data.
        foldername : str
            Name of the folder to be created.
        isdata : bool
            If isdata is True, folder is created in the data tree, if isdata is False, folder is created in the white tree.

        Returns
        -------
        dossier : str
            Path to the created folder.

        """
        dossier_data = f"{dossier_acquis}data/"
        dossier_blanc = f"{dossier_acquis}blanc/"
        if isdata is True:
            if not os.path.exists(f"{dossier_data}{foldername}/"):
                os.makedirs(f"{dossier_data}{foldername}/")
            dossier = f"{dossier_data}{foldername}/"
        else:
            if not os.path.exists(f"{dossier_blanc}{foldername}/"):
                os.makedirs(f"{dossier_blanc}{foldername}/")
            dossier = f"{dossier_data}{foldername}/"
        return dossier