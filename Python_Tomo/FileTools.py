# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""

import os
from linecache import getline
import imageio as im
import numpy as np

def lineinfile(filename):
    """
    Determination of the number of lines in an ASCII file

    Parameters
    ----------
    filename : str
        File to be opened.

    Returns
    -------
    NumberOfLines : int
        Number of lines in the opened file.

    """
    if os.path.isfile(filename):
        File = open(filename, 'r')
        text=File.readlines()
        NumberOfLines = len(text)
        File.seek(0)
        File.close()
    else:
        print("File does not exist")
        File.seek(0)
        File.close()
    return NumberOfLines

def readvalue(filename,keyword):
    """
    Extraction of parameters values

    Parameters
    ----------
    filename : str
        File to be opened.
    keyword : str
        Parameter to be extracted.

    Returns
    -------
    Value : float
        Value of the requested parameter.

    """
    FileContent = []
    nbLignes = lineinfile(filename)
    nbLignesOK = 0
    Value = 0.
    for cpt in range(1,nbLignes+1):
        Ligne = getline(filename,cpt).split()
        if len(Ligne) != 0:
            FileContent.append(Ligne)
            nbLignesOK += 1
    for cpt in range(0,nbLignesOK):
        if FileContent[cpt][0] == keyword:
            Value = float(FileContent[cpt][1])
    return Value

def SAVbin(image,chemin,dim):
    """
    Binary file saving

    Parameters
    ----------
    image : at will !!
        Image to be saved as a binary file.
    chemin : str
        Path + name (without exension) of the file to be saved.
    dim : str
        x or y dimension of the image to be saved.

    Returns
    -------
    None.

    """
    nom_fichier = chemin + "_" + dim + ".bin"
    fid = open(nom_fichier,"w")
    image.tofile(fid)
    fid.close()
    
def SAVtiffCube(Folder,Data):
    """
    Saving data as a multipage Tiff file

    Parameters
    ----------
    Folder : str
        Path of the saved file.
    Data : float64
        Data to be saved.

    Returns
    -------
    None.

    """
    im.volwrite(Folder, Data.transpose((-1, 0, 1)).astype(np.float32))
    
def ReadtiffCube(Folder):
    """
    Reading data from a multipage Tiff file

    Parameters
    ----------
    Folder : str
        Path of the saved file.

    Returns
    -------
    Data : float32
        Extracted data cube.

    """
    Data = im.volread(Folder).transpose(1,-1,0)
    return Data