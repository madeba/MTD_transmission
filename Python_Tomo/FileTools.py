# -*- coding: utf-8 -*-
import os
from linecache import getline

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

def filetolist(filename):
    FileContent = []
    Ligne = []
    nbLignes = lineinfile(filename)
    for cpt in range(1,nbLignes):
        Ligne = getline(filename,cpt).split()
        # Ligne[0]
        FileContent.append(Ligne)  
    return FileContent