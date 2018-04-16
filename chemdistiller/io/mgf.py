# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:26:25 2017

@author: Dr. Ivan Laponogov
"""

import os
from chemdistiller.msspectra.spectrum import MSSpectrum
from MGF2ChemDistiller import MGFParser

def import_mgf(fname, spectral_manager, parameters=None):
    # author : Mingyi Xue
    # haven't been tested
    # when this program occurs error ,
    # remember to delete all code written by Mingyi
    # and roll back to the original one
    parser = MGFParser()
    fname = parser.mgf2cd(fname)#convert mgf to cd chemdistiller input
    for file in fname:
        spectrum = MSSpectrum();
        n, e = os.path.splitext(os.path.basename(fname));
        spectrum.parameters['filename'] = n;
        spectrum.load_from_textfile(fname);
        spectral_manager.ms_spectra.append(spectrum);




