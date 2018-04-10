# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 11:34:09 2017

@author: ilaponog
"""

import sys;
import os;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."));
    sys.path.append(chemdistiller_path);


import time;

from chemdistiller.settings import test_as_single_process;
from chemdistiller.msspectra.manager import SpectralManager;
from chemdistiller.msspectra.adducts import get_adduct_by_name;
from chemdistiller.annotations.peakannotator import PeakAnnotator;
from chemdistiller.scorers.totalscore import total_multiplicative_score;
from chemdistiller.scorers.fragprint import FragPrintScorer;
from chemdistiller.filters.formula import FormulasFilter;

from chemdistiller.io.jsonio import spectra_to_json;

from chemdistiller.annotations.peakannotation import MSPeakAnnotation, merge_annotations;

from operator import itemgetter;

import multiprocessing;



if __name__=='__main__':

    
    test_spectral_database_path='e:/Imperial/BitBucket/temp/annotated_spectra';
    test_spectral_database_path_out='e:/Imperial/BitBucket/temp/annotated_spectra_out';
    
    
    specmanager=SpectralManager();
        
    specmanager.import_textfile_spectra_from_folder(test_spectral_database_path);
    specmanager.export_textfile_spectra_to_folder(test_spectral_database_path_out);
        
             
    specmanager.close();
    
        