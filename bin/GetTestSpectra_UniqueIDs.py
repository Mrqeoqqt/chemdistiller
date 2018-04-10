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



def remove_low_precision_spectra(precision):
        print('Removing low precision spectra');
        
        
        for i in reversed(range(len(specmanager.ms_spectra))):
            spectrum=specmanager.ms_spectra[i];
            
            mass=float(spectrum.parameters['exactmass']);
            mode=spectrum.parameters['mode'];
            keepspectrum=False;
            
            for peakindex in reversed(range(len(spectrum.peaks))):
                keep=False;
                peak=spectrum.peaks[peakindex];
                if hasattr(peak,'parameters'):
                    if 'ion_type' in peak.parameters:            
                        adduct=get_adduct_by_name(peak.parameters['ion_type'],mode);
                        if not (adduct is None):
                            keep=True;
                            mz=adduct.get_mz(mass);
                            minmzd=1.0;
                            for subspectrum in peak.ms_spectra:
                                if abs(mz-float(subspectrum.parameters['precursor_mz']))>precision:
                                    keep=False;
                                    break;
                                else:
                                    for subpeak in subspectrum.peaks:
                                        if abs(subpeak.mz-mz)<minmzd:
                                            minmzd=abs(subpeak.mz-mz);
                            if keep and (minmzd>precision):
                                keep=False;
                if keep==False:
                    del spectrum.peaks[peakindex];
                else:
                    keepspectrum=True;
                    
            if keepspectrum==False:
                del specmanager.ms_spectra[i];
                
        
                        
        print('Finshed removing spectra. Total number left: %s'%len(specmanager.ms_spectra));    
    
    
def correct_precursor_mz(spectral_manager):
        
        for spectrum in specmanager.ms_spectra:
            mass=float(spectrum.parameters['exactmass']);
            mode=spectrum.parameters['mode'];
            for peak in spectrum.peaks:
                if 'ion_type' in peak.parameters:            
                    adduct=get_adduct_by_name(peak.parameters['ion_type'],mode);
                    if not (adduct is None):
                        peak.mz=adduct.get_mz(mass);
    
    

if __name__=='__main__':

    
    test_spectral_database_path='e:/Imperial/TestDB/NIST14_HMDB_MassBank_5fold_CV_FPT3';
    
    ppm=21;

    print('Starting. %s'%time.strftime('%d/%m/%y %H:%M:%S'));
    starttime=time.time();

    precision=0.005;
    
    specmanager=SpectralManager();
        
    specmanager.import_textfile_spectra_from_folder(test_spectral_database_path);
        
    remove_low_precision_spectra(precision);
    correct_precursor_mz(specmanager);    
    
    peaks=[];
    
    for spectrum in specmanager.ms_spectra:
        if spectrum.parameters['charge']=='0':
            if int(spectrum.parameters['crossvalidation_batch_index'])==0:
                for peak in spectrum.peaks:
                    if hasattr(peak,'parameters'):
                        if 'ion_type' in peak.parameters:
                           if (peak.parameters['ion_type'] in FragPrintScorer.supported_adducts):
                               peaks.append(peak);
                
        
    peaks=sorted(peaks,key=lambda peak: peak.mz);
    
    neg=[];
    pos=[];
    for peak in peaks:
        if peak.parent_spectrum.parameters['mode']==-1:
            neg.append(peak);
        else:
            pos.append(peak);
    
    print('Positive mode: %s'%len(pos));
    print('Negative mode: %s'%len(neg));    
    
    
                 
    specmanager.close();
    
        