# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:25:44 2017

@author: Dr. Ivan Laponogov
"""


import sys;
import os;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);


from chemdistiller.msspectra.peak import MSPeak, MSBasicPeak;


def purge_ms1(spectral_manager):
    for j in reversed(range(len(spectral_manager.ms_spectra))):
        spectrum = spectral_manager.ms_spectra[j]
        for i in reversed(range(len(spectrum.peaks))):
            peak = spectrum.peaks[i];
            if isinstance(peak, MSBasicPeak):
                del spectrum.peaks[i];
            elif hasattr(peak,'ms_spectra'):
                if len(peak.ms_spectra) == 0:
                    del spectrum.peaks[i];
            else:
                del spectrum.peaks[i]
            
        if len(spectrum.peaks) == 0:
            del spectral_manager.ms_spectra[j];
            


def merge_ms1(spectral_manager, parameters):
    pass

