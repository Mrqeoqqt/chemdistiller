# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 19:36:52 2016

@author: ilaponog
"""

class MSPeak:
    def __init__(self, basic_peak=None):
        if basic_peak is None:
            self.mz=0.0;
            self.intensity=0.0;
        else:
            self.mz=basic_peak.mz;
            self.intensity=basic_peak.intensity;
        self.parent_spectrum=None;
        

    def close(self):
        if hasattr(self,'ms_spectra'):
            for spectrum in self.ms_spectra:
                spectrum.close();
            
    def __enter__(self):
        return self;

    def __exit__(self, exc_type, exc_value, traceback):
        self.close();
        
    def add_spectrum(self, spectrum):
        self.ms_spectra.append(spectrum);
        spectrum.parent_peak=self;

class MSBasicPeak:
    __slots__ = ['mz', 'intensity'];
    