# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 18:15:32 2016

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



from chemdistiller.msspectra.spectrum import MSSpectrum;
from chemdistiller.utils.sysutils import print_in_the_same_line;

class SpectralManager:
    def __init__(self):
        self.ms_spectra=[];

    def close(self):
        for spectrum in self.ms_spectra:
            spectrum.close();
        
        
    def __enter__(self):
        return self;

    def __exit__(self, exc_type, exc_value, traceback):
        self.close();
        
    def add_spectrum_from_textfile(self, fname):
            spectrum=MSSpectrum();
            n,e=os.path.splitext(os.path.basename(fname));
            spectrum.parameters['filename']=n;
            spectrum.load_from_textfile(fname);
            self.ms_spectra.append(spectrum);
            
    def import_textfile_spectra_from_folder(self, dirname):
        print('Importing from %s'%dirname);
        print('\n');
        onlyfiles = [f for f in os.listdir(dirname) if os.path.isfile(os.path.join(dirname, f))];
        for fname in onlyfiles:
            _,fileext=os.path.splitext(fname);
            if fileext.lower()=='.txt':
                self.add_spectrum_from_textfile(os.path.join(dirname,fname));
                print_in_the_same_line('\r%s spectra loaded'%len(self.ms_spectra));
        print('\n');        
                
    def export_textfile_spectra_to_folder(self, dirname):
        print('Exporting to %s'%dirname)
        print('\n');
        if not os.path.exists(dirname):
            os.makedirs(dirname);
        cc=0;
        for spectrum in self.ms_spectra:
            cc+=1;
            if 'filename' in spectrum.parameters:
                spectrum.save_to_textfile(os.path.join(dirname,'%s.txt'%spectrum.parameters['filename']));
            else:
                spectrum.save_to_textfile(os.path.join(dirname,'%s.txt'%cc));
            print_in_the_same_line('\r%s spectra saved'%cc);
        print('\n');

            
if __name__=='__main__':

    spectral_manager=SpectralManager();
    print(os.path.join(chemdistiller_path,'data/MassBankTestSpectra'));
    spectral_manager.import_textfile_spectra_from_folder(os.path.join(chemdistiller_path,'data/MassBankTestSpectra'));
    print(len(spectral_manager.ms_spectra));
    print(spectral_manager.ms_spectra[0].parameters);
    spectral_manager.close();
    
    