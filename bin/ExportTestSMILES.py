# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 15:51:17 2016

@author: Dr. Ivan Laponogov
"""
import sys;
if __name__=='__main__':

    
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

import os;
chemdistiller_path = os.path.abspath(os.path.join(os.getcwd(), ".."));
sys.path.append(chemdistiller_path);


from chemdistiller.msspectra.manager import SpectralManager;

if __name__=='__main__':
    
    test_spectral_database_path='e:/Imperial/TestDB/NIST14_HMDB_MassBank_5fold_CV_FPT3';
    
    specmanager=SpectralManager();
    print('Reading test spectra...');
    
    specmanager.import_textfile_spectra_from_folder(test_spectral_database_path);
    smilesset=set([]);
    
    
    for spectrum in specmanager.ms_spectra:
        smiles=spectrum.parameters['smiles'];
        smilesset.add(smiles);
        
    fout=open('e:/Imperial/TestDB/TestDBSmiles.txt','w');
    cc=0;
    for smiles in smilesset:
        cc+=1;
        fout.write('%s\t%s\n'%(smiles,cc));
        
        
    fout.close();        
    
    specmanager.close();
