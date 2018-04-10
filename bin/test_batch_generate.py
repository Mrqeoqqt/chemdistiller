# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 15:51:17 2016

@author: Dr. Ivan Laponogov
"""

import sys;
if sys.byteorder!='little':
    print('Only little endian machines currently supported! bye bye ....');
    quit();
import os;
from random import randint;
from chemdistiller.msspectra.manager import SpectralManager;

print('Source path: %s'%os.getcwd());
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),'..')));

test_spectral_database_path='e:/Imperial/TestDB/NIST14_HMDB_MassBank';
test_spectral_database_path_output='e:/Imperial/TestDB/NIST14_HMDB_MassBank_5fold_CV';


nbatch=5;

specmanager=SpectralManager();
print('Reading test spectra...');
print(test_spectral_database_path);
specmanager.import_textfile_spectra_from_folder(test_spectral_database_path);
print('Finshed reading spectra. Total number: %s'%len(specmanager.ms_spectra));
print('Initializing batches');
batch={};
for i in range(len(specmanager.ms_spectra)):
    batch[specmanager.ms_spectra[i].parameters['global_index']]=randint(0,nbatch-1);
         
print('Prepared random batch indexes');
         
print(batch);
         
for i in range(len(specmanager.ms_spectra)):
    specmanager.ms_spectra[i].parameters['crossvalidation_batch_index']=batch[specmanager.ms_spectra[i].parameters['global_index']];

print('Exporting new test set');
         
specmanager.export_textfile_spectra_to_folder(test_spectral_database_path_output);
         
         
specmanager.close();
         
print('Finished');
     
