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
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../chemdistiller")));


import numpy as np;
from chemdistiller.scorers.formula import FormulaScorer;
from chemdistiller.scorers.element import ElementScorer;
from chemdistiller.chemdb.manager import DBManager;
from chemdistiller.filters.formula import FormulasFilter;
from chemdistiller.filters.inchi import InChiFilter;
from chemdistiller.msspectra.manager import SpectralManager;
#%%
from chemdistiller.msspectra.spectrum import MSSpectrum;
#%%
from chemdistiller.utils.base64 import encode_to_base64;

#%%
from chemdistiller.utils.base64 import decode_from_base64;
#%%


print('Source path: %s'%os.getcwd());
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),'..')));

test_spectral_database_path='e:/Imperial/TestDB/NIST14_HMDB_MassBank_5fold_CV';
#%%
test_spectral_database_outpath='e:/Imperial/TestDB/NIST14_HMDB_MassBank_5fold_CV_FPT';
#%%
test_fpt_stat_outfile='e:/Imperial/TestDB/FPT.txt';
#%%
test_chemical_databases=['HMDB','PubChem','ChEBI','MassBank'];

dbfile=os.path.abspath('../tests/databases.list');

nbatch=5;
    
print(dbfile);
if not os.path.isfile(dbfile):
    print('No database file found!');
    quit();
    
    
dbmanager=DBManager(dbfile);
print('Initialized Compound Databases:');
dbnames=dbmanager.get_database_names();
for s in dbnames:
      if s in test_chemical_databases:
           print('%s [*]'%s);
      else:
           print('%s'%s);                 
                 
specmanager=SpectralManager();
print('Reading test spectra...');
print(test_spectral_database_path);
specmanager.import_textfile_spectra_from_folder(test_spectral_database_path);
print('Finshed reading spectra. Total number: %s'%len(specmanager.ms_spectra));

print('Initializing batches');

batch=[];
for i in range(nbatch):
    batch.append([]);

for i in range(len(specmanager.ms_spectra)):
      batch[int(specmanager.ms_spectra[i].parameters['crossvalidation_batch_index'])].append(i);
         
print(batch);


#%%

missing=0;

fpts=np.zeros((11416,),dtype=np.uint32);



for spectrum_index in range(len(specmanager.ms_spectra)):
    spectrum=specmanager.ms_spectra[spectrum_index];
    print('Spectrum %s of %s...'%(spectrum_index,len(specmanager.ms_spectra)));
    shortinchi=spectrum.parameters['shortinchi'];
    mass=float(spectrum.parameters['exactmass']);
    charge=int(spectrum.parameters['charge']);
    if charge!=0:
        mz=abs(mass/charge);
    else:
        mz=mass;
    inchifilter=InChiFilter(shortinchi, use_short_inchi=True, match_type=4);
    
    result=dbmanager.query_by_mz_scored(mz, 20, charge, db_indexes=dbmanager.db_indexes_from_db_names(test_chemical_databases,case_sensitive=False),\
    filters=[inchifilter], scorers=[], required_fields=set(['FPT']), results_limit=-1, save_memory=False);
    ff=[];
    if len(result.mol_list)==0:
        missing+=1;
        print('Missing: %s'%missing);
    else:
        for mol in result.mol_list:
            ff.append(encode_to_base64(mol['FPT']));
        ff=set(ff);
        spectrum.parameters['FPTCount']=len(ff);
        cc=-1;
        for subfpt in ff:
            cc+=1;
            spectrum.parameters['FPT_%s'%cc]=subfpt;
            subfpt=decode_from_base64(subfpt);
            subfpt=np.unpackbits(subfpt);
            fpts=np.add(fpts,subfpt);
#%%            
print('Exporting spectra');
specmanager.export_textfile_spectra_to_folder(test_spectral_database_outpath);         
#%%
#test_fpt_stat_outfile='e:/Imperial/TestDB/FPT.txt';

print('Exporting FPT stats');
fout=open(test_fpt_stat_outfile,'w');
for i in range(11416):
    fout.write('%s\t%s\n'%(i,fpts[i]));
fout.close();




#%%
for spectrum_index in range(len(specmanager.ms_spectra)):
    spectrum=specmanager.ms_spectra[spectrum_index];
    if not ('FPTCount' in spectrum.parameters.keys()):
        formulafilter=FormulasFilter(spectrum.parameters['formula']);

        shortinchi=spectrum.parameters['shortinchi'];
        mass=float(spectrum.parameters['exactmass']);
        charge=int(spectrum.parameters['charge']);
        if charge!=0:
            mz=abs(mass/charge);
        else:
            mz=mass;
            
        print('shortinchi: %s'%shortinchi);
        print('mz: %s'%mz);
        print('charge: %s'%charge);
        
        result=dbmanager.query_by_mz_scored(mz, 2, charge, db_indexes=dbmanager.db_indexes_from_db_names(test_chemical_databases,case_sensitive=False),\
        filters=[formulafilter], scorers=[], required_fields=set(['FPT']), results_limit=-1, save_memory=False);
        #print(result.mol_list);
        for s in result.mol_list:
            print(s['ShortInChi']);
        break;

#%%
totalfpt=0;
for spectrum_index in range(len(specmanager.ms_spectra)):
    spectrum=specmanager.ms_spectra[spectrum_index];
    if 'FPTCount' in spectrum.parameters.keys():
        totalfpt+=spectrum.parameters['FPTCount'];
        
#%%
test_fpt_stat_outfile='e:/Imperial/TestDB/FPT2.txt';
print('Exporting FPT stats');
fout=open(test_fpt_stat_outfile,'w');
for i in range(11416):
    fout.write('%s\t%s\n'%(i,float(fpts[i])/totalfpt));
fout.close();

#%%
test_fpt_stat_outfile='e:/Imperial/TestDB/FPT3.txt';
print('Exporting FPT stats');
expcount=0;
fout=open(test_fpt_stat_outfile,'w');
for i in range(11416):
    if (float(fpts[i])/totalfpt)>=0.05 and (float(fpts[i])/totalfpt)<=0.95:
        fout.write('%s\t%s\n'%(i,float(fpts[i])/totalfpt));
        expcount+=1;
fout.close();
print(expcount)

#%%         
test_spectral_database_path='e:/Imperial/TestDB/NIST14_HMDB_MassBank_5fold_CV_FPT2';
specmanager=SpectralManager();
print('Reading test spectra...');
print(test_spectral_database_path);
specmanager.import_textfile_spectra_from_folder(test_spectral_database_path);
print('Finshed reading spectra. Total number: %s'%len(specmanager.ms_spectra));


#%%
ll=[];
print(len(specmanager.ms_spectra));
for spectrum_index in reversed(range(len(specmanager.ms_spectra))):
    spectrum=specmanager.ms_spectra[spectrum_index];
    if not ('fptcount' in spectrum.parameters.keys()):
        _=specmanager.ms_spectra.pop(spectrum_index);
        print('popped');
    else:
        
        neg=False;
        pos=False;
        
        for peak in spectrum.peaks:
            for subspectrum in peak.ms_spectra:
                if subspectrum.parameters['mode']==-1:
                    neg=True;
                elif subspectrum.parameters['mode']==1:
                    pos=True;
        
        if neg and pos:
            ll.append(spectrum_index);
            specneg=MSSpectrum();
            for key in spectrum.parameters.keys():
                specneg.parameters[key]=spectrum.parameters[key];
            specneg.parameters['mode']=-1;
            spectrum.parameters['mode']=1;            
            for i in reversed(range(len(spectrum.peaks))):
                if len(spectrum.peaks[i].ms_spectra)==0:
                    _=spectrum.peaks.pop(i);
                else:
                    move=False;
                    for spec in spectrum.peaks[i].ms_spectra:
                        if spec.parameters['mode']==-1:
                            move=True;
                            break;
                    if move:
                        specneg.peaks.append(spectrum.peaks.pop(i));
            specmanager.ms_spectra.append(specneg);
                
print(len(ll));
print(len(specmanager.ms_spectra));
#%%

print('Exporting spectra');
test_spectral_database_outpath='e:/Imperial/TestDB/NIST14_HMDB_MassBank_5fold_CV_FPT3';
specmanager.export_textfile_spectra_to_folder(test_spectral_database_outpath);         


#%%
         
specmanager.close();
dbmanager.close();
print('Finished');
     
