# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 12:13:08 2016

@author: Dr. Ivan Laponogov
"""

import sys;
if sys.byteorder!='little':
    print('Only little endian machines currently supported! bye bye ....');
    quit();

import os;
print(os.path.abspath(os.path.join(os.getcwd(),"..")));

sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"..")));

#import time;

import numpy as np;
#from chemdistiller.scorers.formula import FormulaScorer;
#from chemdistiller.scorers.element import ElementScorer;
#from chemdistiller.chemdb.manager import DBManager;
#from chemdistiller.filters.formula import FormulasFilter;
#from chemdistiller.filters.inchi import InChiFilter;
from chemdistiller.msspectra.manager import SpectralManager;
#from chemdistiller.msspectra.spectrum import MSSpectrum;
#from chemdistiller.utils.base64 import encode_to_base64;
from chemdistiller.utils.base64 import decode_from_base64;
#from chemdistiller.annotations.individualpeak import PeakAnnotator;
#from chemdistiller.scorers.totalscore import total_multiplicative_score;
#from chemdistiller.scorers.fragprint import FragPrintScorer;

trainSVM_spectral_database_path='e:/Imperial/TestDB/NIST14_HMDB_MassBank_5fold_CV_FPT3';

trainSVM_output_folder='e:/Imperial/TestDB/SVM_Train';

fpt_mask_file='e:/Imperial/TestDB/FPT3_0.05-0.95_2567.txt';

if not os.path.exists(trainSVM_output_folder):
    os.makedirs(trainSVM_output_folder);

nbatch=5;


print('Loading spectra....');
specmanager=SpectralManager();
specmanager.import_textfile_spectra_from_folder(trainSVM_spectral_database_path);

print('Finished reading spectra!');

#%%
print('Positive mode:');
maxmz=0.0;
for spectrum in specmanager.ms_spectra:
    if spectrum.parameters['mode']==1:
        if float(spectrum.parameters['exactmass'])>maxmz:
            maxmz=float(spectrum.parameters['exactmass']);
    
print(maxmz);

print('Negative mode:');
maxmz=0.0;
for spectrum in specmanager.ms_spectra:
    if spectrum.parameters['mode']==-1:
        if float(spectrum.parameters['exactmass'])>maxmz:
            maxmz=float(spectrum.parameters['exactmass']);
    
print(maxmz);

from operator import itemgetter;
adducts={};
missing=0;
for spectrum in specmanager.ms_spectra:
    for peak in spectrum.peaks:
        if 'ion_type' in peak.parameters:
            adducttype=peak.parameters['ion_type'];
            if adducttype in adducts.keys():
                adducts[adducttype]+=1;
            else:
                adducts[adducttype]=1;
        elif peak.ms_spectra:
            adducttype=peak.ms_spectra[0].parameters['precursor_ion'];
            if adducttype in adducts.keys():
                adducts[adducttype]+=1;
            else:
                adducts[adducttype]=1;
            
            
lst=[];            
for adducttype in adducts.keys():
    lst.append([adducttype,adducts[adducttype]]);

lst=sorted(lst,key=itemgetter(1),reverse=True);

for s in lst:
    print('%s\t%s'%(s[0],s[1]));
    
#print('Missing: %s'%missing);


#%%
positivefpt=[];
negativefpt=[];

positivespec=[];
negativespec=[];

positivebatch=[];
negativebatch=[];

scount=0;
for spectrum in specmanager.ms_spectra:
    scount+=1;
    print('%s of %s'%(scount,len(specmanager.ms_spectra)));
    fpts=[];
    fptcount=int(spectrum.parameters['fptcount']);
    ff=set();

    for i in range(fptcount):
        ff.add(spectrum.parameters['fpt_%s'%i]);

    for fpt in ff:    
        fpts.append(np.unpackbits(decode_from_base64(fpt)));
        
    batch=int(spectrum.parameters['crossvalidation_batch_index']);
    for peak in spectrum.peaks:
        ion_type='';
        if 'ion_type' in peak.parameters:
            ion_type=peak.parameters['ion_type'];
        elif peak.ms_spectra:
            ion_type=peak.ms_spectra[0].parameters['precursor_ion'];
        if ion_type=='[M+H]+' or ion_type=='[M-H]-':
            subcount=0;
            vector=np.zeros((20000,1),dtype=np.float32);
            values=[];
            for subspectrum in peak.ms_spectra:
                 if subspectrum.parameters['level']==2:
                     subspectrum.normalize_to_one();
                     subcount+=1;
                     for peak in subspectrum.peaks:
                         mzi=int(peak.mz*10);
                         vector[mzi]+=peak.intensity;
            if subcount>0:
                for i in range(0,20000):
                    if vector[i][0]>0:
                        values.append('%s:%s'%(i,vector[i][0]/subcount));
                values.append('20000:0');
                if spectrum.parameters['mode']==1:
                    positivefpt.append(fpts);
                    positivespec.append(values);
                    positivebatch.append(batch);
                elif spectrum.parameters['mode']==-1:
                    negativefpt.append(fpts);
                    negativespec.append(values);
                    negativebatch.append(batch);

print(len(positivefpt));
print(len(positivespec));
print(len(positivebatch));
print(len(negativefpt));
print(len(negativespec));
print(len(negativebatch));
#%%
fpt_mask=[];
finp=open(fpt_mask_file,'r');
for s in finp:
    s=s.rstrip('\n').split('\t');
    fpt_mask.append([int(s[0]),float(s[1])]);
finp.close();    
print(len(fpt_mask));
#%%
def generate_output(outpath, batch, mode, fpts, specs, batches):
    outpath=os.path.join(outpath,str(batch),mode);
    if not os.path.exists(outpath):
        os.makedirs(outpath);
    for index in range(len(fpt_mask)):
        print('%s, %s, %s of %s'%(mode, batch, index,len(fpt_mask)));
        fileset=set();
        for j in range(len(batches)):
            if batches[j]!=batch:
                for fpt in fpts[j]:
                    fvalue=fpt[fpt_mask[index][0]];
                    spectrum=specs[j];
                    ss='%s'%fvalue;
                    for value in spectrum:
                        ss+=' %s'%value;
                    fileset.add(ss);
        fout=open(os.path.join(outpath,'%s.model'%fpt_mask[index][0]),'w');
        for s in fileset:
            fout.write('%s\n'%s);
        fout.close();
#%%
#generate_output(trainSVM_output_folder, -1, 'Negative', negativefpt, negativespec, negativebatch);
#generate_output(trainSVM_output_folder, -1, 'Positive', positivefpt, positivespec, positivebatch);

#%%
    
#for batch in range(nbatch):
#    generate_output(trainSVM_output_folder, batch, 'Negative', negativefpt, negativespec, negativebatch);
#    generate_output(trainSVM_output_folder, batch, 'Positive', positivefpt, positivespec, positivebatch);
    
#%%
def generate_output_test(outpath, batch, mode, fpts, specs, batches):
    outpath=os.path.join(outpath,str(batch),mode);
    if not os.path.exists(outpath):
        os.makedirs(outpath);
    for index in range(len(fpt_mask)):
        print('%s, %s, %s of %s'%(mode, batch, index,len(fpt_mask)));
        #fileset=set();
        fout=open(os.path.join(outpath,'%s.test'%fpt_mask[index][0]),'w');
        for j in range(len(batches)):
            if batches[j]==batch:
                for fpt in fpts[j]:
                    fvalue=fpt[fpt_mask[index][0]];
                    spectrum=specs[j];
                    ss='%s'%fvalue;
                    for value in spectrum:
                        ss+=' %s'%value;
                    #fileset.add(ss);
                    fout.write('%s\n'%ss);
        fout.close();


#%%    

for batch in range(nbatch):
    generate_output_test(trainSVM_output_folder, batch, 'Negative', negativefpt, negativespec, negativebatch);
    generate_output_test(trainSVM_output_folder, batch, 'Positive', positivefpt, positivespec, positivebatch);

#%%    
    
def generate_scripts(outpath, batch, mode, fpts, batches):
    countfile=-1;
    outpath=os.path.join(outpath,str(batch),mode);
    if not os.path.exists(outpath):
        os.makedirs(outpath);
    for index in range(len(fpt_mask)):
        fptoccup=0.0;
        bitindex=fpt_mask[index][0];
        fptcount=0;
        fvalue=0;
        for j in range(len(batches)):
            if batches[j]!=batch:
                fptcount+=len(fpts[j]);
                for fpt in fpts[j]:
                    fvalue+=fpt[bitindex];
        if fptcount>0:
            fptoccup=float(fvalue)/fptcount;

        
        if fptoccup>=0.05 and fptoccup<=0.95:
            datafile='%s.model'%fpt_mask[index][0];
            c1=1/(1-fptoccup);
            c2=1/(fptoccup);
            countfile+=1;
            print('%s, %s, %s of %s'%(mode, batch, index,len(fpt_mask)));
            fout=open(os.path.join(outpath,'processlinear_%s.bat'%countfile),'w');
            fout.write('#!/bin/bash\n');
            fout.write('#Occupancy of bit %s is %s\n'%(bitindex,fptoccup));
            
            
            
            fout.write('cp $WORK/$JOBDIR/%s.gz $TMPDIR/$PBS_ARRAY_INDEX/\n'%datafile);
            fout.write('gunzip %s.gz\n'%datafile);
    
            fout.write('echo Linear Classifier Weighted >linear_weighted.out\n');
            fout.write('echo Linear Classifier Weighted\n');
            for i in range(-20,20,3):
                fout.write('echo %s\n'%(i));
                fout.write('echo Setting -c=%s >>linear_weighted.out\n\n'%(2**i));
                fout.write('echo ./svm-train -t 0 -c %s -v 5 -w0 %s -w1 %s %s >>linear_weighted.out\n'%(2**i,c1,c2,datafile));
                fout.write('./svm-train -t 0 -c %s -v 5 -w0 %s -w1 %s %s >>linear_weighted.out\n'%(2**i,c1,c2,datafile));
                
    
            fout.write('echo Linear Classifier NotWeighted >linear_notweighted.out\n');
            fout.write('echo Linear Classifier NotWeighted\n');
            for i in range(-20,20,3):
                fout.write('echo %s\n'%(i));
                fout.write('echo Setting -c=%s >>linear_notweighted.out\n\n'%(2**i));
                fout.write('echo ./svm-train -t 0 -c %s -v 5 %s >>linear_notweighted.out\n'%(2**i,datafile));
                fout.write('./svm-train -t 0 -c %s -v 5 %s >>linear_notweighted.out\n'%(2**i,datafile));
    
            fout.write('python parsesvmoutput_linearrefine.py\n');
            fout.write('source ./do_models_linearrefine.bat\n');
    
            fout.write('python parsesvmoutput_linear.py\n');
            fout.write('source ./do_models_linear.bat\n');
            fout.write('gzip %s.linear_weight\n'%datafile);
            fout.write('gzip %s.linear_noweight\n'%datafile);
            fout.write('gzip linear_notweighted.out\n');
            fout.write('gzip linear_weighted.out\n');        
            
            fout.write('cp do_models_linear.bat $WORK/$JOBDIR/%s.linear.do.bat\n'%datafile);
            fout.write('cp %s.linear_weight.gz $WORK/$JOBDIR/\n'%datafile);
            fout.write('cp %s.linear_noweight.gz $WORK/$JOBDIR/\n'%datafile);
            fout.write('cp linear_notweighted.out.gz $WORK/$JOBDIR/%s.linear_notweighted.out.gz\n'%datafile);
            fout.write('cp linear_weighted.out.gz $WORK/$JOBDIR/%s.linear_weighted.out.gz\n'%datafile);
    
            fout.close();
    
            fout=open(os.path.join(outpath,'processnonlinear_%s.bat'%countfile),'w');
            fout.write('#!/bin/bash\n');
            fout.write('#Occupancy of bit %s is %s\n'%(bitindex,fptoccup));
            fout.write('echo Radial Basis Classifier Weighted >radial_weighted.out\n');
            fout.write('echo Radial Basis Classifier Weighted\n');
            
            fout.write('cp $WORK/$JOBDIR/%s.gz $TMPDIR/$PBS_ARRAY_INDEX/\n'%datafile);
            fout.write('gunzip %s.gz\n'%datafile);
    
            
            for i in range(-20,20,3):
                for j in range(-20,20,3):
                    fout.write('echo %s %s\n'%(i,j));
                    fout.write('echo Setting -c=%s -g=%s >>radial_weighted.out\n\n'%(2**i,2**j));
                    fout.write('echo ./svm-train -t 2 -c %s -g %s -v 5 -w0 %s -w1 %s %s >>radial_weighted.out\n'%(2**i,2**j,c1,c2,datafile));
                    fout.write('./svm-train -t 2 -c %s -g %s -v 5 -w0 %s -w1 %s %s >>radial_weighted.out\n'%(2**i,2**j,c1,c2,datafile));
            
            fout.write('echo Radial Basis Classifier NotWeighted >radial_notweighted.out\n');
            fout.write('echo Radial Basis Classifier NotWeighted\n');
            for i in range(-20,20,3):
                for j in range(-20,20,3):
                    fout.write('echo %s %s\n'%(i,j));
                    fout.write('echo Setting -c=%s -g=%s >>radial_notweighted.out\n\n'%(2**i,2**j));
                    fout.write('echo ./svm-train -t 2 -c %s -g %s -v 5 %s >>radial_notweighted.out\n'%(2**i,2**j,datafile));
                    fout.write('./svm-train -t 2 -c %s -g %s -v 5 %s >>radial_notweighted.out\n'%(2**i,2**j,datafile));
    
            fout.write('python parsesvmoutput_nonlinearrefine.py\n');
            fout.write('source ./do_models_nonlinearrefine.bat\n');
    
                    
            fout.write('python parsesvmoutput_nonlinear.py\n');
            fout.write('source ./do_models_nonlinear.bat\n');
            fout.write('gzip %s.radial_weight\n'%datafile);
            fout.write('gzip %s.radial_noweight\n'%datafile);
            fout.write('gzip radial_notweighted.out\n');
            fout.write('gzip radial_weighted.out\n');   
            
            fout.write('cp do_models_nonlinear.bat $WORK/$JOBDIR/%s.radial.do.bat\n'%datafile);
            fout.write('cp %s.radial_weight.gz $WORK/$JOBDIR/\n'%datafile);
            fout.write('cp %s.radial_noweight.gz $WORK/$JOBDIR/\n'%datafile);
            fout.write('cp radial_notweighted.out.gz $WORK/$JOBDIR/%s.radial_notweighted.out.gz\n'%datafile);
            fout.write('cp radial_weighted.out.gz $WORK/$JOBDIR/%s.radial_weighted.out.gz\n'%datafile);
                    
            fout.close();
            #break;

#%%
#generate_scripts(trainSVM_output_folder, -1, 'Negative', negativefpt, negativebatch);
#generate_scripts(trainSVM_output_folder, -1, 'Positive', positivefpt, positivebatch);


#%%

        
#for batch in range(nbatch):
#    generate_scripts(trainSVM_output_folder, batch, 'Negative', negativefpt, negativebatch);
#    generate_scripts(trainSVM_output_folder, batch, 'Positive', positivefpt, positivebatch);
    
    
    
#%%

def generate_scripts_test(outpath, batch, mode, fpts, batches):
    countfile=-1;
    outpath=os.path.join(outpath,str(batch),mode);
    if not os.path.exists(outpath):
        os.makedirs(outpath);
    for index in range(len(fpt_mask)):
        fptoccup=0.0;
        bitindex=fpt_mask[index][0];
        fptcount=0;
        fvalue=0;
        for j in range(len(batches)):
            if batches[j]==batch:
                fptcount+=len(fpts[j]);
                for fpt in fpts[j]:
                    fvalue+=fpt[bitindex];
        if fptcount>0:
            fptoccup=float(fvalue)/fptcount;

        
        
            datafile='%s.test'%fpt_mask[index][0];
            
            c1=1/(1-fptoccup);
            c2=1/(fptoccup);
            countfile+=1;
            print('%s, %s, %s of %s'%(mode, batch, index,len(fpt_mask)));
            
            fout=open(os.path.join(outpath,'testlinear_%s.bat'%countfile),'w');
            fout.write('#!/bin/bash\n');
            
            fout.write('cp $WORK/$JOBDIR/%s.gz $TMPDIR/$PBS_ARRAY_INDEX/\n'%datafile);
            fout.write('gunzip %s.gz\n'%datafile);
            
            modelfile='%s.test'%fpt_mask[index][0];
    
            fout.write('./svm-predict %s %s %s >>%s.log\n'%(datafile));
            
            
                
    
            fout.write('echo Linear Classifier NotWeighted >linear_notweighted.out\n');
            fout.write('echo Linear Classifier NotWeighted\n');
            for i in range(-20,20,3):
                fout.write('echo %s\n'%(i));
                fout.write('echo Setting -c=%s >>linear_notweighted.out\n\n'%(2**i));
                fout.write('echo ./svm-train -t 0 -c %s -v 5 %s >>linear_notweighted.out\n'%(2**i,datafile));
                fout.write('./svm-train -t 0 -c %s -v 5 %s >>linear_notweighted.out\n'%(2**i,datafile));
    
            fout.write('python parsesvmoutput_linearrefine.py\n');
            fout.write('source ./do_models_linearrefine.bat\n');
    
            fout.write('python parsesvmoutput_linear.py\n');
            fout.write('source ./do_models_linear.bat\n');
            fout.write('gzip %s.linear_weight\n'%datafile);
            fout.write('gzip %s.linear_noweight\n'%datafile);
            fout.write('gzip linear_notweighted.out\n');
            fout.write('gzip linear_weighted.out\n');        
            
            fout.write('cp do_models_linear.bat $WORK/$JOBDIR/%s.linear.do.bat\n'%datafile);
            fout.write('cp %s.linear_weight.gz $WORK/$JOBDIR/\n'%datafile);
            fout.write('cp %s.linear_noweight.gz $WORK/$JOBDIR/\n'%datafile);
            fout.write('cp linear_notweighted.out.gz $WORK/$JOBDIR/%s.linear_notweighted.out.gz\n'%datafile);
            fout.write('cp linear_weighted.out.gz $WORK/$JOBDIR/%s.linear_weighted.out.gz\n'%datafile);
    
            fout.close();
    
            fout=open(os.path.join(outpath,'processnonlinear_%s.bat'%countfile),'w');
            fout.write('#!/bin/bash\n');
            fout.write('#Occupancy of bit %s is %s\n'%(bitindex,fptoccup));
            fout.write('echo Radial Basis Classifier Weighted >radial_weighted.out\n');
            fout.write('echo Radial Basis Classifier Weighted\n');
            
            fout.write('cp $WORK/$JOBDIR/%s.gz $TMPDIR/$PBS_ARRAY_INDEX/\n'%datafile);
            fout.write('gunzip %s.gz\n'%datafile);
    
            
            for i in range(-20,20,3):
                for j in range(-20,20,3):
                    fout.write('echo %s %s\n'%(i,j));
                    fout.write('echo Setting -c=%s -g=%s >>radial_weighted.out\n\n'%(2**i,2**j));
                    fout.write('echo ./svm-train -t 2 -c %s -g %s -v 5 -w0 %s -w1 %s %s >>radial_weighted.out\n'%(2**i,2**j,c1,c2,datafile));
                    fout.write('./svm-train -t 2 -c %s -g %s -v 5 -w0 %s -w1 %s %s >>radial_weighted.out\n'%(2**i,2**j,c1,c2,datafile));
            
            fout.write('echo Radial Basis Classifier NotWeighted >radial_notweighted.out\n');
            fout.write('echo Radial Basis Classifier NotWeighted\n');
            for i in range(-20,20,3):
                for j in range(-20,20,3):
                    fout.write('echo %s %s\n'%(i,j));
                    fout.write('echo Setting -c=%s -g=%s >>radial_notweighted.out\n\n'%(2**i,2**j));
                    fout.write('echo ./svm-train -t 2 -c %s -g %s -v 5 %s >>radial_notweighted.out\n'%(2**i,2**j,datafile));
                    fout.write('./svm-train -t 2 -c %s -g %s -v 5 %s >>radial_notweighted.out\n'%(2**i,2**j,datafile));
    
            fout.write('python parsesvmoutput_nonlinearrefine.py\n');
            fout.write('source ./do_models_nonlinearrefine.bat\n');
    
                    
            fout.write('python parsesvmoutput_nonlinear.py\n');
            fout.write('source ./do_models_nonlinear.bat\n');
            fout.write('gzip %s.radial_weight\n'%datafile);
            fout.write('gzip %s.radial_noweight\n'%datafile);
            fout.write('gzip radial_notweighted.out\n');
            fout.write('gzip radial_weighted.out\n');   
            
            fout.write('cp do_models_nonlinear.bat $WORK/$JOBDIR/%s.radial.do.bat\n'%datafile);
            fout.write('cp %s.radial_weight.gz $WORK/$JOBDIR/\n'%datafile);
            fout.write('cp %s.radial_noweight.gz $WORK/$JOBDIR/\n'%datafile);
            fout.write('cp radial_notweighted.out.gz $WORK/$JOBDIR/%s.radial_notweighted.out.gz\n'%datafile);
            fout.write('cp radial_weighted.out.gz $WORK/$JOBDIR/%s.radial_weighted.out.gz\n'%datafile);
                    
            fout.close();
    
for batch in range(nbatch):
    generate_scripts_test(trainSVM_output_folder, batch, 'Negative', negativefpt, negativebatch);
    generate_scripts_test(trainSVM_output_folder, batch, 'Positive', positivefpt, positivebatch);
    
    
    
    
    