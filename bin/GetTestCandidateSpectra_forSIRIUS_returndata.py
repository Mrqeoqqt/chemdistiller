# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 11:34:09 2017

@author: ilaponog
"""

import sys;
import os;
import numpy as np;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."));
    sys.path.append(chemdistiller_path);


import time;

from chemdistiller.msspectra.manager import SpectralManager;
from chemdistiller.msspectra.adducts import get_adduct_by_name;
from chemdistiller.scorers.fragprint import FragPrintScorer;
from chemdistiller.chemdb.manager import DBManager;
from chemdistiller.utils.inchi import inchikey_from_inchi;
from chemdistiller.utils.base64 import decode_from_base64, encode_to_base64;


from operator import itemgetter;




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
    
def get_mean_mz_weighted(spectrum):
    sum_mz=0.0;
    sum_i=0.0;
    for peak in spectrum.peaks:
        sum_mz+=peak.mz*peak.intensity;
        sum_i+=peak.intensity;
    if sum_i>0.0:
        spectrum.mean_mz_weighted=sum_mz/sum_i;


def find_and_store_candidates(dbpath, mass, dbmanager, dbindexes):
    print(mass);
    result=dbmanager.query_by_mz_scored(mass, 20.0, 0, dbindexes, required_fields=set(['InChI','IDs','InChIKey']), results_limit=0, save_memory=False, revisitable=False);
    print('%s candidates'%len(result.mol_list));
    with open(dbpath,'w') as fout:
        fout.write('MonoisotopicMass|InChI|Identifier|InChIKey2|InChIKey1|MolecularFormula\n');
        for candidate in result.mol_list:
            key=candidate['InChIKey'].split('-');
            formula=candidate['InChI'].split('/')[0];
            fout.write('%s|InChI=1S/%s|%s|%s|%s|%s\n'%(candidate['Mass'], candidate['InChI'], candidate['IDs'], key[1], key[0], formula));
    


def save_inchi(peaklist, outdir, mode):
    
    if mode==1:
        outpath='%s/positive'%outdir;
    else:
        outpath='%s/negative'%outdir;
    
    if not os.path.exists(outpath):
        os.makedirs(outpath);
        
    print(outpath);    
    for i in range(len(peaklist)):
        peak=peaklist[i];
        with open('%s/%s.inchi'%(outpath,i),'w') as fout:
            inchi=peak.parent_spectrum.parameters['inchi'];
            inchikey=inchikey_from_inchi(inchi);
            fout.write('%s\t%s\n'%(inchi, inchikey))
            
    
    
def process_output(peaklist, outdir, mode, dbmanager, dbindexes):
    
    if mode==1:
        outpath='%s/positive'%outdir;
    else:
        outpath='%s/negative'%outdir;
    
    if not os.path.exists(outpath):
        os.makedirs(outpath);
        
    print(outpath);
    
    fout=open('%s/process.bat'%outpath,'w');
    for i in range(len(peaklist)):
        print('%s of %s'%(i, len(peaklist)))
        peak=peaklist[i];
        mass=float(peak.parent_spectrum.parameters['exactmass']);
        fout.write('java -jar G:\\MetFrag\\MetFrag2.3-CL.jar params_%s.txt\n'%i);
        dbpath='%s/db_%s.psv'%(outpath, i);
        with open('%s/params_%s.txt'%(outpath, i), 'w') as fparam:
            fparam.write('#\n');
            fparam.write('# data file containing mz intensity peak pairs (one per line)\n');
            fparam.write('#\n');
            fparam.write('PeakListPath = %s.txt\n'%i);
            fparam.write('\n');
            fparam.write('#\n');
            fparam.write('# database parameters -> how to retrieve candidates\n');
            fparam.write('# \n');
            fparam.write('#\n');
            fparam.write('MetFragDatabaseType = LocalPSV\n');
            fparam.write('LocalDatabasePath=%s\n'%dbpath);
            fparam.write('\n');
            fparam.write('NeutralPrecursorMass = %.5f\n'%mass);
            fparam.write('#IonizedPrecursorMass = 349.93356\n');
            fparam.write('\n');
            fparam.write('#\n');
            fparam.write('# peak matching parameters\n');
            fparam.write('#\n');
            fparam.write('FragmentPeakMatchAbsoluteMassDeviation = 0.005\n');
            fparam.write('FragmentPeakMatchRelativeMassDeviation = 20\n');
            if mode==1:
                fparam.write('PrecursorIonMode = 1\n');
                fparam.write('IsPositiveIonMode = True\n');
            else:
                fparam.write('PrecursorIonMode = -1\n');
                fparam.write('IsPositiveIonMode = False\n');
                
            fparam.write('#\n');
            fparam.write('# scoring parameters\n');
            fparam.write('#\n');
            fparam.write('MetFragScoreTypes = FragmenterScore\n');
            fparam.write('MetFragScoreWeights = 1.0\n');
            fparam.write('\n');
            fparam.write('#\n');
            fparam.write('# output\n');
            fparam.write('# SDF, XLS, CSV, ExtendedXLS, ExtendedFragmentsXLS\n');
            fparam.write('#\n');
            fparam.write('MetFragCandidateWriter = CSV\n');
            fparam.write('SampleName = %s\n'%i);
            fparam.write('ResultsPath = results\n');
            fparam.write('#\n');
            fparam.write('# following parameteres can be kept as they are\n');
            fparam.write('#\n');
            fparam.write('MaximumTreeDepth = 2\n');
            fparam.write('MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter\n');
            fparam.write('MetFragPostProcessingCandidateFilter = InChIKeyFilter\n');
            fparam.write('# NumberThreads = 1\n');
    
        with open('%s/%s.txt'%(outpath, i),'w') as fpeaks:
            for j in peak.merged_spectrum:
                fpeaks.write('%s\t%s\n'%(j[0],j[1]));
        find_and_store_candidates(dbpath, mass, dbmanager, dbindexes);
    
    fout.close();
    
    

def normalize_and_merge_spectra(speclist, precision):
    result=[];
    for i in range(len(speclist)):
        spec=speclist[i];
        spec.normalize_to_one();
        for peak in spec.peaks:
            result.append([peak.mz, peak.intensity/len(speclist)]);
        result=sorted(result,key=itemgetter(0));
    
    FragPrintScorer.condense_peak_list(result, precision, 1);
    
    return result;


def generate_fpt(index, peak, subpath, mode, mask):
    fptcount=int(peak.parent_spectrum.parameters['fptcount']);
    fpt=np.zeros((len(mask),), dtype=np.float32);
    for i in range(fptcount):
        f=peak.parent_spectrum.parameters['fpt_%s'%i];
        f=decode_from_base64(f);
        f=np.array(np.unpackbits(f), dtype=np.float32);
        fpt=np.add(fpt,f);
    fpt=np.divide(fpt, fptcount);
    fpt=np.rint(fpt);
    fpt=np.subtract(np.multiply(fpt, 2),1);
    
    batch=int(peak.parent_spectrum.parameters['crossvalidation_batch_index']);
    subpath=subpath+'/%s'%batch;
    if not os.path.exists(subpath):
        os.makedirs(subpath);
    fname=subpath+'/%s.fpt'%index;
    
    with open(fname,'w') as fout:
        for i in range(len(mask)):
            if mask[i]>0:
                fout.write('%s '%fpt[i]);
    
def generate_list(index, peak, subpath, mode):
    batch=int(peak.parent_spectrum.parameters['crossvalidation_batch_index']);
    if not os.path.exists(subpath):
        os.makedirs(subpath);
    fname=subpath+'/ids.txt';
    with open(fname,'a') as fout:
        fout.write('%s %s\n'%(index, batch));
        

    
    
    
def generate_mgf(index, peak, subpath, mode):
    batch=int(peak.parent_spectrum.parameters['crossvalidation_batch_index']);
    #if batch!=0:
    #    return
    subpath=subpath+'/%s'%batch;
    if not os.path.exists(subpath):
        os.makedirs(subpath);
    fname=subpath+'/%s.mgf'%index;
    
    
    formula=peak.parent_spectrum.parameters['inchi'].split('/')[0];
    
    with open(fname,'w') as fout:
        fout.write('BEGIN IONS\n');
        fout.write('PEPMASS=%s\n'%peak.mz);
        fout.write('CHARGE=%s\n'%mode);
        fout.write('MSLEVEL=2\n');
        if mode>0:
            fout.write('IONMODE=positive\n');
        else:
            fout.write('IONMODE=negative\n');
        if peak.parameters['ion_type']=='[M+H]+':
            adduct='M+H';
        elif  peak.parameters['ion_type']=='[M-H]-':
            adduct='M-H';
            
        fout.write('NAME=%s %s\n'%('compound',adduct));
        
        fout.write('INCHI=InChI=1S/%s\n'%peak.parent_spectrum.parameters['inchi']);
        for speak in peak.merged_spectrum:
            fout.write('%s	%s\n'%(speak[0], speak[1]));
        fout.write('END IONS\n');
    
    fname2=subpath+'/process.bat';
    
    if batch==0:
        ff='';
        for i in formula:
            if i<'0' or i>'9':
                ff+=i;
        with open(fname2,'a') as fout:
            fout.write('e:/sirius3-win64-3.4.1/sirius3-console-64.exe --ppm-max 20 --elements CHNOPSFIClBrSi--output %s/%s %s\n'%( subpath, index, fname))
        
    else:
        with open(fname2,'a') as fout:
                fout.write('e:/sirius3-win64-3.4.1/sirius3-console-64.exe --ppm-max 20 --formula %s %s --output %s/%s\n'%(formula, fname, subpath, index))

def loadresults(peaks, data_path):
    fname = data_path + '/fpt_mask.default';
    with open(fname,'rb') as finp:
        fpt_mask = np.fromfile(finp, dtype=np.uint8);
        print(fpt_mask.shape);
        fpt_mask = encode_to_base64(fpt_mask);    
        print(len(fpt_mask))
        print(fpt_mask)
        time.sleep(15)
    
    
    
    pred_test = np.loadtxt(data_path + '/pred_test.txt')
    pred_train = np.loadtxt(data_path + '/pred_train.txt')

    inchi_test = [];
    formulas_test = [];
    scores_test = [];

    inchi_train = [];
    formulas_train = [];
    scores_train = [];
    
    fname = data_path + '/train_inchi.txt';
    with open(fname,'r') as finp:
        for s in finp:
            inchi_train.append(inchikey_from_inchi(s.rstrip('\n')).split('-')[0]);
            
    fname = data_path + '/train_formulas.txt';
    with open(fname,'r') as finp:
        for s in finp:
            formulas_train.append(s.rstrip('\n'));

    fname = data_path + '/train_formula_scores.txt';
    with open(fname,'r') as finp:
        for s in finp:
            scores_train.append(float(s.rstrip('\n')));



    fname = data_path + '/test_inchi.txt';
    with open(fname,'r') as finp:
        for s in finp:
            inchi_test.append(inchikey_from_inchi(s.rstrip('\n')).split('-')[0]);
            
    fname = data_path + '/test_formulas.txt';
    with open(fname,'r') as finp:
        for s in finp:
            formulas_test.append(s.rstrip('\n'));

    fname = data_path + '/test_formula_scores.txt';
    with open(fname,'r') as finp:
        for s in finp:
            scores_test.append(float(s.rstrip('\n')));

    print(len(pred_test))
    print(len(inchi_test))
    print(len(formulas_test))
    print(len(scores_test))

    print(len(pred_train))
    print(len(inchi_train))
    print(len(formulas_train))
    print(len(scores_train))
    
    if len(pred_test)!=len(inchi_test) or len(inchi_test) != len(formulas_test) or len(formulas_test) != len(scores_test):
        raise ValueError('Unequal test sizes!');

    if len(pred_train)!=len(inchi_train) or len(inchi_train) != len(formulas_train) or len(formulas_train) != len(scores_train):
        raise ValueError('Unequal train sizes!');


    inchiindex = {};
    
    for i in range(len(peaks)):
        pi = peaks[i].parent_spectrum.parameters['inchi'];
        pi = inchikey_from_inchi(pi).split('-')[0]
        inchiindex[pi] = i;
    
    #with open(data_path + '/listinchi.txt', 'w') as fout:
    #    for key in inchiindex.keys():
    #        fout.write('%s\n'%key)
    

    inchi = inchi_test + inchi_train
    formulas = formulas_test + formulas_train
    scores = scores_test + scores_train
    pred = np.vstack((pred_test, pred_train))
    
    lost_count = 0;    
    for i in range(len(inchi)):
        if inchi[i] in inchiindex:
            j = inchiindex[inchi[i]];
            peak = peaks[j];
            if 'csifingerid_count' in peak.parameters:
                peak.parameters['csifingerid_count'] += 1;
            else:
                peak.parameters['csifingerid_count'] = 1;
            
            subi = peak.parameters['csifingerid_count'];
            p = pred[i, :].tolist();
            for k in range(len(p)):
                p[k] = str(p[k])
            peak.parameters['csifingerid_predfpt_%s'%subi] = ','.join(p);
            peak.parameters['csifingerid_fptmask_%s'%subi] = fpt_mask;
            peak.parameters['csifingerid_formula_%s'%subi] = formulas[i];
            peak.parameters['csifingerid_score_%s'%subi] = scores[i];
            
        else:
            print('Error! Inchi not found! %s'%inchi[i]);
            lost_count += 1;
    print('Total Missing: %s'%lost_count)        
    
    
    notpresent_count = len(peaks);
    for peak in peaks:
        if 'csifingerid_count' in peak.parameters:
            notpresent_count -= 1;
    
    
    return lost_count, notpresent_count

if __name__=='__main__':
    
    #test_chemical_databases_large=['PubChem','MassBank','HMDB','ChEBI'];
    #test_chemical_databases_small=['TestDB','MassBank','HMDB','ChEBI'];
    
    test_spectral_database_path='e:/Imperial/TestDB/NIST14_HMDB_MassBank_5fold_CV_FPT3';
    outpath='i:/CSI_FingerID_20ppm_elements_testready_withbothtesttrain_new';
    
    result_path = 'e:/Linux/Exchange/CSI_FingerID/examples/ChemDistillerData/kernels/';
    
    #fname=os.path.join(outpath,'negative','fpt_mask.default');
    #with open(fname,'rb') as finp:
    #    negative_fpt_mask=np.unpackbits(np.fromfile(finp, dtype=np.uint8));    

    #fname=os.path.join(outpath,'positive','fpt_mask.default');
    #with open(fname,'rb') as finp:
    #    positive_fpt_mask=np.unpackbits(np.fromfile(finp, dtype=np.uint8));    
        
    
    ppm=21;

    print('Starting. %s'%time.strftime('%d/%m/%y %H:%M:%S'));
    starttime=time.time();

    precision=0.005;
    
    specmanager=SpectralManager();
        
    specmanager.import_textfile_spectra_from_folder(test_spectral_database_path);
    print('Reading complete')
        
    remove_low_precision_spectra(precision);
    print('Removed low precision spectra')
    correct_precursor_mz(specmanager);   
    print('Corrected MZ precursor')

    #print(positive_fpt_mask.shape)    
    #print(negative_fpt_mask.shape)    
    
    peaks=[];
    #dbfile=os.path.abspath('../tests/databases_HDF5.list');
    #print(dbfile);
    #dbmanager=DBManager(dbfile);
    #db_indexes_large=dbmanager.db_indexes_from_db_names(test_chemical_databases_large,case_sensitive=True);
    #db_indexes_small=dbmanager.db_indexes_from_db_names(test_chemical_databases_small,case_sensitive=True);
    na = 0;
    i = 0;
    
    for spectrum in specmanager.ms_spectra:
        if spectrum.parameters['charge']=='0':
                #if int(spectrum.parameters['crossvalidation_batch_index'])==0:
                for peak in spectrum.peaks:
                        if hasattr(peak,'parameters'):
                            if 'ion_type' in peak.parameters:
                               if (peak.parameters['ion_type'] in FragPrintScorer.supported_adducts):
                                   peaks.append(peak);
                    
    for i in reversed(range(len(peaks))):
        peak=peaks[i];
        if len(peak.ms_spectra)<1:
            del peaks[i];
        else:
            peak.merged_spectrum=normalize_and_merge_spectra(peak.ms_spectra, precision);
    neg=[];
    pos=[];
    for peak in peaks:
        if peak.parent_spectrum.parameters['mode']==-1:
            neg.append(peak);
        else:
            pos.append(peak);
    
    print('Positive mode: %s'%len(pos));
    print('Negative mode: %s'%len(neg));    
    
    ii, nna = loadresults(pos, result_path+'positive')    
    
    
    i += ii;
    na += nna;
    
    ii, nna = loadresults(neg, result_path+'negative')

    i += ii;
    na += nna;
    
    #for i in range(len(pos)):
    #    generate_mgf(i, pos[i], outpath+'/positive', 1);
    #    generate_fpt(i, pos[i], outpath+'/positive', 1, positive_fpt_mask);
    #    generate_list(i, pos[i], outpath+'/positive', 1);

    #for i in range(len(neg)):
    #    generate_mgf(i, neg[i], outpath+'/negative', -1);
    #    generate_fpt(i, neg[i], outpath+'/negative', -1, negative_fpt_mask);
    #    generate_list(i, neg[i], outpath+'/negative', -1);
    
    
    
    #process_output(neg,'g:/MetFrag/smallDB0',-1, dbmanager, db_indexes_small);
    #process_output(pos,'g:/MetFrag/smallDB0',1, dbmanager, db_indexes_small);
    
    #process_output(neg,'g:/MetFrag/largeDB0',-1, dbmanager, db_indexes_large);

    #process_output(pos,'g:/MetFrag/largeDB0',1, dbmanager, db_indexes_large);

    #save_inchi(neg,'g:/MetFrag/smallDB0',-1);
    #save_inchi(pos,'g:/MetFrag/smallDB0',1);
    #save_inchi(neg,'g:/MetFrag/largeDB0',-1);
    #save_inchi(pos,'g:/MetFrag/largeDB0',1);
    
    
    #dbmanager.close();
    specmanager.export_textfile_spectra_to_folder(outpath)
    specmanager.close();
    print('Lost :%s'%i)
    print('Peaks without CSI:FingerID: %s'%na)
    
        