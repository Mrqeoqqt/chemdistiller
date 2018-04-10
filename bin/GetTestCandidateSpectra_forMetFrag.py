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

from chemdistiller.msspectra.manager import SpectralManager;
from chemdistiller.msspectra.adducts import get_adduct_by_name;
from chemdistiller.scorers.fragprint import FragPrintScorer;
from chemdistiller.chemdb.manager import DBManager;
from chemdistiller.utils.inchi import inchikey_from_inchi;


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
    

if __name__=='__main__':
    
    test_chemical_databases_large=['PubChem','MassBank','HMDB','ChEBI'];
    test_chemical_databases_small=['TestDB','MassBank','HMDB','ChEBI'];
    
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
    dbfile=os.path.abspath('../tests/databases_HDF5.list');
    print(dbfile);
    dbmanager=DBManager(dbfile);
    db_indexes_large=dbmanager.db_indexes_from_db_names(test_chemical_databases_large,case_sensitive=True);
    db_indexes_small=dbmanager.db_indexes_from_db_names(test_chemical_databases_small,case_sensitive=True);
    
    
    for spectrum in specmanager.ms_spectra:
        if spectrum.parameters['charge']=='0':
            if int(spectrum.parameters['crossvalidation_batch_index'])==0:
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
    #process_output(neg,'g:/MetFrag/smallDB0',-1, dbmanager, db_indexes_small);
    #process_output(pos,'g:/MetFrag/smallDB0',1, dbmanager, db_indexes_small);
    
    #process_output(neg,'g:/MetFrag/largeDB0',-1, dbmanager, db_indexes_large);

    #process_output(pos,'g:/MetFrag/largeDB0',1, dbmanager, db_indexes_large);

    save_inchi(neg,'g:/MetFrag/smallDB0',-1);
    save_inchi(pos,'g:/MetFrag/smallDB0',1);
    save_inchi(neg,'g:/MetFrag/largeDB0',-1);
    save_inchi(pos,'g:/MetFrag/largeDB0',1);
    
    
    dbmanager.close();
                 
    specmanager.close();
    
        