# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 15:03:03 2016

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
    
import numpy as np;
from chemdistiller.settings import test_as_single_process, registered_scorers, registered_filters;
from chemdistiller.scorers.svm import SVM;
from chemdistiller.msspectra.adducts import global_supported_adducts;
from chemdistiller.utils.sysutils import print_in_the_same_line;


if not test_as_single_process:
    import multiprocessing;


class FingerPrintScorer:
    supported_adducts=SVM.supported_adducts;
    supported_isotopes=SVM.supported_isotopes;
    required_fields=['FPT'];
    name='FingerPrintScorer';
    
    def __init__(self, process_id, settings):
        #self.supported_adducts=FingerPrintScorer.supported_adducts;
        #self.required_fields=FingerPrintScorer.required_fields;
        #self.supported_isotopes=FingerPrintScorer.supported_isotopes;
        self.settings=settings;
        self.process_id=process_id;

    @staticmethod
    def __generate_fpt_vectors(peak, overwrite=False):
        if overwrite==False:
            if hasattr(peak,'fpt_sparse_vector'):
                return peak.fpt_sparse_vector;
        
        subcount=0;
        vector=np.zeros((20000,1),dtype=np.float32);
        values=[];
        for subspectrum in peak.ms_spectra:
            if subspectrum.parameters['level']==2:
                 subspectrum.normalize_to_one();
                 subcount+=1;
                 for subpeak in subspectrum.peaks:
                     mzi=int(subpeak.mz*10);
                     vector[mzi]+=subpeak.intensity;
        if subcount>0:
            for i in range(0,20000):
                if vector[i][0]>0:
                    values.append((i,vector[i][0]/subcount));
        
        values.append((20000,0.0));
                 
        return values;
        

    @staticmethod
    def tidy_up_peaks(peaks, required_fields):
        for peak in peaks:
            if hasattr(peak,'fpt_sparse_vector'):
                del peak.fpt_sparse_vector;
        if not ('FPT') in required_fields:
            for peak in peaks:
                for annotation in peak.annotations:
                    if hasattr(annotation,'fpt_mask'):
                        del annotation.fpt_mask;
                    if hasattr(annotation,'predicted_fpt'):
                        del annotation.predicted_fpt;
                    if not (annotation.mol_candidates is None):
                        for candidate in annotation.mol_candidates.mol_list:
                            try:
                                del candidate['FPT'];
                            except:
                                print('No "FPT" found');

    @staticmethod
    def predict_fingerprints(peaks, fingerprint_model_path, settings,  nproc=1, overwrite=False):
        global test_as_single_process;
        annotations={};
        print('\n');
        i=0;
        maxi=len(peaks);
        for peak in peaks:
            i+=1;
            print_in_the_same_line('\rPeak %s of %s'%(i,maxi));
            if peak.annotations:
                peak.fpt_sparse_vector=FingerPrintScorer.__generate_fpt_vectors(peak, overwrite);
                for annotation in peak.annotations:
                    if annotation.adduct.definition in FingerPrintScorer.supported_adducts:
                        if annotation.isotope in FingerPrintScorer.supported_isotopes:
                            if (overwrite==True) or (not hasattr(annotation,'predicted_fpt')) or (not hasattr(annotation,'fpt_mask')):
                                desc='%s#%s'%(annotation.adduct.definition, annotation.isotope);
                                if desc in annotations:
                                    annotations[desc].append(annotation);
                                else:
                                    annotations[desc]=[annotation];
        print('\n');
        for desc in annotations.keys():
            svm=SVM(fingerprint_model_path, settings, nproc);
            svm.process_annotations(annotations[desc]);



        
    def close(self):
        return 0;


    def configure_scorer(self, peak, annotation):
        if hasattr(annotation,'predicted_fpt'):
            self.reference_fpt=annotation.predicted_fpt;
            if hasattr(annotation,'fpt_mask'):
                self.fpt_mask=annotation.fpt_mask;
                self.search_fpt=np.bitwise_and(self.reference_fpt,self.fpt_mask);
                return 0;
            else:
                return -1;
        else:
            return -1;
        
        
    def process_molecular_candidate_record(self, source_data_base, record):
        testvector=record['FPT'];
        testvector=np.bitwise_and(testvector,self.fpt_mask);
        v1=np.bitwise_and(testvector,self.search_fpt);
        v2=np.bitwise_or(testvector,self.search_fpt);
        c1=np.sum(np.unpackbits(v1));
        c2=np.sum(np.unpackbits(v2));
        if c2>0:
            record['Scores']['FPT']=float(c1)/c2;
        else:
            record['Scores']['FPT']=1.0;
            
for adduct in FingerPrintScorer.supported_adducts:
    global_supported_adducts.append(adduct);

registered_scorers[FingerPrintScorer.__name__]=FingerPrintScorer;

if __name__ == '__main__':
    if not test_as_single_process:
        multiprocessing.freeze_support();
    