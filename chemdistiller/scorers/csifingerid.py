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
from chemdistiller.utils.base64 import decode_from_base64, encode_to_base64;
from chemdistiller.filters.formula import FormulasFilter;

if not test_as_single_process:
    import multiprocessing;


class CSIFingerIDScorer:
    supported_adducts=SVM.supported_adducts;
    supported_isotopes=SVM.supported_isotopes;
    required_fields=['FPT', 'Formula', 'FormulaVector'];
    name='CSIFingerIDScorer';
    
    def __init__(self, process_id, settings):
        #self.supported_adducts=FingerPrintScorer.supported_adducts;
        #self.required_fields=FingerPrintScorer.required_fields;
        #self.supported_isotopes=FingerPrintScorer.supported_isotopes;
        self.settings=settings;
        self.process_id=process_id;

        
    def close(self):
        return 0;


    def configure_scorer(self, peak, annotation):
        peak = annotation.parent_peak;
        
        if 'csifingerid_count' in peak.parameters:
            self.csicount = int(peak.parameters['csifingerid_count']);
            self.csi_fpts = [];
            self.csi_formulas = [];
            self.csi_formulascores = [];
            self.csi_fpt_masks = [];
            for i in range(1, self.csicount + 1):
                fscore = float(peak.parameters['csifingerid_score_%s'%i]);
                if fscore > 0.0:
                    pred_fpt = peak.parameters['csifingerid_predfpt_%s'%i];
                    formula = peak.parameters['csifingerid_formula_%s'%i];
                    mask = peak.parameters['csifingerid_fptmask_%s'%i];
                    
                    #print(len(mask))
                    #print(mask)
                    mask = np.unpackbits(decode_from_base64(mask));

                    fpt = pred_fpt.split(',');
                    for i in range(len(fpt)):
                        fpt[i] = float(fpt[i]);
                    fpt = np.array(fpt, dtype = np.float64);
                    
                    full_fpt = np.zeros(mask.shape, dtype = np.float64);
                    full_fpt[mask > 0] = fpt[:];
                    
                    self.csi_fpts.append(full_fpt);
                    self.csi_formulas.append(FormulasFilter(formula));
                    self.csi_fpt_masks.append(mask);
                    self.csi_formulascores.append(fscore);
            self.csicount = len(self.csi_fpts)
            return 0;
        else:
            '''if os.path.isfile('e:/Imperial/TestDB/HDF5_CSIFingerID_SmallDB_new/errorlist.txt'):
                with open('e:/Imperial/TestDB/HDF5_CSIFingerID_SmallDB_new/errorlist.txt', 'a') as fout:
                    fout.write('%s,%s\n'%(peak.parent_spectrum.parameters['inchi'],peak.parent_spectrum.parameters['mode']));
            else:
                with open('e:/Imperial/TestDB/HDF5_CSIFingerID_SmallDB_new/errorlist.txt', 'w') as fout:
                    fout.write('%s,%s\n'%(peak.parent_spectrum.parameters['inchi'],peak.parent_spectrum.parameters['mode']));
            '''    
            
            return -1;
        
        
    def process_molecular_candidate_record(self, source_data_base, record):
        
        score = 0.0;
        
        for i in range(self.csicount):
            if not self.csi_formulas[i].rejected(record):
                
                testvector=record['FPT'];
                testvector = np.unpackbits(testvector)
                
                testvector = np.multiply(testvector, self.csi_fpt_masks[i]);
                
                v = np.vstack((testvector, self.csi_fpts[i]));
                
                c1=np.sum(np.min(v, axis = 0));
                c2=np.sum(np.max(v, axis = 0));

                if c2>0:
                    score = c1/c2;
                else:
                    score = 1.0;

                if self.settings[0]:
                    score = score * self.csi_formulascores[i];
            
        record['Scores']['CSI_FingerID'] = score;
            
for adduct in CSIFingerIDScorer.supported_adducts:
    global_supported_adducts.append(adduct);

#registered_scorers[CSIFingerIDScorer.__name__] = CSIFingerIDScorer;

if __name__ == '__main__':
    if not test_as_single_process:
        multiprocessing.freeze_support();
    