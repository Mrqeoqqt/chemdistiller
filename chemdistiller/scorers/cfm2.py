# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 15:03:03 2016

@author: Dr. Ivan Laponogov
"""

import sys;
import os;
from math import log, exp, sqrt;
from copy import deepcopy;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);


#import numpy as np;
from operator import itemgetter;

#from chemdistiller.settings import test_as_single_process;
from chemdistiller.utils.periodictable import proton_mass, electron_mass, hydrogen_mass;
from chemdistiller.msspectra.adducts import global_supported_adducts;
import multiprocessing;
from chemdistiller.settings import test_as_single_process, registered_scorers, registered_filters;


class CFMScorer2:
    supported_adducts=set(['[M+H]+','[M-H]-']);
    supported_isotopes=set([0]);
    required_fields=['CFM2'];
    delta_mz=0.01;
    delta_mz_type=1;
    score_types=['CFM2_DotProduct'];# CFM2_DotProduct, CFM2_Jaccard, CFM2_WPrecision, CFM2_WRecall, CFM2_Precision, CFM2_Recall
    treat_spectra=0;# 0 - merge, 1 - sqrt(e0^2+e1^2+e2^2), 2 - e0*e1*e2
    use_spectral_precision=True;
    name='CFMScorer2';
    #delta_mz_type=0 - for ppm, 1 - for mass units.
    
    def __init__(self, process_id, settings):
        global test_as_single_process;
        
        self.supported_adducts=CFMScorer2.supported_adducts;
        self.required_fields=CFMScorer2.required_fields;
        self.supported_isotopes=CFMScorer2.supported_isotopes;
        if settings:
            self.delta_mz, self.delta_mz_type, self.score_types, self.treat_spectra = settings;
        self.score_types=set(self.score_types);
        self.process_id=process_id;
        self.empty_record_count = 0;
        self.non_empty_record_count = 0;

    
    @staticmethod
    def condense_peak_list(peaks, delta_mz):
        counts=[1]*len(peaks);
        condensed=False;
        
        for i in reversed(range(1,len(peaks))):
            item1=peaks[i];
            item2=peaks[i-1];
            if abs(item1[0]/counts[i]-item2[0])<=delta_mz:#for binning test                            
                item2[0]+=item1[0];
                item2[1]+=item1[1];
                counts[i-1]+=counts[i];
                del peaks[i];
                del counts[i];
                condensed=True;
                    
        if condensed:
            for i in range(len(peaks)):
                peaks[i][0]=peaks[i][0]/counts[i];
    
    
    @staticmethod
    def merge_peak_lists(peak_lists):
        if peak_lists:
            merged=[];
            for peak_list in peak_lists:
                merged+=deepcopy(peak_list);
            peak_lists[0]=sorted(merged, key=itemgetter(0));
            del peak_lists[1:];
    
    @staticmethod
    def normalize_peaks(data):
        sumintens=0.0;
        for i in range(len(data)):
            sumintens+=data[i][1];
        if sumintens>0.0:
            for i in range(len(data)):
                data[i][1]=data[i][1]/sumintens;


    @staticmethod
    def _mean_mz_weighted(peaks):
        sum_mz=0.0;
        sum_i=0.0;
        for peak in peaks:
            sum_mz+=peak[0]*peak[1];
            sum_i+=peak[1];
        if sum_i>0.0:
            return sum_mz/sum_i;
        return 0.0

    @staticmethod
    def _normalize_and_merge_spectra(peak_list, precision):
        result=[];
        for i in range(len(peak_list)):
            peaks=peak_list[i];
            
            for peak in peaks:
                result.append([peak[0], peak[1]/len(peak_list)]);
            result=sorted(result,key=itemgetter(0));
        CFMScorer2.condense_peak_list(result, precision);
        return result;
        
    
    @staticmethod
    def _np_spectrum_to_list(vector):
        result=[];
        for i in range(len(vector)):
            result.append([float(vector[i,0]), float(vector[i,1])]);
        return result;

    
    def _generate_CFM2_peak_list(self, peak, overwrite=False):
        if overwrite==False:
            if hasattr(peak,'CFM2_peak_list'):
                return;
        
        self.CFM2_peak_list=[];
        for subspectrum in peak.ms_spectra:
            if subspectrum.parameters['level']==2:
                 values=[];
                 subspectrum.normalize_to_one();
                 for subpeak in subspectrum.peaks:
                     values.append([subpeak.mz, subpeak.intensity]);
                 self.CFM2_peak_list.append(values);
        
        if self.delta_mz_type==0:
            CFMScorer2._CFM_peak_list_to_ppm(self.CFM2_peak_list);
                 
        if self.treat_spectra==0:
            CFMScorer2.merge_peak_lists(self.CFM2_peak_list);
            CFMScorer2.condense_peak_list(self.CFM2_peak_list[0], self.delta_mz);
            
            
        self.CFM2_mean_mz_weighted=[];    
        for i in self.CFM2_peak_list:
            CFMScorer2.normalize_peaks(i);
            self.CFM2_mean_mz_weighted.append(CFMScorer2._mean_mz_weighted(i));
        
        
        
                 
        
        
        
    @staticmethod
    def _CFM_peak_list_to_ppm(CFM_peak_list):
        for peak_list in CFM_peak_list:
            for peak in peak_list:
                peak[0]=log(peak[0]/1.00794)*1.0e6;

    @staticmethod
    def tidy_up_peaks(peaks, required_fields):
        for peak in peaks:
            if hasattr(peak,'CFM2_peak_list'):
                del peak.CFM2_peak_list;
            if hasattr(peak,'CFM2_mean_mz_weighted'):
                del peak.CFM2_mean_mz_weighted;
            
        if not ('CFM2' in required_fields):
            for peak in peaks:
                for annotation in peak.annotations:
                    if not (annotation.mol_candidates is None):
                        for candidate in annotation.mol_candidates.mol_list:
                            if     'CFM2_DotProduct' in  candidate:
                                del candidate['CFM2_DotProduct'];
                            elif   'CFM2_Jaccard' in  candidate:
                                del candidate['CFM2_Jaccard'];
                            elif   'CFM2_WPrecision' in  candidate:
                                del candidate['CFM2_WPrecision'];
                            elif   'CFM2_WRecall' in  candidate:
                                del candidate['CFM2_WRecall'];
                            elif   'CFM2_Precision' in  candidate:
                                del candidate['CFM2_Precision'];
                            elif   'CFM2_Recall' in  candidate:
                                del candidate['CFM2_Recall'];
                            elif   'CFM2' in  candidate:
                                del candidate['CFM2'];    
                            
    
        
    def close(self):
        return 0;
    
    @staticmethod    
    def _qmatch(mzref, mzcalc, tolerance):
        idx1=[];
        idx2=[];
        refcount=len(mzref);
        calccount=len(mzcalc);
        i1=0;
        i2=0;
        while (i1<refcount) and (i2<calccount):
            if abs(mzref[i1][0]-mzcalc[i2][0])<=tolerance:
                idx1.append(i1);
                idx2.append(i2);
                i1+=1;
                i2+=1;
            elif mzref[i1][0]<mzcalc[i2][0]:
                idx1.append(i1);
                idx2.append(-1);
                i1+=1;
            else:
                idx1.append(-1);
                idx2.append(i2);
                i2+=1;
                
        return idx1, idx2
        

    def configure_scorer(self, peak, annotation): 
        if self.use_spectral_precision:
            if hasattr(peak,'parameters'):
                if 'subspectra_precision' in peak.parameters:
                    self.delta_mz=float(peak.parameters['subspectra_precision']);
                if 'subspectra_precision_type' in peak.parameters:
                    self.delta_mz_type=int(peak.parameters['subspectra_precision_type']);

        self.exactmass=peak.parent_spectrum.parameters['exactmass'];
        self.parent_ion_mz=peak.mz;
        self.mode=peak.parent_spectrum.parameters['mode'];
        
        self._generate_CFM2_peak_list(peak, True);
        
        if (not self.supported_adducts) or (annotation.adduct.definition in self.supported_adducts):
            if (not self.supported_isotopes) or (annotation.isotope in self.supported_isotopes):
                self.adduct_type=annotation.adduct.definition;
                self.isotope=annotation.isotope;
                return 0;
            else:
                return -2; #Unsupported isotope;
        else:
            return -1; #Unsupported adduct;
        
        
        #return 0 if accepted, error code otherwise
        
    def process_molecular_candidate_record(self, source_data_base, record):
        global test_as_single_process;
        global debug_log;
        if self.mode==1:
            testvector=record['CFM2'][1];
        else:
            testvector=record['CFM2'][0];
        
        if testvector:
            testvector[0]=CFMScorer2._np_spectrum_to_list(testvector[0]);
            testvector[1]=CFMScorer2._np_spectrum_to_list(testvector[1]);
            testvector[2]=CFMScorer2._np_spectrum_to_list(testvector[2]);
        
            for i in testvector:
                CFMScorer2.normalize_peaks(i);
                
            if self.delta_mz_type==0:
                CFMScorer2._CFM_peak_list_to_ppm(testvector);
            self.non_empty_record_count += 1;                
        else:
            self.empty_record_count += 1;        
        
        if len(testvector)==0:
            CFM_DotProduct=[0.0];
            CFM_Jaccard=[0.0];
            CFM_WPrecision=[0.0];
            CFM_WRecall=[0.0];
            CFM_Precision=[0.0];
            CFM_Recall=[0.0];
        else:
            CFM_DotProduct=[];
            CFM_Jaccard=[];
            CFM_WPrecision=[];
            CFM_WRecall=[];
            CFM_Precision=[];
            CFM_Recall=[];
            
            
            
            mz0=CFMScorer2._mean_mz_weighted(testvector[0]);
            mz1=CFMScorer2._mean_mz_weighted(testvector[1]);
            mz2=CFMScorer2._mean_mz_weighted(testvector[2]);
            
            for i in range(len(self.CFM2_peak_list)):
                mean_mz=self.CFM2_mean_mz_weighted[i];
                refdata=self.CFM2_peak_list[i];
                
                bestindex=0;
                bestmz=abs(mean_mz-mz0);
                if abs(mean_mz-mz1)<bestmz:
                    bestmz=abs(mean_mz-mz1);
                    bestindex=1;
                if abs(mean_mz-mz2)<bestmz:
                    #bestmz=abs(mean_mz-mz2);
                    bestindex=2;
                calcdata=testvector[bestindex];    
            
                match=CFMScorer2._qmatch(refdata, calcdata, self.delta_mz);
                
                a=0;
                b=0;
                dotproduct=0.0;
                wrecall=0.0;
                wprecision=0.0;
                recall=0;
                precision=0;
                for i in range(len(match[0])):    
                    if match[0][i]>-1 and match[1][i]>-1:
                        a+=1;
                        dotproduct+=refdata[match[0][i]][1]*calcdata[match[1][i]][1];
                        wrecall+=refdata[match[0][i]][1];
                        recall+=1;
                        wprecision+=calcdata[match[1][i]][1];
                        precision+=1;
                        
                    if match[0][i]>-1 or match[1][i]>-1:
                        b+=1;
            
                if b!=0:
                    jaccard=float(a)/b;
                else:
                    jaccard=1.0;
                    
                if len(refdata)>0:
                    recall=float(recall)/len(refdata);
            
                if len(calcdata)>0:
                    precision=float(precision)/len(calcdata);
    
                CFM_DotProduct.append(dotproduct);
                CFM_Jaccard.append(jaccard);
                CFM_WPrecision.append(wprecision);
                CFM_WRecall.append(wrecall);
                CFM_Precision.append(precision);
                CFM_Recall.append(recall);

            
            if self.treat_spectra==1:
                    #sqrt
                    CFM_DotProduct[0]=CFM_DotProduct[0]*CFM_DotProduct[0];
                    CFM_Jaccard[0]=CFM_Jaccard[0]*CFM_Jaccard[0];
                    CFM_WPrecision[0]=CFM_WPrecision[0]*CFM_WPrecision[0];
                    CFM_WRecall[0]=CFM_WRecall[0]*CFM_WRecall[0];
                    CFM_Precision[0]=CFM_Precision[0]*CFM_Precision[0];
                    CFM_Recall[0]=CFM_Recall[0]*CFM_Recall[0];
                    for i in range(1, len(CFM_DotProduct)):
                        CFM_DotProduct[0]+=CFM_DotProduct[i]*CFM_DotProduct[i];
                        CFM_Jaccard[0]+=CFM_Jaccard[i]*CFM_Jaccard[i];
                        CFM_WPrecision[0]+=CFM_WPrecision[i]*CFM_WPrecision[i];
                        CFM_WRecall[0]+=CFM_WRecall[i]*CFM_WRecall[i];
                        CFM_Precision[0]+=CFM_Precision[i]*CFM_Precision[i];
                        CFM_Recall[0]+=CFM_Recall[i]*CFM_Recall[i];
                    CFM_DotProduct[0]=sqrt(CFM_DotProduct[0]); 
                    CFM_Jaccard[0]=sqrt(CFM_Jaccard[0]); 
                    CFM_WPrecision[0]=sqrt(CFM_WPrecision[0]); 
                    CFM_WRecall[0]=sqrt(CFM_WRecall[0]); 
                    CFM_Precision[0]=sqrt(CFM_Precision[0]); 
                    CFM_Recall[0]=sqrt(CFM_Recall[0]); 
                
            elif self.treat_spectra==2:
                    #mult
                    for i in range(1, len(CFM_DotProduct)):
                        CFM_DotProduct[0]*=CFM_DotProduct[i];
                        CFM_Jaccard[0]*=CFM_Jaccard[i];
                        CFM_WPrecision[0]*=CFM_WPrecision[i];
                        CFM_WRecall[0]*=CFM_WRecall[i];
                        CFM_Precision[0]*=CFM_Precision[i];
                        CFM_Recall[0]*=CFM_Recall[i];
            
        if 'CFM2_DotProduct' in self.score_types:
            record['Scores']['CFM2_DotProduct']=CFM_DotProduct[0];
            
        if 'CFM2_Jaccard' in self.score_types:
            record['Scores']['CFM2_Jaccard']=CFM_Jaccard[0];
            
        if 'CFM2_WPrecision' in self.score_types:
            record['Scores']['CFM2_WPrecision']=CFM_WPrecision[0];
            
        if 'CFM2_WRecall' in self.score_types:
            record['Scores']['CFM2_WRecall']=CFM_WRecall[0];
            
        if 'CFM2_Precision' in self.score_types:
            record['Scores']['CFM2_Precision']=CFM_Precision[0];
            
        if 'CFM2_Recall' in self.score_types:
            record['Scores']['CFM2_Recall']=CFM_Recall[0];                      
                

            

        
        
            
for adduct in CFMScorer2.supported_adducts:
    global_supported_adducts.append(adduct);
            
#registered_scorers[CFMScorer2.__name__]=CFMScorer2;
        
if __name__=='__main__':

    multiprocessing.freeze_support();


    

    