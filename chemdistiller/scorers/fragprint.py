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


#import numpy as np;
from operator import itemgetter;

#from chemdistiller.settings import test_as_single_process;
from chemdistiller.utils.periodictable import proton_mass, electron_mass;
from chemdistiller.msspectra.adducts import global_supported_adducts;
import multiprocessing;
from chemdistiller.settings import test_as_single_process, registered_scorers, registered_filters;


class FragPrintScorer:
    supported_adducts=set(['[M+H]+','[M-H]-']);
    supported_isotopes=set([0]);
    required_fields=['Frag'];
    delta_H=[0,1,-1];
    delta_mz=0.01;
    delta_mz_type=1;
    use_spectral_precision=True;
    pre_treat_type=1;
    name='FragPrintScorer';
    #delta_mz_type=0 - for ppm, 1 - for mass units.
    
    def __init__(self, process_id, settings):
        global test_as_single_process;
        
        self.supported_adducts=FragPrintScorer.supported_adducts;
        self.required_fields=FragPrintScorer.required_fields;
        self.supported_isotopes=FragPrintScorer.supported_isotopes;
        if settings:
            self.delta_H, self.delta_mz, self.delta_mz_type, self.pre_treat_type = settings;
        self.process_id=process_id;


    @staticmethod
    def condense_peak_list(peaks, delta_mz, delta_mz_type):
        counts=[1]*len(peaks);
        condensed=False;
                    
        if delta_mz_type==1:
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
                    
        elif delta_mz_type==0:
            for i in reversed(range(1,len(peaks))):
                item1=peaks[i];
                item2=peaks[i-1];
                if abs(item1[0]/counts[i]-item2[0])<=delta_mz*item2[0]/1e6:
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
    def __generate_subpeak_list(peak, overwrite=False):
        if overwrite==False:
            if hasattr(peak,'subpeak_list'):
                return;
        
        subcount=0;
        values=[];
        for subspectrum in peak.ms_spectra:
            if subspectrum.parameters['level']==2:
                 subspectrum.normalize_to_one();
                 subcount+=1;
                 for subpeak in subspectrum.peaks:
                     values.append([subpeak.mz, subpeak.intensity]);
                     
    
        if subcount>0:
            for value in values:
                value[1]=value[1]/subcount;
                 
        peak.subpeak_list=values;
        



    @staticmethod
    def tidy_up_peaks(peaks, required_fields):
        for peak in peaks:
            if hasattr(peak,'subpeak_list'):
                del peak.subpeak_list;
        if not ('Frag' in required_fields):
            for peak in peaks:
                for annotation in peak.annotations:
                    if not (annotation.mol_candidates is None):
                        for candidate in annotation.mol_candidates.mol_list:
                            try:
                                del candidate['Frag'];
                            except:
                                print('No "Frag" found');
        
    @staticmethod
    def generate_fragment_lists_for_peakblock(peaks, progress_report_queue, overwrite, blockindex):
        peakcount=len(peaks);
        for i in range(peakcount):
            FragPrintScorer.__generate_subpeak_list(peaks[i], overwrite);
            if not (progress_report_queue is None):
                    progress_report_queue.put((blockindex,i+1,peakcount));
                        
        
    @staticmethod    
    def generate_fragment_lists(peaks, nproc=1, overwrite=False):
        global test_as_single_process;
        for peak in peaks:
            FragPrintScorer.__generate_subpeak_list(peak, overwrite);

        
    def close(self):
        return 0;

    def configure_scorer(self, peak, annotation): 
        hydrogen_mass=proton_mass+electron_mass;
        if self.use_spectral_precision:
            if hasattr(peak,'parameters'):
                if 'subspectra_precision' in peak.parameters:
                    self.delta_mz=float(peak.parameters['subspectra_precision']);
                if 'subspectra_precision_type' in peak.parameters:
                    self.delta_mz_type=int(peak.parameters['subspectra_precision_type']);
                
        self.frag_test_vector=[];
        
        #self.correct_short_inchi=peak.parent_spectrum.parameters['shortinchi'];
        self.exactmass=peak.parent_spectrum.parameters['exactmass'];
        self.parent_ion_mz=peak.mz;
        delta_H_count=len(self.delta_H);
        if (not self.supported_adducts) or (annotation.adduct.definition in self.supported_adducts):
            if (not self.supported_isotopes) or (annotation.isotope in self.supported_isotopes):
                #prepare for querying
                self.adduct_type=annotation.adduct.definition;
                self.isotope=annotation.isotope;
                if annotation.isotope==0:
                    if annotation.adduct.definition=='[M+H]+': 
                        if self.delta_mz_type==1:
                            for dH in self.delta_H:
                                for subpeakindex in range(0,len(peak.subpeak_list)):
                                    subpeak=peak.subpeak_list[subpeakindex];
                                    if self.pre_treat_type==0:
                                        if abs(subpeak[0]-self.parent_ion_mz)<=self.delta_mz:
                                            self.frag_test_vector.append([annotation.adduct.get_mass(subpeak[0])+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                        else:    
                                            self.frag_test_vector.append([subpeak[0]+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                    elif self.pre_treat_type==1:
                                        self.frag_test_vector.append([annotation.adduct.get_mass(subpeak[0])+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                    else:
                                        self.frag_test_vector.append([subpeak[0]+dH*hydrogen_mass,subpeak[1]/delta_H_count]);

                        elif self.delta_mz_type==0:
                            for dH in self.delta_H:
                                for subpeakindex in range(0,len(peak.subpeak_list)):
                                    subpeak=peak.subpeak_list[subpeakindex];
                                    if self.pre_treat_type==0:
                                        if abs(subpeak[0]-self.parent_ion_mz)<=self.delta_mz*self.parent_ion_mz/1e6:
                                            self.frag_test_vector.append([annotation.adduct.get_mass(subpeak[0])+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                        else:    
                                            self.frag_test_vector.append([subpeak[0]+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                    elif self.pre_treat_type==1:
                                        self.frag_test_vector.append([annotation.adduct.get_mass(subpeak[0])+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                    else:
                                        self.frag_test_vector.append([subpeak[0]+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                
                    elif annotation.adduct.definition=='[M-H]-':

                        if self.delta_mz_type==1:
                            for dH in self.delta_H:
                                for subpeakindex in range(0,len(peak.subpeak_list)):
                                    subpeak=peak.subpeak_list[subpeakindex];
                                    if self.pre_treat_type==0:
                                        if abs(subpeak[0]-self.parent_ion_mz)<=self.delta_mz:
                                            self.frag_test_vector.append([annotation.adduct.get_mass(subpeak[0])+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                        else:    
                                            self.frag_test_vector.append([subpeak[0]+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                    elif self.pre_treat_type==1:
                                        self.frag_test_vector.append([annotation.adduct.get_mass(subpeak[0])+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                    else:
                                        self.frag_test_vector.append([subpeak[0]+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                    
                        elif self.delta_mz_type==0:
                            for dH in self.delta_H:
                                for subpeakindex in range(0,len(peak.subpeak_list)):
                                    subpeak=peak.subpeak_list[subpeakindex];
                                    if self.pre_treat_type==0:                                    
                                        if abs(subpeak[0]-self.parent_ion_mz)<=self.delta_mz*self.parent_ion_mz/1e6:
                                            self.frag_test_vector.append([annotation.adduct.get_mass(subpeak[0])+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                        else:    
                                            self.frag_test_vector.append([subpeak[0]+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                    elif self.pre_treat_type==1:
                                        self.frag_test_vector.append([annotation.adduct.get_mass(subpeak[0])+dH*hydrogen_mass,subpeak[1]/delta_H_count]);
                                    else:
                                        self.frag_test_vector.append([subpeak[0]+dH*hydrogen_mass,subpeak[1]/delta_H_count]);

                        
                    elif annotation.adduct.definition=='[2M+H]+':
                        print('TODO');
                        #Do not forget charge if charged adduct!
                    
                    
                self.frag_test_vector=sorted(self.frag_test_vector, key=itemgetter(0), reverse=True);
                self.condense_peak_list(self.frag_test_vector, self.delta_mz, self.delta_mz_type);
                
                return 0;

            else:
                return -2; #Unsupported isotope;
        else:
            return -1; #Unsupported adduct;
        
        
        #return 0 if accepted, error code otherwise
        
    def process_molecular_candidate_record(self, source_data_base, record):
        global test_as_single_process;
        global debug_log;
        testvector=record['Frag'];

        #if source_data_base.charged:
        #    if test_as_single_process:
        #        print('Charged');
        #        chargesvector=record['FragCharge'];
        #    #Expand charged fragments here

        counts=[1]*len(testvector);
        
        condensed=False;
        
        if self.delta_mz_type==1:
                        for i in reversed(range(1,len(testvector))):
                            item1=testvector[i]/counts[i];
                            item2=testvector[i-1];
                            if abs(item1-item2)<=self.delta_mz:
                                testvector[i-1]+=item1;
                                counts[i-1]+=counts[i];
                                del testvector[i];
                                del counts[i];
                                condensed=True;
                
        elif self.delta_mz_type==0:

                        for i in reversed(range(1,len(testvector))):
                            item1=testvector[i]/counts[i];
                            item2=testvector[i-1];
                            if abs(item1-item2)<=self.delta_mz*item1/1e6:                            
                                testvector[i-1]+=item1;
                                counts[i-1]+=counts[i];
                                del testvector[i];
                                del counts[i];
                                condensed=True;
                
        if condensed:
            for i in range(len(testvector)):
                    testvector[i]=testvector[i]/counts[i];
        
            
        score=0.0;
        # Calculate score multiplication here
        i=0;
        j=0;
        max_i=len(testvector);
        max_j=len(self.frag_test_vector);
        
        
        if self.delta_mz_type==1:
            #absolute delta MZ mass units
            while i<max_i and j<max_j:
                if abs(testvector[i]-self.frag_test_vector[j][0])<=self.delta_mz: 
                    score+=self.frag_test_vector[j][1];
                    i+=1;
        
                elif testvector[i]>self.frag_test_vector[j][0]:
                    i+=1;
                else:                    
                    j+=1;
        elif self.delta_mz_type==0:
            #ppm delta MZ
            while i<max_i and j<max_j:
                if abs(testvector[i]-self.frag_test_vector[j][0])<=self.delta_mz*self.frag_test_vector[j][0]/1e6:
                    score+=self.frag_test_vector[j][1];
                    i+=1;
                elif testvector[i]>self.frag_test_vector[j][0]:
                    i+=1;
                else:                    
                    j+=1;
        
        

        record['Scores']['Frag']=score;
        
            
for adduct in FragPrintScorer.supported_adducts:
    global_supported_adducts.append(adduct);
            
registered_scorers[FragPrintScorer.__name__]=FragPrintScorer;
        
if __name__=='__main__':

    multiprocessing.freeze_support();


    

    