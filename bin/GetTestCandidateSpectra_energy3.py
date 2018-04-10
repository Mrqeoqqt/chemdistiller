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

from chemdistiller.settings import test_as_single_process;
from chemdistiller.msspectra.manager import SpectralManager;
from chemdistiller.msspectra.adducts import get_adduct_by_name;
from chemdistiller.annotations.peakannotator import PeakAnnotator;
from chemdistiller.scorers.totalscore import total_multiplicative_score;
from chemdistiller.scorers.fragprint import FragPrintScorer;
from chemdistiller.filters.formula import FormulasFilter;

from chemdistiller.io.jsonio import spectra_to_json;

from chemdistiller.annotations.peakannotation import MSPeakAnnotation, merge_annotations;

from operator import itemgetter;

import multiprocessing;



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
    
def process_output(peaklist, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir);
    
    foutst=open('%s.txt'%outdir,'w');                          
    foutst.write('%s\n'%len(peaklist));
    

    for i in range(len(peaklist)):
        foutst.write('%s\tInChI=1S/%s\t0\n'%(i,peaklist[i].parent_spectrum.parameters['inchi']));
        peak=peaklist[i];
        fout=open('%s/%s.txt'%(outdir,i),'w');                          
        fout.write('energy0\n');
        for j in peak.subpeak_list_low:
            fout.write('%.5f %.5f\n'%(j[0],j[1]));
        fout.write('energy1\n');
        for j in peak.subpeak_list_med:
            fout.write('%.5f %.5f\n'%(j[0],j[1]));    
        fout.write('energy2\n');
        for j in peak.subpeak_list_high:
            fout.write('%.5f %.5f\n'%(j[0],j[1]));    
        fout.close();
        
    foutst.close();

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
    
def make_three_energies(peak, precision):
    if len(peak.ms_spectra)==1:
        peak.subpeak_list_low=normalize_and_merge_spectra(peak.ms_spectra, precision);
        peak.subpeak_list_med=normalize_and_merge_spectra(peak.ms_spectra, precision);
        peak.subpeak_list_high=normalize_and_merge_spectra(peak.ms_spectra, precision);
    elif len(peak.ms_spectra)==2:
        peak.subpeak_list_low=normalize_and_merge_spectra(peak.ms_spectra[0:1], precision);
        peak.subpeak_list_med=normalize_and_merge_spectra(peak.ms_spectra, precision);
        peak.subpeak_list_high=normalize_and_merge_spectra(peak.ms_spectra[-1:], precision);
    else:
        group1=[];
        group2=[];
        group3=[];
        mz1=peak.ms_spectra[0].mean_mz_weighted;
        mz3=peak.ms_spectra[-1].mean_mz_weighted;
        mz_low=(mz3-mz1)/3+mz1;
        mz_high=(mz3-mz1)*2/3+mz1;
        for spectrum in peak.ms_spectra:
            if spectrum.mean_mz_weighted<=mz_low:
                group1.append(spectrum);
            elif spectrum.mean_mz_weighted>=mz_high:
                group3.append(spectrum);                
            else:
                group2.append(spectrum);
        peak.subpeak_list_low=normalize_and_merge_spectra(group3, precision);
        if not group2:
            group2=group1+group3;
        peak.subpeak_list_med=normalize_and_merge_spectra(group2, precision);
        peak.subpeak_list_high=normalize_and_merge_spectra(group1, precision);                
    

if __name__=='__main__':

    
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
        for j in reversed(range(len(peak.ms_spectra))):
            if not 'collision_energy' in peak.ms_spectra[j].parameters:
                del peak.ms_spectra[j];
            else:
                get_mean_mz_weighted(peak.ms_spectra[j]);
        if len(peak.ms_spectra)<1:
            del peaks[i];
        else:
            peak.ms_spectra=sorted(peak.ms_spectra, key=lambda spec: spec.mean_mz_weighted);
            make_three_energies(peak, 0.01);
        
    peaks=sorted(peaks,key=lambda peak: peak.mz);
    
    neg=[];
    pos=[];
    for peak in peaks:
        if peak.parent_spectrum.parameters['mode']==-1:
            neg.append(peak);
        else:
            pos.append(peak);
    
    print('Positive mode: %s'%len(pos));
    print('Negative mode: %s'%len(neg));    
    
    #neg=neg[:1150];
    #pos=pos[:2500];
    #neg=neg[1150:];
    #pos=pos[2500:];
    
    #print('Positive mode: %s'%len(pos));
    #print('Negative mode: %s'%len(neg));    
    '''
    fout=open('energy_mzw_correlation_neg.txt','w');
    for peak in neg:
        for spectrum in peak.ms_spectra:
            if spectrum.parameters['collision_energy']>0.0:
                fout.write('%s\t%s\t%s\n'%(peak.parent_spectrum.parameters['inchi'], spectrum.parameters['collision_energy'], spectrum.mean_mz_weighted));
    fout.close();

    fout=open('energy_mzw_correlation_pos.txt','w');
    for peak in pos:
        for spectrum in peak.ms_spectra:
            if spectrum.parameters['collision_energy']>0.0:
                fout.write('%s\t%s\t%s\n'%(peak.parent_spectrum.parameters['inchi'], spectrum.parameters['collision_energy'], spectrum.mean_mz_weighted));
    fout.close();
    '''
    #for i in range(10):
    #    count=int(len(neg)/10.0*(i+1));
    #    process_output(neg[:count],'g:/CFM/learning_courve/%s/negative/data'%(i+1));
    #    
    #    count=int(len(pos)/10.0*(i+1));
    #    process_output(pos[:count],'g:/CFM/learning_courve/%s/positive/data'%(i+1));
    
    process_output(neg,'g:/CFM/control/negative_0');
    process_output(pos,'g:/CFM/control/positive_0');
    
                 
    specmanager.close();
    
        