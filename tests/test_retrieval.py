# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 15:51:17 2016

@author: Dr. Ivan Laponogov
"""
import sys;
if __name__=='__main__':

    
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

import os;
#dir_path = os.path.dirname(os.path.realpath(__file__));
chemdistiller_path = os.path.abspath(os.path.join(os.getcwd(), ".."));
#print('Expected chemdistiller path: %s'%chemdistiller_path);
sys.path.append(chemdistiller_path);

import time;

from chemdistiller.settings import test_as_single_process;
from chemdistiller.msspectra.manager import SpectralManager;
from chemdistiller.msspectra.adducts import get_adduct_by_name;
from chemdistiller.annotations.peakannotator import PeakAnnotator;
from chemdistiller.scorers.totalscore import total_multiplicative_score;
from chemdistiller.scorers.fragprint import FragPrintScorer;
from chemdistiller.filters.formula import FormulasFilter;
from chemdistiller.filters.inchi import InChiFilter;

from chemdistiller.annotations.peakannotation import MSPeakAnnotation, merge_annotations;

import multiprocessing;

if __name__=='__main__':
    #test_spectral_database_path='e:/Imperial/TestDB/smalltest';
    #test_spectral_database_path='e:/Imperial/TestDB/Selected80';
    
    
    def logit(message):
        global logfile;
        print(message);
        logfile.write('%s\n'%message);
    
    def logtime():
        global starttime;
        global logfile;
        dtime=time.time()-starttime;
        logit('+ %s:%s:%s'%(int(dtime)//3600,
                        int(dtime)%3600//60,
                        int(dtime)%60));
    
    
    def load_inchi_list(fname):
        finp=open(fname,'r');
        result=[];
        for s in finp:
            s=s.rstrip('\n').split('\t');
            inchi=s.pop();
            inchi=inchi.replace('1S/','');
            result.append(inchi);
        finp.close();
        return(result);
    
    
    def remove_spectra_not_in_inchi_list(selected_InChi,dblist):
        logit('Removing spectra with InChi-s not in the list');
        
        '''
        for i in reversed(range(len(specmanager.ms_spectra))):
            spectrum=specmanager.ms_spectra[i];
            
            inchi=spectrum.parameters['inchi'];
            if not (inchi in selected_InChi):
                del specmanager.ms_spectra[i];
                
        '''
        
        for i in reversed(range(len(specmanager.ms_spectra))):
            spectrum=specmanager.ms_spectra[i];
            
            #mass=float(spectrum.parameters['exactmass']);
            inchi=spectrum.parameters['inchi'];
            #mode=spectrum.parameters['mode'];
            
            keepspectrum=False;
            if (inchi in selected_InChi):
                for peakindex in reversed(range(len(spectrum.peaks))):
                                    
                    peak=spectrum.peaks[peakindex];
                    for subspectrumindex in reversed(range(len(peak.ms_spectra))):
                        subspectrum=peak.ms_spectra[subspectrumindex];
                        if not (subspectrum.parameters['dbsource'] in dblist):
                            del peak.ms_spectra[subspectrumindex];
    
                    keep=len(peak.ms_spectra)>0;
                    if keep==False:
                         del spectrum.peaks[peakindex];
                    else:
                        keepspectrum=True;
                        
            if keepspectrum==False:
                del specmanager.ms_spectra[i];
                
        
                        
    
                        
        logit('Finshed removing spectra. Total number left: %s'%len(specmanager.ms_spectra));    
    
    
    
    def remove_low_precision_spectra(precision):
        logit('Removing low precision spectra');
        
        
        for i in reversed(range(len(specmanager.ms_spectra))):
            spectrum=specmanager.ms_spectra[i];
            
            mass=float(spectrum.parameters['exactmass']);
            mode=spectrum.parameters['mode'];
            keepspectrum=False;
            
            for peakindex in reversed(range(len(spectrum.peaks))):
                keep=False;
                peak=spectrum.peaks[peakindex];
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
                
        
                        
        logit('Finshed removing spectra. Total number left: %s'%len(specmanager.ms_spectra));    
    
    
    class ResBlock:
    
        def __init__(self):
            global max_results_per_query;
            self.count=1;
            self.correct=0;
            self.best_results=[0]*max_results_per_query;
            self.worst_results=[0]*max_results_per_query;
            self.correct_top_mean=0;
            self.correct_top_max=0;
        
        def percent(self):
            if self.count>0:
                self.correct=self.correct*100/self.count;
                for i in range(len(self.best_results)):
                    self.best_results[i]=self.best_results[i]*100/self.count;
                for i in range(len(self.worst_results)):
                    self.worst_results[i]=self.worst_results[i]*100/self.count;
                self.correct_top_max=self.correct_top_max*100/self.count;
                self.correct_top_mean=self.correct_top_mean*100/self.count;
    
    
    
    
    
    def correct_precursor_mz(spectral_manager):
        
        for spectrum in specmanager.ms_spectra:
            mass=float(spectrum.parameters['exactmass']);
            mode=spectrum.parameters['mode'];
            for peak in spectrum.peaks:
                if 'ion_type' in peak.parameters:            
                    adduct=get_adduct_by_name(peak.parameters['ion_type'],mode);
                    if not (adduct is None):
                        peak.mz=adduct.get_mz(mass);
    
    
    
    
    def testrun(scorers_list=[], scorers_settings=[], batchindex=-1):
    
            #for batchindex in range(5):
            
            logit('Assessing Fragprint based search. Batch %s...'%batchindex);
            
            total_peaks=0;
    
            best_results_neg=[0]*max_results_per_query;
            worst_results_neg=[0]*max_results_per_query;
            return_number_neg=[0]*max_results_per_query;
            neg_mode_count=0;
            
            best_results_pos=[0]*max_results_per_query;
            worst_results_pos=[0]*max_results_per_query;
            return_number_pos=[0]*max_results_per_query;
            pos_mode_count=0;
            
                
            logtime();
        
            logit('Assuming correct adduct info and 0 isotope.');
            logit('No formula assumption. Assuming Mass precision: %s ppm'%ppm);
            logit('Only considering [M+H]+ and [M-H]- adducts for now');
            
            peaks=[];
            for spectrum_index in range(len(specmanager.ms_spectra)):
                print('%s of %s spectra preprocessed'%(spectrum_index,len(specmanager.ms_spectra)));
                spectrum=specmanager.ms_spectra[spectrum_index];
                if batchindex==-1 or int(spectrum.parameters['crossvalidation_batch_index'])==batchindex:
                    for peak in spectrum.peaks:
                        if 'ion_type' in peak.parameters:
                            if peak.parameters['ion_type'] in FragPrintScorer.supported_adducts:
                                peak.ppm=ppm; #Assuming 0 isotope states only.
                                shortinchi=peak.parent_spectrum.parameters['shortinchi'];                            
                                #formula=peak.parent_spectrum.parameters['formula'];                            
                                #correct_inchi=InChiFilter(shortinchi,True,4);
                                #correct_formula=FormulasFilter(formula);
                                
                                annotation=MSPeakAnnotation(adduct=get_adduct_by_name(peak.parameters['ion_type'],\
                                    spectrum.parameters['mode']), isotope=0, formula_scorer=None,\
                                    filters=[], scores={'AdductIsotopeScore':1.0});
                                if not (annotation.adduct is None):
                                    peak.annotations=[annotation];
                                    annotation.parent_peak=peak;
                                    peaks.append(peak);
                                else:
                                    logit('Unsupported Adduct: %s'%peak.parameters['ion_type']);
                        
            logit('Peaks to annotate: %s'%len(peaks));
    
            #
            peak_annotator.annotate_peaks(peaks, test_chemical_databases, \
                scorers_list=scorers_list, scorers_settings=scorers_settings, total_score=total_multiplicative_score,\
                results_limit=max_results_per_query, save_memory=False, batch_index=batchindex, ppm=ppm, overwrite=True);
            
            for i in reversed(range(len(peaks))):
                if not peaks[i].annotations[0].mol_candidates:
                    del peaks[i];
            
            logit('Finished annotating... Annotated peaks :%s'%len(peaks));    
    
    
            logtime();
            
            total_peaks=len(peaks);
            
            for peakindex in range(total_peaks):
                peak=peaks[peakindex];
                peakmode=0;
                if peak.parent_spectrum.parameters['mode']==1:
                    pos_mode_count+=1;
                    peakmode=1;
                elif peak.parent_spectrum.parameters['mode']==-1:
                    neg_mode_count+=1;
                    peakmode=-1;
                    
                print('%s of %s '%(peakindex,total_peaks));
                shortinchi=peak.parent_spectrum.parameters['shortinchi'];
    
                for annotation in peak.annotations:
                    annotation.min_correct=-1;
                    annotation.max_correct=-1;
                    
                    #Get annotation score
                    
                    annotation.mean_score=0.0;
                    annotation.max_score=0.0;
                    for index in range(len(annotation.mol_candidates.mol_list)):
                        total_score=annotation.mol_candidates.mol_list[index]['TotalScore'];
                        #print(annotation.mol_candidates.mol_list[index]);
                        annotation.mean_score+=total_score;
                        if total_score>annotation.max_score:
                            annotation.max_score=total_score;
                    if len(annotation.mol_candidates.mol_list)>0:
                        annotation.mean_score=annotation.mean_score/len(annotation.mol_candidates.mol_list);
                        if len(annotation.mol_candidates.mol_list)<=max_results_per_query:
                            if peakmode==1:
                                return_number_pos[len(annotation.mol_candidates.mol_list)-1]+=1;
                            elif peakmode==-1:
                                return_number_neg[len(annotation.mol_candidates.mol_list)-1]+=1;                            
                        
                    
                    
                    
    
                    for index in reversed(range(1,len(annotation.mol_candidates.mol_list))):
                        if annotation.mol_candidates.mol_list[index-1]['ShortInChi']==annotation.mol_candidates.mol_list[index]['ShortInChi']:
                            del annotation.mol_candidates.mol_list[index]; # For the purpose of statistics condensing sequential identical Short InChi-s
                            
                    #print(len(annotation.mol_candidates.mol_list));
                    #print('Correct: %s'%shortinchi);
                    
                    for index in range(len(annotation.mol_candidates.mol_list)):
                        #print(annotation.mol_candidates.mol_list[index]['ShortInChi']);
                        if annotation.mol_candidates.mol_list[index]['ShortInChi']==shortinchi:
                            #print('Got it!');
                            if annotation.min_correct==-1:
                                annotation.min_correct=index;
                            annotation.max_correct=index;
    
                    if annotation.min_correct>-1:
                        if peakmode==1:
                            for i in range(len(best_results_pos)):
                                if annotation.min_correct<=i:
                                    best_results_pos[i]+=1;
                                if annotation.max_correct<=i:
                                    worst_results_pos[i]+=1;
                        elif peakmode==-1:
                            for i in range(len(best_results_neg)):
                                if annotation.min_correct<=i:
                                    best_results_neg[i]+=1;
                                if annotation.max_correct<=i:
                                    worst_results_neg[i]+=1;
                    #else:
                    #    logit('No candidate for %s (Mass: %s, peak mode: %s, peak_mz %s, deltaM %s)'%(shortinchi,peak.parent_spectrum.parameters['exactmass'],peakmode, peak.mz,abs(peak.mz-float(peak.parent_spectrum.parameters['exactmass']))));
            logit('Finished Calculating Retrieval Stats.');
            logit('Positive Mode Count: %s'%pos_mode_count);
            logit('Negative Mode Count: %s'%neg_mode_count);
            logit('Total Peaks: %s'%total_peaks);
            logtime();
    
            logit('Positive Mode (Total: %s ):'%pos_mode_count);
            if pos_mode_count>0:
                for i in range(max_results_per_query):
                    logit('Correct within first\t%s:\tBest:\t%s%%\tWorst:\t%s%%\tCandidateCount:\t%s'%(i+1,best_results_pos[i]*100/pos_mode_count,worst_results_pos[i]*100/pos_mode_count,return_number_pos[i]));
    
            logit('Negative Mode (Total: %s ):'%neg_mode_count);
            if neg_mode_count>0:
                for i in range(max_results_per_query):
                    logit('Correct within first\t%s:\tBest:\t%s%%\tWorst:\t%s%%\tCandidateCount:\t%s'%(i+1,best_results_neg[i]*100/neg_mode_count,worst_results_neg[i]*100/neg_mode_count,return_number_neg[i]));
    
            logit('Both Modes (Total: %s ):'%(pos_mode_count+neg_mode_count));
            if neg_mode_count>0 or pos_mode_count>0:
                for i in range(max_results_per_query):
                    logit('Correct within first\t%s:\tBest:\t%s%%\tWorst:\t%s%%\tCandidateCount:\t%s'%(i+1,(best_results_neg[i]+best_results_pos[i])*100/(neg_mode_count+pos_mode_count),(worst_results_neg[i]+worst_results_pos[i])*100/(neg_mode_count+pos_mode_count),return_number_pos[i]+return_number_neg[i]));
            
            #logit('Missed Peaks due to unknown mode: %s'%(total_peaks-neg_mode_count-pos_mode_count));
    
            logtime();        
            '''
            logit('Annotation score stats:');
    
            correct_top_mean_neg=0;
            correct_top_max_neg=0;
            
            
            correct_top_mean_pos=0;
            correct_top_max_pos=0;
            
    
            
            for peakindex in range(total_peaks):
                peak=peaks[peakindex];
                peakmode=0;
                if peak.parent_spectrum.parameters['mode']==1:
                    peakmode=1;
                elif peak.parent_spectrum.parameters['mode']==-1:
                    peakmode=-1;
                
                maxscore_mean=0;
                maxscore_max=0;
                correct_score_mean=0;
                correct_score_max=0;
                for annotation in peak.annotations:
                    if annotation.mean_score>maxscore_mean:
                        maxscore_mean=annotation.mean_score;
                    if annotation.max_score>maxscore_max:
                        maxscore_max=annotation.max_score;
                    if annotation.min_correct>-1:
                        correct_score_mean=annotation.mean_score;
                        correct_score_max=annotation.max_score;
                if correct_score_mean>=maxscore_mean:
                    correct_top_mean_pos+=1;
                if correct_score_max>=maxscore_max:
                    correct_top_max_pos+=1;
                    
            logit('Positive Mode');
            if pos_mode_count>0:
                logit('Correct best mean score annotation\t%s%%'%(correct_top_mean_pos*100/pos_mode_count));
                logit('Correct best max score annotation\t%s%%'%(correct_top_max_pos*100/pos_mode_count));
    
            logit('Negative Mode');
            if neg_mode_count>0:
                logit('Correct best mean score annotation\t%s%%'%(correct_top_mean_neg*100/neg_mode_count));
                logit('Correct best max score annotation\t%s%%'%(correct_top_max_neg*100/neg_mode_count));
                
            logit('Both Modes');
            if neg_mode_count>0 or pos_mode_count>0:
                logit('Correct best mean score annotation\t%s%%'%((correct_top_mean_pos+correct_top_mean_neg)*100/(pos_mode_count+neg_mode_count)));
                logit('Correct best max score annotation\t%s%%'%((correct_top_max_pos+correct_top_max_neg)*100/(pos_mode_count+neg_mode_count)));
    
            logtime();
            logit('Adduct Retrieval Stats');
    
            adduct_results={};        
            for peakindex in range(total_peaks):
                peak=peaks[peakindex];
                correct_adduct=peak.parameters['ion_type'];
    
                if not (correct_adduct in adduct_results):
                    adduct_results[correct_adduct]=ResBlock();
                else:
                    adduct_results[correct_adduct].count+=1;
                
                maxscore_mean=0;
                maxscore_max=0;
                correct_score_mean=0;
                correct_score_max=0;
                
                for annotation in peak.annotations:
                    if annotation.adduct.definition==correct_adduct:
                        if annotation.min_correct>-1:
                            adduct_results[correct_adduct].correct+=1;
                            for i in range(len(adduct_results[correct_adduct].best_results)):
                                if annotation.min_correct<=i:
                                    adduct_results[correct_adduct].best_results[i]+=1;
                                if annotation.max_correct<=i:
                                    adduct_results[correct_adduct].worst_results[i]+=1;
                
                
                    if annotation.mean_score>maxscore_mean:
                        maxscore_mean=annotation.mean_score;
                    if annotation.max_score>maxscore_max:
                        maxscore_max=annotation.max_score;
                    if annotation.min_correct>-1:
                        correct_score_mean=annotation.mean_score;
                        correct_score_max=annotation.max_score;
                if correct_score_mean>=maxscore_mean:
                    adduct_results[correct_adduct].correct_top_mean+=1;
                if correct_score_max>=maxscore_max:
                    adduct_results[correct_adduct].correct_top_max+=1;        
                        
            for adduct in adduct_results.keys():
                logit('Adduct type: %s'%adduct);    
                adduct_result=adduct_results[adduct];
                adduct_result.percent();
                logit('TotalCount: %s'%adduct_result.count);
                logit('Correct molecule in adduct: %s%%'%(adduct_result.correct));
                logit('Correct adduct in top mean: %s%%'%(adduct_result.correct_top_mean));
                logit('Correct adduct in top max: %s%%'%(adduct_result.correct_top_max));
                for i in range(max_results_per_query):
                    logit('Best/Worst for within the first\t%s:\t%s%%/%s%%'%(i+1,\
                        adduct_result.best_results[i],adduct_result.worst_results[i]));
            
            logit('Merging and sorting different adduct annotations...');
            logtime();
    
            best_results_pos=[0]*max_results_per_query;
            worst_results_pos=[0]*max_results_per_query;
            best_results_neg=[0]*max_results_per_query;
            worst_results_neg=[0]*max_results_per_query;
            
            for peakindex in range(total_peaks):
                peak=peaks[peakindex];
                peakmode=0;
                if peak.parent_spectrum.parameters['mode']==1:
                    peakmode=1;
                elif peak.parent_spectrum.parameters['mode']==-1:
                    peakmode=-1;
                shortinchi=peak.parent_spectrum.parameters['shortinchi'];
    
                merge_annotations(peak, remove_old_annotations=True, results_limit=max_results_per_query, total_score=total_multiplicative_score);
                annotation=peak.merged_annotations;
            
                min_correct=-1;
                max_correct=-1;
                for index in reversed(range(1,len(annotation.mol_candidates))):
                    if annotation.mol_candidates[index-1]['ShortInChi']==annotation.mol_candidates[index]['ShortInChi']:
                         del annotation.mol_candidates[index]; # For the purpose of statistics condensing sequential identical Short InChi-s
                            
                for index in range(len(annotation.mol_candidates)):
                        if annotation.mol_candidates[index]['ShortInChi']==shortinchi:
                            if min_correct==-1:
                                min_correct=index;
                            max_correct=index;
    
                if min_correct>-1:
                        if peakmode==1:
                            for i in range(len(best_results_pos)):
                                if min_correct<=i:
                                    best_results_pos[i]+=1;
                                if max_correct<=i:
                                    worst_results_pos[i]+=1;
                        elif peakmode==-1:
                            for i in range(len(best_results_neg)):
                                if min_correct<=i:
                                    best_results_neg[i]+=1;
                                if max_correct<=i:
                                    worst_results_neg[i]+=1;
    
            logit('Positive Mode (Total: %s ):'%pos_mode_count);
            if pos_mode_count>0:
                for i in range(max_results_per_query):
                    logit('Merged Correct within first\t%s:\tBest:\t%s%%\tWorst:\t%s%%'%(i+1,best_results_pos[i]*100/pos_mode_count,worst_results_pos[i]*100/pos_mode_count));
    
            logit('Negative Mode (Total: %s ):'%neg_mode_count);
            if neg_mode_count>0:
                for i in range(max_results_per_query):
                    logit('Merged Correct within first\t%s:\tBest:\t%s%%\tWorst:\t%s%%'%(i+1,best_results_neg[i]*100/neg_mode_count,worst_results_neg[i]*100/neg_mode_count));
    
            logit('Both Modes (Total: %s ):'%(pos_mode_count+neg_mode_count));
            if neg_mode_count>0 or pos_mode_count>0:
                for i in range(max_results_per_query):
                    logit('Merged Correct within first\t%s:\tBest:\t%s%%\tWorst:\t%s%%'%(i+1,(best_results_neg[i]+best_results_pos[i])*100/(neg_mode_count+pos_mode_count),(worst_results_neg[i]+worst_results_pos[i])*100/(neg_mode_count+pos_mode_count)));
            
            #logit('Finished analysing FragPrint only retrieval');
            logtime();
            '''
    




if __name__=='__main__':

    
    test_spectral_database_path='e:/Imperial/TestDB/NIST14_HMDB_MassBank_5fold_CV_FPT3';
    
    test_chemical_databases=['PubChem','MassBank','HMDB','ChEBI'];
    
    test_output_folder='e:/Imperial/TestDB/Output_multi_precision0.001_proton';
    
    fingerprint_model_path='e:/Imperial/TestDB/SVM_Train/Ported_csc/linear_noweight';    
    
    if not os.path.exists(test_output_folder):
        os.makedirs(test_output_folder);
    
    dbfile=os.path.abspath('../tests/databases.list');
    
    nbatch=5;
    nproc=6;
    ppm=20;
    max_results_per_query=100;

    
    
    print(dbfile);
    if not os.path.isfile(dbfile):
        print('No database file found!');
        quit();
        
        
    logfile=open(os.path.join(test_output_folder,'log.txt'),'w');


    print('Starting. %s'%time.strftime('%d/%m/%y %H:%M:%S'));
    starttime=time.time();


    logit('Starting. %s'%time.strftime('%d/%m/%y %H:%M:%S'));
    logit('Spectral path: %s'%test_spectral_database_path);
    logit('Chemical Databases Included: %s'%test_chemical_databases);
                     
    specmanager=SpectralManager();
    logit('Reading test spectra...');
    #print(test_spectral_database_path);
    specmanager.import_textfile_spectra_from_folder(test_spectral_database_path);
    
    logit('Finshed reading spectra. Total number: %s'%len(specmanager.ms_spectra));
    
    logtime();

    selected_InChi=load_inchi_list('e:/Imperial/HMDB/PubChem_Test/InChi_Selected100.list');

    #remove_spectra_not_in_inchi_list(selected_InChi,['HMDB','MassBank']);

    remove_low_precision_spectra(0.001);
    correct_precursor_mz(specmanager);    
    specmanager.export_textfile_spectra_to_folder('e:/Imperial/TestDB/Precise_0.001');
                    
    logit('Initializing batches');

    batch=[];
    for i in range(nbatch):
        batch.append([]);
    
    for i in range(len(specmanager.ms_spectra)):
          batch[int(specmanager.ms_spectra[i].parameters['crossvalidation_batch_index'])].append(i);
             
    print(batch);

    multiprocessing.freeze_support();

    with PeakAnnotator(dbfile, fingerprint_model_path, nproc) as peak_annotator:
        
            
        logit('dH for Fragmenter: [0], precision 0.001, analysis type 0');
        logtime();
        testrun(['FragPrintScorer'],[([0], 0.01, 1, 0)]);
    
        logit('dH for Fragmenter: [0,1], precision 0.001, analysis type 0');
        logtime();
        testrun(['FragPrintScorer'],[([0,1], 0.01, 1, 0)]);

        logit('dH for Fragmenter: [0,-1], precision 0.001, analysis type 0');
        logtime();
        testrun(['FragPrintScorer'],[([0,-1], 0.01, 1, 0)]);

        logit('dH for Fragmenter: [0,-1,1], precision 0.001, analysis type 0');
        logtime();
        testrun(['FragPrintScorer'],[([0,-1,1], 0.01, 1, 0)]);


        logit('dH for Fragmenter: [0], precision 0.001, analysis type 1');
        logtime();
        testrun(['FragPrintScorer'],[([0], 0.01, 1, 1)]);
    
        logit('dH for Fragmenter: [0,1], precision 0.001, analysis type 1');
        logtime();
        testrun(['FragPrintScorer'],[([0,1], 0.01, 1, 1)]);

        logit('dH for Fragmenter: [0,-1], precision 0.001, analysis type 1');
        logtime();
        testrun(['FragPrintScorer'],[([0,-1], 0.01, 1, 1)]);

        logit('dH for Fragmenter: [0,-1,1], precision 0.001, analysis type 1');
        logtime();
        testrun(['FragPrintScorer'],[([0,-1,1], 0.01, 1, 1)]);
        
        
        logit('dH for Fragmenter: [0], precision 0.001, analysis type 2');
        logtime();
        testrun(['FragPrintScorer'],[([0], 0.01, 1, 2)]);
    
        logit('dH for Fragmenter: [0,1], precision 0.001, analysis type 2');
        logtime();
        testrun(['FragPrintScorer'],[([0,1], 0.01, 1, 2)]);

        logit('dH for Fragmenter: [0,-1], precision 0.001, analysis type 2');
        logtime();
        testrun(['FragPrintScorer'],[([0,-1], 0.01, 1, 2)]);

        logit('dH for Fragmenter: [0,-1,1], precision 0.001, analysis type 2');
        logtime();
        testrun(['FragPrintScorer'],[([0,-1,1], 0.01, 1, 2)]);

             
    specmanager.close();
    
    logtime();
    
    logit('Finished');
    
    logfile.close();
    
