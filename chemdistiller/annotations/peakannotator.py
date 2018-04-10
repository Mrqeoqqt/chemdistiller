# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 13:10:14 2016

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

from chemdistiller.settings import test_as_single_process;
from chemdistiller.scorers.fragprint import FragPrintScorer;
from chemdistiller.scorers.fingerprint import FingerPrintScorer;
from chemdistiller.scorers.cfm import CFMScorer;
from chemdistiller.scorers.cfm2 import CFMScorer2;
from chemdistiller.chemdb.manager import DBManager;
from chemdistiller.utils.progress import ProgressDisplay;
from chemdistiller.scorers.csifingerid import CSIFingerIDScorer;

if not test_as_single_process:
    import multiprocessing;


from chemdistiller.scorers.totalscore import total_multiplicative_score;


def annotate_peakblock(args):
    global test_as_single_process;
    peaks, database_list_file, test_chemical_databases, scorers_list, scorers_settings, total_score, results_limit, save_memory, ppm, blockindex, progress_report_queue, required_fields=args;
    #Initialize chemical DB manager
    result=[];
    with DBManager(database_list_file) as dbmanager:
    
        db_indexes=dbmanager.db_indexes_from_db_names(test_chemical_databases,case_sensitive=True);
        
        #Initialize selected scorers
        available_scorers=[];
    
        if 'FragPrintScorer' in scorers_list:
            if scorers_settings:
                fragprint_scorer=FragPrintScorer(blockindex, scorers_settings[scorers_list.index('FragPrintScorer')]);
            else:
                fragprint_scorer=FragPrintScorer(blockindex, []);                
            available_scorers.append(fragprint_scorer);
    
        if 'FingerPrintScorer' in scorers_list:
            if scorers_settings:            
                fingerprint_scorer=FingerPrintScorer(blockindex, scorers_settings[scorers_list.index('FingerPrintScorer')]);
            else:
                fingerprint_scorer=FingerPrintScorer(blockindex, []);
                
            available_scorers.append(fingerprint_scorer);
            
        if 'CFMScorer' in scorers_list:
            if scorers_settings:
                CFM_scorer=CFMScorer(blockindex, scorers_settings[scorers_list.index('CFMScorer')]);
            else:
                CFM_scorer=CFMScorer(blockindex, []);                
            available_scorers.append(CFM_scorer);

        if 'CFMScorer2' in scorers_list:
            if scorers_settings:
                CFM_scorer2=CFMScorer2(blockindex, scorers_settings[scorers_list.index('CFMScorer2')]);
            else:
                CFM_scorer2=CFMScorer2(blockindex, []);                
            available_scorers.append(CFM_scorer2);
        
        if 'CSIFingerIDScorer' in scorers_list:
            if scorers_settings:
                CSI_scorer=CSIFingerIDScorer(blockindex, scorers_settings[scorers_list.index('CSIFingerIDScorer')]);
            else:
                CSI_scorer=CSIFingerIDScorer(blockindex, []);                
            available_scorers.append(CSI_scorer);

        # Score annotations in peaks
        peakcount=len(peaks);    
        for peakindex in range(peakcount):
            peak=peaks[peakindex];
            res=[];
            result.append(res);
            
            if test_as_single_process:
                print('Annotating %s of %s'%(peakindex+1,peakcount));
            elif not (progress_report_queue is None):
                progress_report_queue.put((blockindex,peakindex+1,peakcount));
                
            for annotation in peak.annotations:
                #Prepare scorers
                accepted=True;
                selected_scorers=[];
                for scorer in available_scorers:
                    if hasattr(scorer,'supported_adducts'):
                        if not (annotation.adduct.definition in scorer.supported_adducts):
                            accepted=False;
                            print('Unsupported adduct %s in scorer %s. Skipping annotation.'%(annotation.adduct.definition, scorer.name));
                            break;
                        
                    if scorer.configure_scorer(peak, annotation)==0:
                        selected_scorers.append(scorer);
                    else:
                        print('Error initializing scorer %s. Skipping annotation.'%scorer.name);
                        accepted=False;
                        break;
                #Prepare molecular search
                if accepted:
                    mz=annotation.adduct.get_mass(peak.mz);
                    charge=annotation.adduct.charge;
                    if not (annotation.formula_scorer is None):  #If formula scorer pre-defined - use it.
                        selected_scorers.append(annotation.formula_scorer);
                    if not (annotation.element_scorer is None):  #If element scorer pre-defined - use it.
                        selected_scorers.append(annotation.element_scorer);    
                    #Query chemical database
                    res.append(dbmanager.query_by_mz_scored(mz, ppm, charge, db_indexes,\
                        filters=annotation.filters, scorers=selected_scorers, required_fields=required_fields, results_limit=results_limit, save_memory=save_memory));
        for scorer in available_scorers:
            if isinstance(scorer, CFMScorer) or isinstance(scorer, CFMScorer2):
                print(scorer.empty_record_count, scorer.non_empty_record_count);
            scorer.close();
    return result;

class PeakAnnotator:

    def __init__(self, database_list_files, fingerprint_model_path, nproc=1):
        self.database_list_files=database_list_files;
        self.fingerprint_model_path=fingerprint_model_path;
        self.nproc=nproc;
        

    def annotate_peaks(self, peaks, test_chemical_databases, \
            scorers_list=[], scorers_settings=[], total_score=total_multiplicative_score,\
            required_fields=set(), results_limit=100, save_memory=False, ppm=20, overwrite=False, tidy_up=True):

        global test_as_single_process;

        if test_as_single_process:
            print('Preparing peak blocks...');

        self.peakblocks=[];
        for i in range(self.nproc):
            self.peakblocks.append([]);
        for i in range(len(peaks)):
            self.peakblocks[i%self.nproc].append(peaks[i]);
                
        self.test_chemical_databases=test_chemical_databases;
        self.scorers_list=scorers_list;
        self.scorers_settings=scorers_settings;
        self.total_score=total_score;
        self.results_limit=results_limit;
        self.save_memory=save_memory;
        self.ppm=ppm;
        self.overwrite=overwrite;
        self.peaks=peaks;
        self.required_fields=required_fields;

        if 'FingerPrintScorer' in self.scorers_list:
            #if test_as_single_process:
            print('Preparing FingerPrintScorer...');
            if self.scorers_settings:
                FingerPrintScorer.predict_fingerprints(self.peaks, self.fingerprint_model_path, self.scorers_settings[self.scorers_list.index('FingerPrintScorer')], self.nproc, self.overwrite);
            else:
                FingerPrintScorer.predict_fingerprints(self.peaks, self.fingerprint_model_path, [], self.nproc, self.overwrite);                
        
        if 'FragPrintScorer' in self.scorers_list:
            print('Preparing FragPrintScorer...');
            FragPrintScorer.generate_fragment_lists(self.peaks, self.nproc, self.overwrite);
            print('Completed FragPrintScorer preparations');
        
        #if 'CFMScorer' in self.scorers_list:
        #    print('Preparing CFMScorer...');
        #    CFMScorer.generate_fragment_lists(self.peaks, self.nproc, self.overwrite);
        #    print('Completed CFM preparations');
        
            
        if test_as_single_process:
            cc=0;
            for block in self.peakblocks:
                cc+=1;
                print('Block %s of %s'%(cc,len(self.peakblocks)));
                result=annotate_peakblock((block, self.database_list_files[0], self.test_chemical_databases,\
                    self.scorers_list, self.scorers_settings, self.total_score, self.results_limit, self.save_memory, \
                    self.ppm, cc, None, self.required_fields.copy()));
                for i in range(len(result)):
                    res=result[i];
                    for j in range(len(res)):
                        block[i].annotations[j].mol_candidates=res[j];
                
        else:
            
            progress_display=ProgressDisplay('Retrieving and scoring candidates...', self.nproc);
            pool=multiprocessing.Pool(self.nproc);

            args=[];
            m = multiprocessing.Manager();
            progress_report_queue = m.Queue();

            for blockindex in range(self.nproc):
                block=self.peakblocks[blockindex];
                db_index=blockindex%len(self.database_list_files);
                args.append((block, self.database_list_files[db_index], self.test_chemical_databases,\
                    self.scorers_list, self.scorers_settings, self.total_score, self.results_limit, self.save_memory, \
                    self.ppm, blockindex, progress_report_queue, self.required_fields));
                    
            poolresult=pool.map_async(annotate_peakblock,args);
            
            while not poolresult.ready():
                if not progress_report_queue.empty():
                    progress_display.display_progress_from_queue(progress_report_queue.get());
            
    
            result=poolresult.get();
            pool.close();
            pool.join();
    
            for i in range(len(result)):
                block=result[i];
                peakblock=self.peakblocks[i];
                for j in range(len(block)):
                    respeak=block[j];
                    peak=peakblock[j];
                    for k in range(len(respeak)):
                        peak.annotations[k].mol_candidates=respeak[k];
                
        if tidy_up:
            print('\nTidying up....');
            
            if 'FingerPrintScorer' in self.scorers_list:
                FingerPrintScorer.tidy_up_peaks(peaks, self.required_fields);
            
            if 'FragPrintScorer' in self.scorers_list:
                FragPrintScorer.tidy_up_peaks(peaks, self.required_fields);
                
            if 'CFMScorer' in self.scorers_list:
                CFMScorer.tidy_up_peaks(peaks, self.required_fields);    
               
            if 'CFMScorer2' in self.scorers_list:
                CFMScorer2.tidy_up_peaks(peaks, self.required_fields);       
            print('Finished');
            
        
        
    def close(self):
        return 0;
        
    def __enter__(self):
        return self;

    def __exit__(self, exc_type, exc_value, traceback):
        self.close();


if __name__ == '__main__':
    if not test_as_single_process:
        multiprocessing.freeze_support();
