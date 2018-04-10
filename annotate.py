# -*- coding: utf-8 -*-
"""
==============================================================================
        ChemDistiller Annotator Script
==============================================================================
    This script is designed for a simple independent annotation of tandem MS 
spectra.


Copyright:                      Imperial College London
Project:                        MetaSpace
Funding:                        Horizon 2020
License:                        BSD
Chief project investigator:     Dr. Kirill Veselkov
Lead developer:                 Dr. Ivan Laponogov

"""


import sys;
import os;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.dirname(os.path.realpath(__file__));
    sys.path.append(chemdistiller_path);

import time;

from chemdistiller.settings import test_as_single_process, registered_scorers, registered_filters;
from chemdistiller.msspectra.manager import SpectralManager;
from chemdistiller.msspectra.adducts import get_adduct_by_name, get_positive_mode_adducts, get_negative_mode_adducts,\
                                            global_supported_adducts;
from chemdistiller.annotations.peakannotator import PeakAnnotator;
from chemdistiller.scorers.totalscore import total_multiplicative_score;
from chemdistiller.scorers.fragprint import FragPrintScorer;
from chemdistiller.scorers.formula import FormulaScorer;
from chemdistiller.scorers.element import ElementScorer;
from chemdistiller.scorers.fingerprint import FingerPrintScorer;
from chemdistiller.filters.formula import FormulasFilter;
from chemdistiller.filters.inchi import InChIFilter;
from chemdistiller.filters.element import ElementCompositionFilter;
from chemdistiller.io.jsonio import spectra_to_json;
from chemdistiller.io.html_report import generate_HTML_report;
from chemdistiller.utils.periodictable import Formula;
from chemdistiller.utils.inchi import get_short_inchi_from_full_inchi;

from chemdistiller.annotations.peakannotation import MSPeakAnnotation, merge_annotations;

import multiprocessing;

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
    


if __name__=='__main__':
    
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

    
    
def run_annotation(spectral_input_path, db_file, chemical_databases, output_folder, SVMs_path, ncpu, ppm, max_results_per_query, test_mode=False, ignore_peaks_without_spectra=True):

    global starttime;
    global logfile;

    if not os.path.isfile(db_file):
        print('No database list file found!');
        quit();

    db_files=[db_file];
    
    if not os.path.exists(output_folder):
        try:
            os.makedirs(output_folder);
        except:
            print('Cannot create output path %s !'%output_folder);
            print('Run: python annotate.py -h for help');
            sys.exit(2);
    
    logfile=open(os.path.join(output_folder,'log.txt'),'w');
    starttime=time.time();

    logit('Starting. %s'%time.strftime('%d/%m/%y %H:%M:%S'));
    logit('Input Spectral path: %s'%spectral_input_path);
    logit('Output path: %s'%output_folder)
    if chemical_databases:
        logit('Chemical Databases Included: %s'%chemical_databases);
    else:
        logit('Chemical Databases Included: All');
        
    logit('from chemical database paths file %s'%db_file);
    logit('SVMs from %s'%SVMs_path);
    logit('nCPUs to use: %s'%ncpu);
    logit('Default ppm tolerance for MS1: %s'%ppm);    

    
    specmanager=SpectralManager();
    logit('Reading test spectra...');
        
    specmanager.import_textfile_spectra_from_folder(spectral_input_path);
        
    logit('Finshed reading spectra. Total number: %s'%len(specmanager.ms_spectra));
        
    logtime();
    
    multiprocessing.freeze_support();
    
    scorers_list=['FingerPrintScorer','FragPrintScorer'];
    scorers_settings=[[],[]];
    total_score=total_multiplicative_score;
        
    with PeakAnnotator(db_files, SVMs_path, ncpu) as peak_annotator:
        logit('Beginning annotation...');
        logtime();
        
        total_peaks=0;
    
        best_results_neg=[0]*max_results_per_query;
        worst_results_neg=[0]*max_results_per_query;
        return_number_neg=[0]*max_results_per_query;
        neg_mode_count=0;
            
        best_results_pos=[0]*max_results_per_query;
        worst_results_pos=[0]*max_results_per_query;
        return_number_pos=[0]*max_results_per_query;
        pos_mode_count=0;
        
        peaks=[];
        
        supported_adducts=set(global_supported_adducts);
                                
        #Collect all supported adducts
        #if 'FingerPrintScorer' in scorers_list:
        #    supported_adducts=supported_adducts|FingerPrintScorer.supported_adducts;

        #if 'FragPrintScorer' in scorers_list:
        #    supported_adducts=supported_adducts|FragPrintScorer.supported_adducts;

        #Leave only adducts supported by all selected filters                                
        if 'FingerPrintScorer' in scorers_list:
            supported_adducts=supported_adducts&FingerPrintScorer.supported_adducts;

        if 'FragPrintScorer' in scorers_list:
            supported_adducts=supported_adducts&FragPrintScorer.supported_adducts;

                
        logit('Preparing annotations');                        
        if test_mode:
            logit('Running in test mode.....');            
            
            logit('Assuming correct adduct info and 0 isotope.');
            logit('No formula assumption. Mass precision: %s ppm'%ppm);
            logit('Only considering [M+H]+ and [M-H]- adducts for now');
            
            print('\n');
                
            for spectrum_index in range(len(specmanager.ms_spectra)):
                    #print('\r%s of %s spectra preprocessed'%(spectrum_index,len(specmanager.ms_spectra)));
                    spectrum=specmanager.ms_spectra[spectrum_index];
                    
                    for peak in spectrum.peaks:
                        if (not ignore_peaks_without_spectra) or (hasattr(peak,'ms_spectra') and peak.ms_spectra):
                            if hasattr(peak,'parameters') and ('ion_type' in peak.parameters):
                                if (peak.parameters['ion_type'] in FragPrintScorer.supported_adducts):
                                    peak.ppm=ppm; 
                                    
                                    shortinchi=get_short_inchi_from_full_inchi(peak.parent_spectrum.parameters['inchi']);
                                    
                                    annotation=MSPeakAnnotation(peak, adduct=\
                                        get_adduct_by_name(peak.parameters['ion_type'], \
                                        spectrum.parameters['mode']), isotope=0, \
                                        formula_scorer=None, element_scorer=None, filters=[], scores={});
                                        
                                    #testformulascorer=FormulaScorer();
                                    #testformulascorer.setup_scorer(peak.parent_spectrum.parameters['formula'],1.0);
                                    #testelementscorer=ElementScorer();
                                    #testelementscorer.setup_scorer({'C':0.8,'O':0.9,'Si':0.1});
                                    
                                    #formulafilter=FormulasFilter(formulas=[peak.parent_spectrum.parameters['formula']]);
                                    
                                    #elementfilter=ElementCompositionFilter(peak.parent_spectrum.parameters['formula'], peak.parent_spectrum.parameters['formula']);
                                    
                                    #inchifilter=InChIFilter(ref_inchi=peak.parent_spectrum.parameters['inchi'],use_short_inchi=False, match_type=0)
                                                                        
                                    
                                    #annotation=MSPeakAnnotation(peak, adduct=\
                                    #    get_adduct_by_name(peak.parameters['ion_type'], \
                                    #    spectrum.parameters['mode']), isotope=0, \
                                    #    formula_scorer=testformulascorer, element_scorer=testelementscorer, \
                                    #    filters=[formulafilter, inchifilter, elementfilter], scores={'ExtraScore':25.0, 'AllExtra':1.9});
                                        
                                    
                                    if not (annotation.adduct is None):
                                        peak.annotations=[annotation];
                                        peaks.append(peak);
                                    else:
                                        logit('Unsupported Adduct: %s'%peak.parameters['ion_type']);
        else:
            

            logit('Running in normal mode.....');            
            logit('Mass precision: %s ppm'%ppm);
            
            print('\n');
            
            for spectrum_index in range(len(specmanager.ms_spectra)):
                    #print('\r%s of %s spectra preprocessed'%(spectrum_index,len(specmanager.ms_spectra)));
                    spectrum=specmanager.ms_spectra[spectrum_index];
                    for peak in spectrum.peaks:
                        if (not ignore_peaks_without_spectra) or (hasattr(peak,'ms_spectra') and peak.ms_spectra):
                            peak.ppm=ppm; 
                            generate_annotation=True;
                            if hasattr(peak, 'annotations'):
                                if peak.annotations:
                                    generate_annotation=False;
                                    peaks.append(peak);
                            
                            if generate_annotation:
                                #Assuming isotope 0 by default
                                if spectrum.parameters['mode']==1:
                                    selected_adducts=get_positive_mode_adducts(supported_adducts);
                                else:
                                    selected_adducts=get_negative_mode_adducts(supported_adducts);
                                
                                peak.annotations=[];
                                for adduct in selected_adducts:
                                    annotation=MSPeakAnnotation(peak, adduct=adduct, isotope=0, formula_scorer=None, element_scorer=None, filters=[], scores={});
                                    peak.annotations.append(annotation);
                                    peaks.append(peak);
                                    
        
        logit('Peaks to annotate: %s'%len(peaks));
        
        if test_mode:
            required_fields=set(['ShortInChI','InChI','SMILES','Formula','IDs']);
            #required_fields=set(['ShortInChI','InChI','SMILES','Formula','IDs', 'InChIKeyValues', 'InChiKey', 'FormulaVector', 'ElementVector', 'Frag', 'FPT']);
            
        else:
            required_fields=set(['InChI','SMILES','Formula','IDs']);
        
        peak_annotator.annotate_peaks(peaks, chemical_databases, \
                scorers_list=scorers_list, scorers_settings=scorers_settings, total_score=total_score, required_fields=required_fields,\
                results_limit=max_results_per_query, save_memory=False, ppm=ppm, overwrite=True);
        
        
        
        annotation_count=0;
        total_candidates=0;                    

        total_peaks=len(peaks);
        for peakindex in range(total_peaks):
            peak=peaks[peakindex];
            if peak.annotations:
                for annotation in peak.annotations:
                    if not (annotation.mol_candidates is None):
                        annotation_count+=1;
                        total_candidates+=annotation.mol_candidates.total_candidate_count;
                        
        
        
        logit('Finished annotating... ');
        logit('Annotated peaks :%s'%total_peaks);    
        if total_peaks>0:
            logit('Averaged annotations per peak:%s'%(float(annotation_count)/total_peaks));    
            logit('Averaged candidate molecules per peak:%s'%(float(total_candidates)/total_peaks));    
        
        
        #Calculating retrieval statistics using test data
        #=====================================================================
        
        if test_mode:
            #logit('Removing peaks with no annotations found from consideration...');
            #for i in reversed(range(len(peaks))):
            #    if not peaks[i].annotations:
            #        del peaks[i];
                
            #logit('Finished removing... Annotated peaks :%s'%len(peaks));    
                
            total_peaks=len(peaks);
            print('Calculating retrieval stats...\n');
                
            for peakindex in range(total_peaks):
                peak=peaks[peakindex];
                peakmode=0;
                if peak.parent_spectrum.parameters['mode']==1:
                    pos_mode_count+=1;
                    peakmode=1;
                elif peak.parent_spectrum.parameters['mode']==-1:
                    neg_mode_count+=1;
                    peakmode=-1;
                        
                #print('\r%s of %s '%(peakindex,total_peaks));
                shortinchi=get_short_inchi_from_full_inchi(peak.parent_spectrum.parameters['inchi']);
        
                for annotation in peak.annotations:
                    annotation.min_correct=-1;
                    annotation.max_correct=-1;
                    annotation.mean_score=0.0;
                    annotation.max_score=0.0;
                    for index in range(len(annotation.mol_candidates.mol_list)):
                        total_score=annotation.mol_candidates.mol_list[index]['TotalScore'];
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
                        if annotation.mol_candidates.mol_list[index-1]['ShortInChI']==annotation.mol_candidates.mol_list[index]['ShortInChI']:
                            del annotation.mol_candidates.mol_list[index]; # For the purpose of statistics condensing sequential identical Short InChI-s
                        
                    for index in range(len(annotation.mol_candidates.mol_list)):
                        if annotation.mol_candidates.mol_list[index]['ShortInChI']==shortinchi:
                            annotation.mol_candidates.mol_list[index]['Annotation']='Correct';
                            if annotation.min_correct==-1:
                                annotation.min_correct=index;
                            annotation.max_correct=index;
                        else:
                            annotation.mol_candidates.mol_list[index]['Annotation']='Wrong';
                            
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
    
            logit('Finished Calculating Retrieval Stats.');
            logit('Positive Mode Count: %s'%pos_mode_count);
            logit('Negative Mode Count: %s'%neg_mode_count);
            logit('Total Peaks: %s'%total_peaks);
        
            logit('Positive Mode (Total: %s ):'%pos_mode_count);
            if pos_mode_count>0:
                for i in range(max_results_per_query):
                    logit('Correct within first\t%s:\tBest:\t%.2f%%\tWorst:\t%.2f%%\tCCount:\t%s'%(i+1,float(best_results_pos[i])*100/pos_mode_count,worst_results_pos[i]*100/pos_mode_count,return_number_pos[i]));
        
            logit('Negative Mode (Total: %s ):'%neg_mode_count);
            if neg_mode_count>0:
                for i in range(max_results_per_query):
                    logit('Correct within first\t%s:\tBest:\t%.2f%%\tWorst:\t%.2f%%\tCCount:\t%s'%(i+1,float(best_results_neg[i])*100/neg_mode_count,worst_results_neg[i]*100/neg_mode_count,return_number_neg[i]));
        
            logit('Both Modes (Total: %s ):'%(pos_mode_count+neg_mode_count));
            if neg_mode_count>0 or pos_mode_count>0:
                for i in range(max_results_per_query):
                    logit('Correct within first\t%s:\tBest:\t%.2f%%\tWorst:\t%.2f%%\tCCount:\t%s'%(i+1,float(best_results_neg[i]+best_results_pos[i])*100/(neg_mode_count+pos_mode_count),float(worst_results_neg[i]+worst_results_pos[i])*100/(neg_mode_count+pos_mode_count),return_number_pos[i]+return_number_neg[i]));
            for peak in peaks:
                merge_annotations(peak, remove_old_annotations=False);

        #Finished calculating retrieval statistics using test data
        #=====================================================================
        
    logtime();
    logit('Finished Annotation. Exporting results...');
    
    logtime();
    logit('Exporting to JSON...');
        
    spectra_to_json(os.path.join(output_folder,'annotated_spectra.json'), specmanager.ms_spectra);
    
    logtime();
    logit('Exporting to internal text format...');
    
    
    specmanager.export_textfile_spectra_to_folder(os.path.join(output_folder,'annotated_spectra'));
    
    logit('Preparing HTML report...');
    
    generate_HTML_report(os.path.join(output_folder,'Report'), specmanager);
                
    specmanager.close();
    
    logtime();
    
    logit('Finished');
    
    logfile.close();
    
        


if __name__=='__main__':
    spectral_input_path=os.path.join(chemdistiller_path,'data/MassBankTestSpectra');
    #spectral_input_path=os.path.join(chemdistiller_path,'data/short');
    chemical_databases=[];
    output_folder=os.path.join(chemdistiller_path,'testlogs/%s'%time.strftime("[%Y_%m_%d#%H_%M_%S]"));
    SVMs_path=os.path.join(chemdistiller_path,'SVMs');
    db_file=os.path.join(chemdistiller_path,'tests/databases_DEMO.list');
    #print(db_file)
    ncpu=1;
    ppm=20.0;
    max_results_per_query=10;
    max_cpu=multiprocessing.cpu_count();
    
    #print('Available Scorers: %s'%sorted(registered_scorers.keys()));
    #print('Available Filters: %s'%sorted(registered_filters.keys()));
    
    #sys.exit();
    
    program_description='\
==============================================================================\n\
        ChemDistiller Annotator Script\n\
==============================================================================\n\
    This script is designed for a simple independent annotation of tandem MS \n\
spectra. \n\
\n\
Currently implemented scorers:\n\n%s\n\nCurrently implemented filters:\n\n%s'%(\
"\n".join(sorted(registered_scorers.keys())),"\n".join(sorted(registered_filters.keys())));


    
    description_epilog='\n\
==============================================================================\n\
Copyright:                      Imperial College London\n\
Project:                        MetaSpace\n\
Funding:                        Horizon 2020\n\
License:                        BSD\n\
Chief project investigator:     Dr. Kirill Veselkov\n\
Lead developer:                 Dr. Ivan Laponogov\n\
\n\
\n\
'


    import argparse;
    parser = argparse.ArgumentParser(description=program_description, \
                                     formatter_class=argparse.RawDescriptionHelpFormatter, \
                                     epilog=description_epilog);
  
    
    #positional arguments    

    parser.add_argument("input_spectra_folder", \
        help="Path to the folder containing MS spectra to be annotated (internal text format). Default for testing:\n%s"%spectral_input_path, \
        type=str, \
        nargs='?', \
        default='');
        
    parser.add_argument("output_folder", \
        help="Path to the output folder for annotated spectra. Default for testing:\n%s"%os.path.join(chemdistiller_path,'testlogs/[current date_time]'),\
        type=str, \
        nargs='?', \
        default='');
    
    #optional arguments/flags
    parser.add_argument("-t", "--test", \
        help='Run the script in test mode using default parameters and \
            supplemented test spectra and databases. In this case input spectra and \
            output folder default to "<Main>/data/MassBankTestSpectra" and \
            "<Main>/testlogs/[current date_time]", where <Main> is the folder where your \
            ChemDistiller is located (Currently: %s). \
            Options that are supplied with --test option can be used to override \
            default database selection and SVM selection. Note! Before using this option if relying on the default SVMs \
            make sure that folder "<Main>/SVMs" contains a selected unpacked SVM model \
            (see "Install" section of the main "Readme" file of ChemDistiller)!\
            Packed SVM models are currently available from: \
            https://www.mediafire.com/folder/v4lb8s2nns9c6/SVMs\
            They are not included by default as a part of this package due to their \
            large size. Example: python annotate.py --test'%chemdistiller_path,\
        action='store_true', \
        default=False);
        
    parser.add_argument('-d', '--database_list_file', \
        help='override the default database \
            list file "<Main>/tests/databases_DEMO.list" (text file listing folders \n\
            containing compound databases) with the provided DATABASE_LIST_FILE',\
        type=str, \
        default=db_file);

    parser.add_argument('-s', '--svm_folder', \
        help='Override the default folder for the SVM models used\n\
            by Fingerprint filter (default - "<Main>/SVMs" ).\n\
            Make sure this folder contains unpacked SVM model set prior to running the \n\
            script!',\
        type=str, \
        default=SVMs_path);

    parser.add_argument('--databases', \
        help='Provide selective list of the chemical \n\
            databases to be used (default - all available). List is case sensitive and comma-\n\
            separated. Example: \n\
            python annotate.py --databases ChEBI,HMDB,MassBank <input_spectra_folder> <output_folder>\n',\
        type=str, \
        default='');

    parser.add_argument('-p', '--ncpu', \
        help="Provide the number of parallel processes to use (number\n\
            of CPUs). Default=%s, Maximum for this PC: %s\n"%(ncpu, max_cpu),\
        type=int, \
        default=ncpu);

    parser.add_argument('-m', '--delta_mz', \
        help="Provide 1D MS m/z tolerance for chemical DB \n\
            search, ppm. Default=%s ppm\n"%ppm,\
        type=float, \
        default=ppm);

    parser.add_argument('-r', '--max_results', \
        help="Provide maximum number of results to be returned\n\
            per annotation candidate. Default=%s (suitable for small chemical DBs)\n\
            Use -r 0 or for no limit.\n"%max_results_per_query,\
        type=int, \
        default=max_results_per_query);

        
    args = parser.parse_args();


    if args.test:
        if args.input_spectra_folder=='':
            args.input_spectra_folder=spectral_input_path;
        if args.output_folder=='':
            args.output_folder=output_folder;

    if args.input_spectra_folder=='' or args.output_folder=='':
        print('Error: not enough positional arguments!');
        parser.print_help();
        sys.exit(2);

    db_file=os.path.abspath(args.database_list_file);
    if not os.path.isfile(db_file):
        print('Error! No database list file found: %s!'%db_file);
        sys.exit(2);
            
    SVMs_path=os.path.abspath(args.svm_folder);
    if not os.path.exists(SVMs_path):
        print('No SVM folder found: %s!'%SVMs_path);
        sys.exit(2);

    if args.databases=='':
        chemical_databases=[];
    else:
        chemical_databases=args.databases.split(',');
    for i in range(len(chemical_databases)):
        chemical_databases[i]=chemical_databases[i].strip();

    ncpu=args.ncpu;
    if ncpu>max_cpu or ncpu<1:
        print('Error! NCPU can be between 1 and %s only!'%max_cpu);
        sys.exit(2);

    ppm=args.delta_mz;
    if ppm<=0:
        print('Error! DELTA_MZ should be greater than 0.0!');
        sys.exit(2);
        

    max_results_per_query=args.max_results;

    if args.test:
        print('Running in Test Mode');
        run_annotation(args.input_spectra_folder, db_file, chemical_databases, args.output_folder, SVMs_path, ncpu, ppm, max_results_per_query, True);
    else:
        spectral_input_path=os.path.abspath(args.input_spectra_folder);
        output_folder=os.path.abspath(args.output_folder);
        if not os.path.exists(spectral_input_path):
            print('Error! Input spectral folder %s is not found!'%spectral_input_path);
            sys.exit(2);

        run_annotation(spectral_input_path, db_file, chemical_databases, output_folder, SVMs_path, ncpu, ppm, max_results_per_query, False);
        
    print('Exiting');    
    
