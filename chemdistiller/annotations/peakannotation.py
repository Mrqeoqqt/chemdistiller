# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:44:01 2016

@author: Dr. Ivan Laponogov
"""

import os;
import sys;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);

from operator import itemgetter;

from chemdistiller.msspectra.adducts import get_adduct_by_name;
from chemdistiller.scorers.formula import FormulaScorer;
from chemdistiller.scorers.element import ElementScorer;
from chemdistiller.filters.manager import get_filter_by_name;
from chemdistiller.chemdb.results import DBQueryResult;

class MSPeakAnnotation:
    def __init__(self, parent_peak, adduct=None, isotope=0, formula_scorer=None, element_scorer=None, filters=None, scores=None):
        self.mol_candidates=None;
        self.adduct=adduct;
        self.isotope=isotope;
        self.isotope_extra_mass=0.0;
        self.formula_scorer=formula_scorer;
        self.element_scorer=element_scorer;
        if not (filters is None):
            self.filters=filters;
        else:
            self.filters=[];
        if not (scores is None):
            self.scores=scores;
        else:
            self.scores={};
        self.parent_peak=parent_peak;
        self.parameters={};
    
    def _pipe_from_textfile(self, finp):
        while True:
            s=finp.readline();
            if s=='':
                 return;
            s=s.rstrip('\n').lstrip();
            if '##' in s:
                s=s[:s.index('##')];
            
            if '=' in s:
                s=s.split('=', 1);
                if s[0].lower().startswith('adduct'):
                    self.adduct=get_adduct_by_name(s[1]);
                elif s[0].lower().startswith('isotope_extra_mass'):
                    self.isotope_extra_mass=float(s[1]);
                elif s[0].lower().startswith('isotope'):
                    self.isotope=int(s[1]);
                else:
                    self.parameters[s[0].lower()]=s[1];
            elif s.lower().startswith('formula_scorer'):
                self.formula_scorer=FormulaScorer();
                self.formula_scorer._pipe_from_textfile(finp);

            elif s.lower().startswith('element_scorer'):
                self.element_scorer=ElementScorer();
                self.element_scorer._pipe_from_textfile(finp);

            elif s.lower().startswith('mol_candidates'):
                self.mol_candidates=DBQueryResult();
                self.mol_candidates._pipe_from_textfile(finp);

            elif s.lower().startswith('scores'):
                finished=False;
                while not finished:
                    s=finp.readline();
                    if s=='':
                         return;
                    s=s.rstrip('\n').lstrip();
                    if '##' in s:
                        s=s[:s.index('##')];
                    if '=' in s:
                        s=s.split('=', 1);
                        self.scores[s[0]]=float(s[1]);
                    elif (s.lower().startswith('end')):
                        finished=True;

            elif s.lower().startswith('filters'):
                finished=False;
                while not finished:
                    s=finp.readline();
                    if s=='':
                         return;
                    s=s.rstrip('\n').lstrip();
                    if '##' in s:
                        s=s[:s.index('##')];
                    if '=' in s:
                        s=s.split('=', 1);
                        chemfilter=get_filter_by_name(s[1]);
                        chemfilter._pipe_from_textfile(finp);
                        self.filters.append(chemfilter);
                    elif (s.lower().startswith('end')):
                        finished=True;

            elif s.lower().startswith('end'):
                return;
    
    def _pipe_to_textfile(self, fout, indent=''):

        subindent='%s\t'%indent;
        for key in sorted(self.parameters.keys()):
            fout.write('%s%s=%s\n'%(indent,key,self.parameters[key]));
        
        if not (self.adduct is None):
            fout.write('%s%s=%s\n'%(indent,'adduct',self.adduct.definition));

        fout.write('%s%s=%s\n'%(indent,'isotope',self.isotope));
        fout.write('%s%s=%s\n'%(indent,'isotope_extra_mass',self.isotope_extra_mass));

        if not (self.formula_scorer is None):
            fout.write('%sformula_scorer\n'%indent);
            self.formula_scorer._pipe_to_textfile(fout, '%s\t'%indent);
            fout.write('%send\n'%indent);

        if not (self.element_scorer is None):
            fout.write('%selement_scorer\n'%indent);
            self.element_scorer._pipe_to_textfile(fout, '%s\t'%indent);
            fout.write('%send\n'%indent);
        
        if self.scores:
            fout.write('%sscores\n'%indent);
            for key in sorted(self.scores.keys()):
                fout.write('%s%s=%s\n'%(subindent,key,self.scores[key]));
            fout.write('%send\n'%indent);
                
        if self.filters:    
            fout.write('%sfilters\n'%indent);
            for filterobject in self.filters:
               fout.write('%s%s=%s\n'%(subindent,'filter',filterobject.name));
               filterobject._pipe_to_textfile(fout,'%s\t\t\t'%subindent);
               fout.write('%send\n'%subindent);
            fout.write('%send\n'%indent);
        
        if not (self.mol_candidates is None):
            fout.write('%s%s\n'%(indent,'mol_candidates'));
            self.mol_candidates._pipe_to_textfile(fout,'%s\t'%indent);
            fout.write('%send\n'%indent);


        
        
def merge_annotations(peak, remove_old_annotations=True, overwrite_merged_annotations=True, results_limit=100, total_score=None):
    #check if mergeable, i.e. all scores are available for all candidates.
    setofscores=set();
    for annotation in peak.annotations:
        if not (annotation.mol_candidates is None):
            if annotation.mol_candidates.mol_list:
                setofscores=setofscores|set(annotation.mol_candidates.mol_list[0]['Scores'].keys());
            setofscores=setofscores|set(annotation.scores.keys());
    mergeable=True;
    for annotation in peak.annotations:
        if not (annotation.mol_candidates is None):
            if annotation.mol_candidates.mol_list:
                combiscore=set(annotation.mol_candidates.mol_list[0]['Scores'].keys())|set(annotation.scores.keys());
                if setofscores!=combiscore:
                    mergeable=False;
                    break;

    if mergeable==False:
        raise TypeError('Annotations are not mergeable - different annotations have different sets of scores!');
    
    if (not hasattr(peak,'merged_annotations')) or overwrite_merged_annotations:
        peak.merged_annotations=[];
    
    for annotation in peak.annotations:
        if not (annotation.mol_candidates is None):
            for candidate in annotation.mol_candidates.mol_list:
                peak.merged_annotations.append(candidate);
                candidate['Adduct']=annotation.adduct.definition;
                candidate['Isotope']=annotation.isotope;
                candidate['IsotopeExtraMass']=annotation.isotope_extra_mass;
                target_scores=candidate['Scores'];
                for key in annotation.scores.keys():
                    target_scores[key]=annotation.scores[key];
                if not (total_score is None):
                    candidate['TotalScore']=total_score(target_scores);
    peak.merged_annotations=sorted(peak.merged_annotations, key=itemgetter('TotalScore'), reverse=True);
    
    if len(peak.merged_annotations)>results_limit and results_limit>0:
        del peak.merged_annotations[results_limit:];
    if remove_old_annotations:
        del peak.annotations;