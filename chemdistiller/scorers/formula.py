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

from chemdistiller.utils.periodictable import parse_formula, encode_formula_to_array, decode_formula_from_array;
import numpy as np;
from chemdistiller.settings import test_as_single_process, registered_scorers, registered_filters;

class FormulaScorer:
    name='FormulaScorer';

    def __init__(self, use_vector_form=True):
        self.vector_form=use_vector_form;
        if self.vector_form:
            self.required_fields=['FormulaVector'];
            #self.process_molecular_candidate_record=self.__process_molecular_candidate_record_vector;
        else:
            self.required_fields=['Formula'];
            #self.process_molecular_candidate_record=self.__process_molecular_candidate_record_formula;
        self.formulas=[];
        self.scores=[];
        self.unknown_score=0.0;
        
    def process_molecular_candidate_record(self, source_data_base, record):
        if self.vector_form:
            self.__process_molecular_candidate_record_vector(source_data_base, record);
        else:
            self.__process_molecular_candidate_record_formula(source_data_base, record);
            

    def _pipe_from_textfile(self, finp):
         while True:
            s=finp.readline();
            if s=='':
                 return;
            s=s.rstrip('\n').lstrip();
            if '##' in s:
                s=s[:s.index('##')];
            if '=' in s:
                s=s.split('=',1);
                if s[0].lower().startswith('formula'):
                    s=s[1].split(',');
                    self.formulas.append(encode_formula_to_array(parse_formula(s[0])));
                    self.scores.append(float(s[1]));
                elif s[0].lower().startswith('unknown_score'):
                    self.unknown_score=float(s[1]);
                    
            elif s.lower().startswith('dict_form'):
                self.vector_form=False;
            
            elif s.lower().startswith('end'):
                if not self.vector_form:
                    self.required_fields=['Formula'];
                    for i in range(len(self.formulas)):
                        self.formulas[i]=decode_formula_from_array(self.formulas[i]);
                    self.process_molecular_candidate_record=self.__process_molecular_candidate_record_formula;

                return;
    
    def _pipe_to_textfile(self, fout, indent=''):                
        if self.vector_form:
            fout.write('%svector_form\n'%(indent));
        else:
            fout.write('%sdict_form\n'%(indent));
        
        for i in range(len(self.formulas)):
            if self.vector_form:
                fout.write('%s%s=%s,%s\n'%(indent,'formula', decode_formula_from_array(self.formulas[i]).formula_to_string(), self.scores[i]));
            else:
                fout.write('%s%s=%s,%s\n'%(indent,'formula', self.formulas[i].formula_to_string()), self.scores[i]);

        fout.write('%s%s=%s\n'%(indent,'unknown_score',self.unknown_score));


    def __repr__(self):
        l=[];
        for i in range(len(self.formulas)):
            if self.vector_form:
                print(self.scores[i]);
                
                l.append('"%s"=%s'%(decode_formula_from_array(self.formulas[i]).formula_to_string(), self.scores[i]));
            else:
                l.append('"%s"=%s'%(self.formulas[i].formula_to_string(), self.scores[i]));
        s=', '.join(l);
        return 'FormulaScorer(%s)'%s;
    

    def close(self):
        return 0;
        
    def configure_scorer(self, peak, annotation): 
        #This scorer is annotation-independent
        return 0;

    def setup_scorer(self, formulas, scores, unknown_score=0.0):
        self.formulas=[];
        self.scores=[];
        self.unknown_score=unknown_score;
        if isinstance(formulas, str):
            formulas=formulas.split(',');
            for formula in formulas:
                self.formulas.append(encode_formula_to_array(parse_formula(formula)));
        elif isinstance(formulas, list):
            for formula in formulas:
                if isinstance(formula,dict):
                    self.formulas.append(encode_formula_to_array(formula));
                elif isinstance(formula, np.ndarray):
                    self.formulas.append(formula);
                elif isinstance(formula,str):
                    fs=formula.split(',');
                    for f in fs:
                        self.formulas.append(encode_formula_to_array(parse_formula(f)));
                else:
                    raise TypeError('Wrong type argument (formulas) for FormulaVectors initialization! str, formula, list of (formula, formulavector or str) supported only!');
        elif isinstance(formulas, dict):
            self.formulas.append(encode_formula_to_array(formulas));
        elif isinstance(formulas, np.ndarray):
            self.formulas.append(formulas);
        else:
            raise TypeError('Wrong type argument (formulas) for FormulaVectors initialization! str, formula, list of (formula, formulavector or str) supported only!');
        
        if isinstance(scores, str):
            scores=scores.split(',');
            for score in scores:
                self.scores.append(float(score));
        elif  isinstance(scores, float):
            self.scores.append(scores);
        elif  isinstance(scores, int):
            self.scores.append(float(scores));
        elif isinstance(scores, list):
            for score in scores:
                if isinstance(score,float):
                    self.scores.append(score);                    
                elif isinstance(score,dict):
                    self.scores.append(score);                    
                elif isinstance(score,int):
                    self.scores.append(float(score));                    
                elif isinstance(score,str):
                    score=score.split(',');
                    for s in score:
                        self.scores.append(float(s));
                else:
                    raise TypeError('Wrong type argument (scores) for FormulaScorer initialization! str, float, int, list of (float, int or str) or dictionary supported only!');
        else:
            raise TypeError('Wrong type argument (scores) for FormulaScorer initialization! str, float, int, list of (float, int or str) supported only!');
        if len(self.scores)!=len(self.formulas):
            raise TypeError('Number of formulas and number of scores supplied do not match!');
        
        if not self.vector_form:
            for i in range(len(self.formulas)):
                self.formulas[i]=decode_formula_from_array(self.formulas[i]);
            
    
    
    
    def __process_molecular_candidate_record_vector(self, source_data_base, record):
        testformula=record['FormulaVector'];
        found=False;
        for index in range(len(self.formulas)):
            formula=self.formulas[index];
            if np.array_equal(formula, testformula):
                    found=True;
                    if isinstance(self.scores[index],float):
                        record['Scores']['Formula']=self.scores[index];
                    elif isinstance(self.scores[index],dict):
                        for key in self.scores[index].keys():
                            record['Scores'][key]=self.scores[index][key];
                    break;
        if not found:
            record['Scores']['Formula']=self.unknown_score;


    def __process_molecular_candidate_record_formula(self, source_data_base, record):
        testformula=record['Formula'];
        found=False;
        for index in range(len(self.formulas)):
            formula=self.formulas[index];
            if len(testformula.keys())==len(formula.keys()):
                accepted=True;
                for atom in formula.keys():
                    if atom in testformula:
                        if formula[atom]!=testformula[atom]:
                            accepted=False;
                            break;
                    else:
                        accepted=False;
                        break;
                if accepted:
                    found=True;
                    if isinstance(self.scores[index],float):
                        record['Scores']['Formula']=self.scores[index];
                    elif isinstance(self.scores[index],dict):
                        for key in self.scores[index].keys():
                            record['Scores'][key]=self.scores[index][key];
                    break;
        if not found:
            record['Scores']['Formula']=self.unknown_score;
        
registered_scorers[FormulaScorer.__name__]=FormulaScorer;

if __name__=='__main__':
    
    test=FormulaScorer(use_vector_form=False);
    test.setup_scorer('C2H5OH,CH4,PO4','0.3,0.1,0.5',1.0);
    
    record={'Formula':parse_formula('CH4'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'Formula':parse_formula('PO4'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'Formula':parse_formula('C2H5OH'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'Formula':parse_formula('CH3'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    test.setup_scorer('C2H5OH',0.35,1.0);

    record={'Formula':parse_formula('C2H6O'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'Formula':parse_formula('PO4'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    test.setup_scorer([parse_formula('C2H5OH'),'PO4,CH4'],[0.35,'0.7',15],0.0);

    record={'Formula':parse_formula('C2H6O'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'Formula':parse_formula('PO4'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'Formula':parse_formula('CH4'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);

    record={'Formula':parse_formula('O4'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
        
    test.setup_scorer([parse_formula('C2H5OH'),'PO4,CH4'],[0.35,'0.7',{'Float1':15.0,'Dedo':3.0,'Formula':35.0}],0.0);
    

    record={'Formula':parse_formula('C2H6O'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'Formula':parse_formula('PO4'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'Formula':parse_formula('CH4'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);

    record={'Formula':parse_formula('O4'),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    print('FormulaVector');
    test=FormulaScorer();
    test.setup_scorer('C2H5OH,CH4,PO4','0.3,0.1,0.5',1.0);
    
    record={'FormulaVector':encode_formula_to_array(parse_formula('CH4')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'FormulaVector':encode_formula_to_array(parse_formula('PO4')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'FormulaVector':encode_formula_to_array(parse_formula('C2H5OH')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'FormulaVector':encode_formula_to_array(parse_formula('CH3')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    test.setup_scorer(encode_formula_to_array(parse_formula('C2H5OH')),0.35,1.0);

    record={'FormulaVector':encode_formula_to_array(parse_formula('C2H6O')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'FormulaVector':encode_formula_to_array(parse_formula('PO4')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    test.setup_scorer([encode_formula_to_array(parse_formula('C2H5OH')),'PO4,CH4'],[0.35,'0.7',15],0.0);
    print(test)

    record={'FormulaVector':encode_formula_to_array(parse_formula('C2H6O')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'FormulaVector':encode_formula_to_array(parse_formula('PO4')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'FormulaVector':encode_formula_to_array(parse_formula('CH4')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);

    record={'FormulaVector':encode_formula_to_array(parse_formula('O4')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
        
    test.setup_scorer([parse_formula('C2H5OH'),'PO4,CH4'],[0.35,'0.7',{'Float1':15.0,'Dedo':3.0,'Formula':35.0}],0.0);

    record={'FormulaVector':encode_formula_to_array(parse_formula('C2H6O')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'FormulaVector':encode_formula_to_array(parse_formula('PO4')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'FormulaVector':encode_formula_to_array(parse_formula('CH4')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);

    record={'FormulaVector':encode_formula_to_array(parse_formula('O4')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    print(test)
    