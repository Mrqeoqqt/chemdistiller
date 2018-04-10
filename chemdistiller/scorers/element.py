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

from chemdistiller.utils.periodictable import parse_formula, elements_dict, formula_to_element_vector, elements_list;
from chemdistiller.settings import test_as_single_process, registered_scorers, registered_filters;
from chemdistiller.utils.utils import numpy_array_to_string, string_to_numpy_float32_array;

import numpy as np;


class ElementScorer:
    name='ElementScorer';

    def __init__(self):
        self.required_fields=['ElementVector'];
        self.element_scores={};
        
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
                if s[0].lower().startswith('element_scores'):
                    s=s[1].split(':');
                    self.element_scores[s[0]]=float(s[1]);
            elif s.lower().startswith('end'):
                return;
        
    def _pipe_to_textfile(self, fout, indent=''):                
        for element in self.element_scores.keys():
            fout.write('%s%s=%s:%s\n'%(indent,'element_scores', element, self.element_scores[element]));
    
    
    
    def __repr__(self):
        l=[];
        for element in self.element_scores.keys():
            l.append('"%s"=%.3f'%(elements_list[element], self.element_scores[element]))
        s=', '.join(l);
        return 'ElementScorer(%s)'%s;
        
    def close(self):
        return 0;

    def configure_scorer(self, peak, annotation): 
        #This scorer is annotation-independent
        return 0;


    def setup_scorer(self, element_scores):
        self.element_scores={};
        for key in element_scores.keys():
            self.element_scores[elements_dict[key]]=element_scores[key];
        
    def process_molecular_candidate_record(self, source_data_base, record):
        testvector=np.unpackbits(record['ElementVector']);
        score=1.0;
        for i in self.element_scores.keys():
            if testvector[i]>0:
                score=score*self.element_scores[i];
        record['Scores']['Elements']=score;

registered_scorers[ElementScorer.__name__]=ElementScorer;


if __name__=='__main__':
    
    test=ElementScorer();
    test.setup_scorer({'C':0.8,'O':0.9,'Si':0.1});
    
    record={'ElementVector':formula_to_element_vector(parse_formula('CH4')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'ElementVector':formula_to_element_vector(parse_formula('CH4Si')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'ElementVector':formula_to_element_vector(parse_formula('CH4CO')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    record={'ElementVector':formula_to_element_vector(parse_formula('Ar')),'Scores':{}};
    test.process_molecular_candidate_record(None,record);
    print(record);
    
    print(test)
    
    
    