# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:16:17 2017

@author: ilaponog
"""
import os;
import csv;
import sys;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."));
    sys.path.append(chemdistiller_path);

from chemdistiller.filters.formula import FormulasFilter;
from chemdistiller.filters.element import ElementCompositionFilter;
from chemdistiller.utils.periodictable import parse_formula, formula_to_element_vector, set_to_element_vector, list_to_element_vector,\
                                               not_within_elements_binary, element_vector_to_list, encode_formula_to_array;


def process_file(fname, correct_key, best_results, worst_results, correct_formula):
    min_correct=-1;
    max_correct=-1;
    
    correct_elements = ElementCompositionFilter(correct_formula, correct_formula)
    correct_formula=FormulasFilter(correct_formula);   
    
    with open(fname, 'rb') as finp:
        results = list(csv.reader(finp));
        del results[0];

        for i in reversed(range(1,len(results))):
            if results[i][9]==results[i-1][9]:
                del results[i];

            elif correct_elements.rejected({'ElementVector':formula_to_element_vector(parse_formula(results[i][6]))}):
                #print('removed wrong elements %s %s'%(results[i][6], correct_formula))
                del results[i];
                
            elif correct_formula.rejected({'FormulaVector':encode_formula_to_array(parse_formula(results[i][6]))}):
                #print('removed wrong formulas %s %s'%(results[i][6], correct_formula))
                del results[i];
                                    
                 
                
                
        for i in range(len(results)):
            if correct_key==results[i][9]:
                if min_correct==-1:
                    min_correct=i;
                max_correct=i;
        if min_correct>-1:
            for i in range(len(best_results)):
                if min_correct<=i:
                    best_results[i]+=1;
                if max_correct<=i:
                    worst_results[i]+=1;

def results_to_csv(fname, best_results, worst_results, header, count):
    with open(fname,'w') as fout:
        fout.write('%s\n'%header);
        fout.write('Correct in the first:\tBest, %\tWorst, %\n');
        for i in range(len(best_results)):
            fout.write('%s\t%s\t%s\n'%(i+1,float(best_results[i])/count*100, float(worst_results[i])/count*100));
            

inpath='G:/MetFrag/smallDB0/positive';

header='MetFrag retrieval, Batch 0, small DB (CheBI,HMDB, MassBank, NIST14), positive mode, formula';

onlyfiles = [ff for ff in os.listdir(inpath) if os.path.isfile(os.path.join(inpath, ff))];

max_results_per_query=100;

#print(onlyfiles);

best_results=[0]*max_results_per_query;
worst_results=[0]*max_results_per_query;

count=0;
for fname in onlyfiles:
    fname,fext= os.path.splitext(fname);
    if fext=='.inchi':
        print(fname)
        with open(os.path.join(inpath, '%s.inchi'%fname),'r') as finp:
            for s in finp:
                correct_formula=s.rstrip('\n').split('\t')[0].split('/')[0]; 
                correct_key=s.rstrip('\n').split('\t')[1].split('-')[0];
                
        ff=os.path.join(inpath, 'results','%s.csv'%fname);
        if os.path.isfile(ff):
            count+=1;
            process_file(ff, correct_key, best_results, worst_results, correct_formula);
            
results_to_csv(os.path.join(inpath,'results/statistics_formula.csv'), best_results, worst_results, header, count);
            
            
        
