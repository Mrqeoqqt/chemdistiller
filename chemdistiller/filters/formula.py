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
from chemdistiller.settings import test_as_single_process, registered_scorers, registered_filters;
import numpy as np;


class FormulasFilter:
    name='FormulasFilter';
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
                    self.filter.append(encode_formula_to_array(parse_formula(s[1])));
            elif s.lower().startswith('dict_form'):
                self.vector_form=False;
            
            elif s.lower().startswith('end'):
                if not self.vector_form:
                    self.required_fields=['Formula'];
                    for i in range(len(self.filter)):
                        self.filter[i]=decode_formula_from_array(self.filter[i]);
                    self.rejected=self.__rejected_formula;

                return;
        
    def _pipe_to_textfile(self, fout, indent=''):
        #subindent='%s\t'%indent;
        if self.vector_form:
            fout.write('%s%s\n'%(indent,'vector_form'));
        else:
            fout.write('%s%s\n'%(indent,'dict_form'));
        
        for filtervalue in self.filter:
            if self.vector_form:
                fout.write('%s%s=%s\n'%(indent,'formula',decode_formula_from_array(filtervalue).formula_to_string()));
            else:
                fout.write('%s%s=%s\n'%(indent,'formula',filtervalue.formula_to_string()));


    def __init__(self, formulas=None, use_vector_form=True):
        
        self.filter=[];
        self.vector_form=use_vector_form;
        #self.supported_adducts=set();
        if self.vector_form:
            self.required_fields=['FormulaVector'];
        else:
            self.required_fields=['Formula'];
            
        

        if formulas is None:
            return
        if isinstance(formulas, str):
            formulas=formulas.split(',');
            for formula in formulas:
                self.filter.append(encode_formula_to_array(parse_formula(formula)));
        elif isinstance(formulas, list):
            for formula in formulas:
                if isinstance(formula,dict):
                    self.filter.append(encode_formula_to_array(formula));
                elif isinstance(formula, np.ndarray):
                    self.filter.append(formula);
                elif isinstance(formula, str):
                    fs=formula.split(',');
                    for f in fs:
                        self.filter.append(encode_formula_to_array(parse_formula(f)));
                else:
                    raise TypeError('Wrong type argument for FormulasFilter initialization! str, dict, list of (dict or str) supported only!');
            
        elif isinstance(formulas, dict):
            self.filter.append(encode_formula_to_array(formulas));
        elif isinstance(formulas, np.ndarray):
            self.filter.append(formulas);    
        else:
            raise TypeError('Wrong type argument for FormulasFilter initialization! str, dict, list of (dict or str) supported only!');

        if not self.vector_form:
            for i in range(len(self.filter)):
                self.filter[i]=decode_formula_from_array(self.filter[i]);

    def __repr__(self):
        flas=[];
        for filtervalue in self.filter:
            if self.vector_form:
                flas.append('"%s"'%decode_formula_from_array(filtervalue).formula_to_string());
            else:
                flas.append('"%s"'%filtervalue.formula_to_string());

        s=', '.join(flas);
        return "FormulasFilter(%s)"%s;
            
    def rejected(self, testrecord):
        if self.vector_form:
            return self.__rejected_vector(testrecord);
        else:
            return self.__rejected_formula(testrecord);
    
    def __rejected_vector(self, testrecord):
        testformulavector=testrecord['FormulaVector'];
        for formulavector in self.filter:
            if np.array_equal(testformulavector,formulavector):
                return False;
        return True;

    def __rejected_formula(self, testrecord):
        testformula=testrecord['Formula'];
        result=True;
        for formula in self.filter:
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
                    result=False;
                    break;
        return result;

registered_filters[FormulasFilter.__name__]=FormulasFilter;

if __name__=='__main__':
    test=FormulasFilter('C2H5OH', use_vector_form=False);
    print(test.filter);
    test=FormulasFilter('C2H5OH,CH4,PO4', use_vector_form=False);
    print(test.filter);
    test=FormulasFilter(parse_formula('C2H5OH'), use_vector_form=False);
    print(test.filter);
    test=FormulasFilter(['C2H5OH,PO4',parse_formula('CH4')], use_vector_form=False);
    print(test.filter);
    print(test.rejected({'Formula':parse_formula('CH4')})); #False
    print(test.rejected({'Formula':parse_formula('PO4')})); #False
    print(test.rejected({'Formula':parse_formula('H4C')})); #False
    print(test.rejected({'Formula':parse_formula('CH4N')})); #True
    print(test.rejected({'Formula':parse_formula('CH3')})); #True
    print(test.rejected({'Formula':parse_formula('C2H4')})); #True
    print(test.rejected({'Formula':parse_formula('C')})); #True
    print(test.rejected({'Formula':parse_formula('CO4')})); #True
    print('FormulasVector');
    test=FormulasFilter('C2H5OH');
    print(test.filter);
    test=FormulasFilter('C2H5OH,CH4,PO4');
    print(test.filter);
    test=FormulasFilter(encode_formula_to_array(parse_formula('C2H5OH')));
    print(test.filter);
    test=FormulasFilter(['C2H5OH,PO4',parse_formula('CH4')]);
    
    print(test.filter);
    print(test.rejected({'FormulaVector':encode_formula_to_array(parse_formula('CH4'))})); #False
    print(test.rejected({'FormulaVector':encode_formula_to_array(parse_formula('PO4'))})); #False
    print(test.rejected({'FormulaVector':encode_formula_to_array(parse_formula('H4C'))})); #False
    print(test.rejected({'FormulaVector':encode_formula_to_array(parse_formula('CH4N'))})); #True
    print(test.rejected({'FormulaVector':encode_formula_to_array(parse_formula('CH3'))})); #True
    print(test.rejected({'FormulaVector':encode_formula_to_array(parse_formula('C2H4'))})); #True
    print(test.rejected({'FormulaVector':encode_formula_to_array(parse_formula('C'))})); #True
    print(test.rejected({'FormulaVector':encode_formula_to_array(parse_formula('CO4'))})); #True
    print(test)
    