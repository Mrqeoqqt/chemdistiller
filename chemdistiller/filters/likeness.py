# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 14:21:37 2016

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

from chemdistiller.utils.periodictable import parse_formula, formula_to_element_vector, set_to_element_vector, list_to_element_vector,\
                                               not_within_elements_binary, element_vector_to_list;
from chemdistiller.utils.utils import numpy_array_to_string, string_to_numpy_byte_array;
from chemdistiller.settings import test_as_single_process, registered_scorers, registered_filters;

                                    
import numpy as np;                                               

class LikenessFilter:
    name='LikenessFilter';
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
                if s[0].lower().startswith('likeness'):
                    self.likeness.append(s[1]);
                if s[0].lower().startswith('unlikeness'):
                    self.unlikeness.append(s[1]);
                
            elif s.lower().startswith('end'):
                return;
        
    def _pipe_to_textfile(self, fout, indent=''):                
        for i in self.likeness:
            fout.write('%s%s=%s\n'%(indent,'likeness', i));
        for i in self.unlikeness:
            fout.write('%s%s=%s\n'%(indent,'unlikeness', i));
                


    def __init__(self, likeness=None, unlikeness=None):
        if likeness is None:
            self.likeness=[];
        else:
            self.likeness=likeness;
        if unlikeness is None:
            self.unlikeness=[];
        else:
            self.unlikeness=unlikeness;
    
        if (not self.likeness) and (not self.unlikeness):
            raise TypeError('Both likeness and unlikeness are empty!');
            
        self.required_fields=self.likeness+self.unlikeness;
    
    def rejected(self, testrecord):
        
        if self.likeness:
            result=True;
            #print('Nonempty likeness')
            for i in self.likeness:
                if i in testrecord:
                    if testrecord[i]>0:
                        result=False;
                        break;
        else:
            #print('Any likeness');
            result=False;
        
        if result==False:
            for i in self.unlikeness:
                if i in testrecord:
                    if testrecord[i]>0:
                        return True;
        return result;
        
    def __repr__(self):
        if self.likeness:        
            s1=', '.join(self.likeness);
        else:
            s1='Any';
        if self.unlikeness:
            s2=', '.join(self.unlikeness);
        else:
            s2='None';
        return 'LikenessFilter(Likeness: %s; Unlikeness: %s)'%(s1,s2);


registered_filters[LikenessFilter.__name__]=LikenessFilter;



if __name__=='__main__':
    
    test=LikenessFilter(['HMDB','LipidMaps'],['Zinc']);
    
    print(test.rejected({'HMDB':1}));
    print(test.rejected({'LipidMaps':1, 'HMDB':-1}));
    print(test.rejected({'LipidMaps':-1, 'HMDB':-1}));
    print(test.rejected({'HMDB':-1}));
    print(test.rejected({'Zinc':1}));
    print(test.rejected({'Zinc':-1}));

    print(test.rejected({'HMDB':1, 'Zinc':-1}));
    print(test.rejected({'HMDB':1, 'Zinc':1}));
    print(test.rejected({'HMDB':-1, 'Zinc':-1}));
    
    print(test)
    test=LikenessFilter([],['Zinc']);
    print(test.rejected({'LipidMaps':1, 'HMDB':-1, 'Zinc':1}));
    print(test.rejected({'LipidMaps':-1, 'HMDB':-1, 'Zinc':-1}));
    print(test)
    test=LikenessFilter(['HMDB'],[]);
    print(test.rejected({'LipidMaps':1, 'HMDB':1, 'Zinc':1}));
    print(test.rejected({'LipidMaps':-1, 'HMDB':-1, 'Zinc':-1}));
    print(test)
    