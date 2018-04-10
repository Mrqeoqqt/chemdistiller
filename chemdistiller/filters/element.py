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


#TODO: Change to string formula representation

class ElementCompositionFilter:
    name='ElementCompositionFilter';
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
                if s[0].lower().startswith('element_filter'):
                    self.element_filter=string_to_numpy_byte_array(s[1]);
                elif s[0].lower().startswith('compulsory_element_filter'):
                    self.compulsory_element_filter=string_to_numpy_byte_array(s[1]);

            elif s.lower().startswith('end'):
                return;
        
    def _pipe_to_textfile(self, fout, indent=''):                
        fout.write('%s%s=%s\n'%(indent,'element_filter', numpy_array_to_string(self.element_filter)));
        fout.write('%s%s=%s\n'%(indent,'compulsory_element_filter', numpy_array_to_string(self.compulsory_element_filter)));
                


    def __init__(self, compulsory_elements='', elements=''):
        
        if isinstance(elements, str):
            self.element_filter=formula_to_element_vector(parse_formula(elements));
        elif isinstance(elements, dict):
            self.element_filter=formula_to_element_vector(elements);
        elif isinstance(elements, list):
            self.element_filter=list_to_element_vector(elements);
        elif isinstance(elements, set):
            self.element_filter=set_to_element_vector(elements);
        else:
            raise TypeError('Wrong type argument for ElementCompositionFilter initialization! str, dict, set and list supported only!');

        if isinstance(compulsory_elements, str):
            self.compulsory_element_filter=formula_to_element_vector(parse_formula(compulsory_elements));
        elif isinstance(compulsory_elements, dict):
            self.compulsory_element_filter=formula_to_element_vector(compulsory_elements);
        elif isinstance(compulsory_elements, list):
            self.compulsory_element_filter=list_to_element_vector(compulsory_elements);
        elif isinstance(compulsory_elements, set):
            self.compulsory_element_filter=set_to_element_vector(compulsory_elements);
        else:
            raise TypeError('Wrong type argument for ElementCompositionFilter initialization! str, dict, set and list supported only!');
        self.element_filter=np.bitwise_or(self.element_filter,self.compulsory_element_filter);
            
        self.required_fields=['ElementVector'];
    
    def rejected(self, testrecord):
        return not_within_elements_binary(testrecord['ElementVector'], self.element_filter) or\
        not_within_elements_binary(self.compulsory_element_filter, testrecord['ElementVector']);
        
    def __repr__(self):
        s1=''.join(element_vector_to_list(self.compulsory_element_filter));
        s2=''.join(element_vector_to_list(self.element_filter));
        return 'ElementCompositionFilter(compulsory="%s", allowed="%s")'%(s1,s2);


registered_filters[ElementCompositionFilter.__name__]=ElementCompositionFilter;
#full_periodic_table_filter=np.full((12,),255,dtype=np.uint8);

CHNOPS_filter=ElementCompositionFilter('','CHNOPS');

CHNOPS_halogens_filter=ElementCompositionFilter('','CHNOPSFClBrIAt');

if __name__=='__main__':
    
    test=ElementCompositionFilter('Cl','CHON');
    
    print(test.rejected({'ElementVector':formula_to_element_vector(parse_formula('C2H5'))}));
    print(test.rejected({'ElementVector':formula_to_element_vector(parse_formula('C2H5Cl'))}));
    print(test.rejected({'ElementVector':formula_to_element_vector(parse_formula('C2H5S'))}));
    print(test.rejected({'ElementVector':formula_to_element_vector(parse_formula('C2H5SCl'))}));

    test=ElementCompositionFilter('Cl','CHONCl');
    
    print(test.rejected({'ElementVector':formula_to_element_vector(parse_formula('C2H5'))}));
    print(test.rejected({'ElementVector':formula_to_element_vector(parse_formula('C2H5Cl'))}));
    print(test.rejected({'ElementVector':formula_to_element_vector(parse_formula('C2H5S'))}));
    print(test.rejected({'ElementVector':formula_to_element_vector(parse_formula('C2H5SCl'))}));
    
    test=ElementCompositionFilter('','CHONCl');
    print(test)
    