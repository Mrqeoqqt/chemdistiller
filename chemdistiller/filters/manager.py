# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 18:39:13 2017

@author: ilaponog
"""

import sys;
import os;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);


from chemdistiller.filters.formula import FormulasFilter;
from chemdistiller.filters.inchi import InChIFilter;
from chemdistiller.filters.element import ElementCompositionFilter;


def get_filter_by_name(filter_name):

    result=None;

    if filter_name==FormulasFilter.name:
        result=FormulasFilter();
    
    elif filter_name==ElementCompositionFilter.name:
        result=ElementCompositionFilter();
   
    elif filter_name==InChIFilter.name:
        result=InChIFilter();

    return result;