# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:27:43 2017

@author: ilaponog
"""
import numpy as np;

def numpy_array_to_string(array):
    s='';
    if len(array)>0:
        s='%s'%array[0];
        for i in range(1, len(array)):
            s='%s,%s'%(s,array[i]);
    return s;
    
def string_to_numpy_byte_array(s):
    s=s.split(',');
    return np.array(s,dtype=np.uint8);

def string_to_numpy_uint16_array(s):
    s=s.split(',');
    return np.array(s,dtype=np.uint16);
    
def string_to_numpy_float32_array(s):
    s=s.split(',');
    return np.array(s,dtype=np.float32);    
    
def string_to_float_list(s):
    s=s.split(',');
    result=[];
    for ss in s:
        result.append(float(ss));
    return result;