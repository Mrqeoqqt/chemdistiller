# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:28:27 2017

@author: Dr. Ivan Laponogov
"""
import os;


def get_files_list(path = os.getcwd(), recursive = True, extension = '.*'):
    result = [];
    if recursive:
        result = [os.path.join(dp, f) for dp, dn, filenames in os.walk(path) for f in filenames if (os.path.splitext(f)[1] == extension) or extension =='.*']
    else:
        result = [os.path.join(path, f) for f in os.listdir(path) if (os.path.isfile(os.path.join(path, f))) and (os.path.splitext(f)[1] == extension or extension =='.*')];
    
    return result;

if __name__=='__main__':
    l = get_files_list(recursive = True, extension = '.pyc')
    for s in l:
        print(s)