# -*- coding: utf-8 -*-
"""
Created on Thur April 12 2018
@author: Mingyi Xue

This module defines the structure of mgf input files,
and read and write method.
"""

from cd_mgf_keywords import MGF_DEFAULT_PARAMS


class MGF:
    """
    This module abstracts the structure of a mgf file.
    """
    def __init__(self, params=MGF_DEFAULT_PARAMS):
        self.parameters = params

    def get_para(self):
        return self.parameters

    def write_mgf(self,file):
        # this function is reserved for conversion
        # from cd_files to mgf files according to parameters
        pass

    def gen_from_cd(self,cd):
        # this function is reserved for parameters to convert
        # from cd_files to mgf
        pass


class MGFfile:
    """
    read chunks from a .mgf file, containing one or more begin ions...end ions chunks
    create an object: mgf=MGFfile(fname)
    generate an array of MGF chunks: mgf.gen_mgfs()
    get MGF chunk array: mgf_list=get_mgfs()
    """
    def __init__(self,fname):
        self.MGFs=[]
        # a list of MGF
        self.fname=fname

    def gen_mgfs(self):
        with open(self.fname) as f:
            while True:
                params={}
                params['peaks']=[]
                s = f.readline()
                if s=='':
                    break
                s= s.lstrip().upper()
                if s.startswith('BEGIN IONS'):
                    while True:
                        s=f.readline().rstrip('\n').strip()
                        # strip all space characters in a line
                        if s.startswith('END IONS'):
                            mgf = MGF(params)
                            self.MGFs.append(mgf)
                            break
                        if '=' in s:
                            s=s.split('=', 1)
                            # divide by '=' only once
                            params[s[0]] = s[1]
                        else:
                            params['peaks'].append(s)
                            # divide by space characters

    def get_mgfs(self):
        if len(self.MGFs)==0:
            self.gen_mgfs()
        return self.MGFs