# -*- coding: utf-8 -*-
"""
Created on Thur April 14 2018
@author: Mingyi Xue

This module defines the structure of chemdistiller input files,
and read and write method.
"""
import os
import sys
from cd_mgf_keywords import CD_DEFAULT_SUB_PARAMS
from cd_mgf_keywords import CD_DEFAULT_PARAMS



if __name__=='__main__':
    if sys.byteorder!='little':
        #print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);


class CD:
    """
    This class is a structure for chemdistiller input.
    Receive parameters to document or generate a chemdistiller input from a MGF file.
    """
    def __init__(self,params=CD_DEFAULT_PARAMS,sub_params=CD_DEFAULT_SUB_PARAMS):
        self.parameters = params
        self.sub_parameters=sub_params

    def get_params(self):
        return self.parameters

    def write_cd(self,file):
        # this function is defined for conversion
        # from  mgf files to cd_files according to parameters
        with open(file,'w') as f:
            for (k,v) in self.parameters.items():
                if k=='db_molecule_name':
                    # if molecule already has a name, use it
                    # else assign the file name for it
                    # name will be used in .html report generation
                    f.write(k + '=')
                    if v=='':
                        n,e=os.path.splitext(os.path.basename(file))
                        f.write(n+'\n')
                    else:
                        f.write(v+'\n')
                    continue
                if k=='peaks':
                    continue
                if k=='ion_type':
                    continue
                f.write(k+'='+v+'\n')
            f.write('peaks'+'\n')
            f.write('\t1,'+self.parameters['peaks']+',100.0'+'\n')
            f.write('\tion_type='+self.parameters['ion_type']+'\n')

            f.write('\t\tspectrum\n')
            for (k,v) in self.sub_parameters.items():
                if k=='peaks':
                    continue
                f.write('\t\t\t'+k + '=' + v + '\n')

            f.write('\t\t\tpeaks\n')
            for i in range(len(self.sub_parameters['peaks'])):
                # deal with the case level=1
                ls=self.sub_parameters['peaks'][i].split()
                if len(ls)>=2:
                    f.write('\t\t\t\t' + str(i + 1) + ',')
                    f.write(ls[0]+','+ls[1]+'\n')
                    #if length>2:
                        # for peak info more than m/z and intensity
                    #    pass
            f.write('\t\t\tend\n')
            f.write('\t\tend\n')
            f.write('end\n')

    def gen_from_mgf(self, mgf):
        # this function is defined for parameters to conevert
        # from cd_files to mgf
        # a reverse process of chemdistiller.msspectra.spectrum._pipe_from_textfile
        params=mgf.get_para()
        for (k,v) in params.items():
            if k=='peaks':
                self.sub_parameters['peaks']=v
            elif k=='PEPMASS':
                # set precursor and precursor_mz
                self.parameters['peaks']=v
                self.sub_parameters['precursor_mz']=v
            else:
                self.parameters[str.lower(k)] = v

        if int(self.parameters['level'])!=1:
            # check 'level'
            # set level always to 1
            self.parameters['level']='1'

        peaks = self.sub_parameters['peaks']


        p1=peaks[0].split()
        p2=peaks[-1].split()
        if float(p1[0]) > float(p2[0]):
            main_peak=p1
        else:
            main_peak = p2
        # find the largest m/z peak

        if self.parameters['peaks']=='':
            # if PEPMASS not exist,
            # then find the largiest m/z and set as precursor
            self.parameters['peaks'] = main_peak[0]
            self.sub_parameters['precursor_mz'] = main_peak[0]

        if self.parameters['ion_type']=='':
            if len(main_peak) == 3:
                self.parameters['ion_type']=main_peak[2]
                self.sub_parameters['precursor_ion']=main_peak[2]
                if '-' in self.parameters['peaks'][2]:
                    self.parameters['mode'] = '-1'
                    self.sub_parameters['mode'] = '-1'
                else:
                    self.parameters['mode'] = '1'
                    self.sub_parameters['mode'] = '1'
            else:
                self.parameters['ion_type'] = r'[M+H]+'
                self.sub_parameters['precursor_ion'] = r'[M+H]+'
                self.parameters['mode'] = '1'
                self.sub_parameters['mode'] = '1'

        self.sub_parameters['charge']=self.parameters['charge']
        self.sub_parameters['exactmass']=self.parameters['exactmass']
        self.sub_parameters['formula']=self.sub_parameters['formula']
        self.sub_parameters['inchi'] = self.sub_parameters['inchi']










