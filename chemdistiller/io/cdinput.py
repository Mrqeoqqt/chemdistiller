# -*- coding: utf-8 -*-
"""
Created on Thur April 14 2018
@author: Mingyi Xue

This module defines the structure of chemdistiller input files,
and read and write method.
"""
import os
import sys
from cd_mgf_keywords import CD_DEFAULT_PARAMS
from cd_mgf_keywords import CD_DEFAULT_SUB_PARAMS
from cd_mgf_keywords import MAP_OF_MGF_CD


if __name__=='__main__':
    if sys.byteorder!='little':
        #print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);


class CD:
    def __init__(self,params=CD_DEFAULT_PARAMS,sub_params=CD_DEFAULT_SUB_PARAMS):
        self.parameters = params
        self.sub_parameters=sub_params

    def get_params(self):
        #if len(self.parameters)==0:
        #    self.make_pairs()
        return self.parameters

    def write_cd(self,file):
        # this function is defined for conversion
        # from  mgf files to cd_files according to parameters
        with open(file,'w') as f:
            for (k,v) in self.parameters.items():
                if k=='peaks':
                    continue
                f.write(k+'='+v+'\n')
            f.write('peaks'+'\n')
            if int(self.parameters['level'])<2:
                # for tandem mass spectra
                pass
            else:
                # for MS
                for i in range(len(self.parameters['peaks'])):
                    # deal with the case level=1
                    f.write('\t'+str(i+1)+',')
                    ls=self.parameters['peaks'][i].split()
                    length=len(ls)
                    if length>=2:
                        f.write(ls[0]+','+ls[1]+'\n')
                        if length>2:
                            # for peak info more than m/z and intensity
                            pass
                    elif length<2:
                        continue
            f.write('end\n')

    def gen_from_mgf(self, mgf):
        # this function is defined for parameters to conevert
        # from cd_files to mgf
        # a reverse process of chemdistiller.msspectra.spectrum._pipe_from_textfile
        params=mgf.get_para()

        for (k,v) in params.items():
            if k=='peaks':
                self.parameters['peaks']=v
            elif k in MAP_OF_MGF_CD.keys():
                self.parameters[MAP_OF_MGF_CD[k]] = v
            else:
                self.parameters[str.lower(k)] = v







if __name__=='__main__':
    params = {'db_molecule_name': '',
                      'exactmass': '',
                      'formula': '',
                      'fpt_0': '',
                      'fptcount': '',
                      'global_index': '',
                      'inchi': '',
                      'level': '',
                      'mode': '',
                      'collision_energy': '',
                      'peaks': []
                      }





