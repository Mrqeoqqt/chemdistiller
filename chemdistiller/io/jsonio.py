# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 15:12:42 2016

@author: Dr. Ivan Laponogov
"""

import json;

import numpy as np;

import sys;
import os;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);

from chemdistiller.msspectra.peak import MSPeak, MSBasicPeak;
from chemdistiller.msspectra.spectrum import MSSpectrum;
from chemdistiller.annotations.peakannotation import MSPeakAnnotation;
from chemdistiller.msspectra.manager import SpectralManager;
from chemdistiller.utils.base64 import encode_to_base64;
from chemdistiller.filters.element import ElementCompositionFilter;
from chemdistiller.filters.inchi import InChIFilter;
from chemdistiller.filters.formula import FormulasFilter;
from chemdistiller.utils.periodictable import element_vector_to_list, Formula, parse_formula;

class MS_Data_Encoder(json.JSONEncoder):
    def default(self, obj):
       #print(type(obj));
       
       if isinstance(obj, Formula):    
           return obj.formula_to_string();
           
       elif isinstance(obj, MSBasicPeak):
           return {"mz":obj.mz, 'intensity':obj.intensity};
           
       elif isinstance(obj, MSPeak):
            result={"mz":obj.mz, 'intensity':obj.intensity, 'number':obj.number};
            if hasattr(obj,'ppm'):
                if not (obj.ppm is None):
                    result['ppm']=obj.ppm;
                    
            if hasattr(obj,'parameters'):
                if obj.parameters:
                    result['parameters']=obj.parameters;

            if hasattr(obj,'annotations'):
                if obj.annotations:
                    result['annotations']=obj.annotations;
            
            if hasattr(obj,'ms_spectra'):
                if obj.ms_spectra:
                    result['ms_spectra']=obj.ms_spectra;


            if hasattr(obj,'subpeak_list'):
                if not (obj.subpeak_list is None):
                    if obj.subpeak_list:
                        result['subpeak_list']=obj.subpeak_list;

            if hasattr(obj,'fpt_sparse_vector'):
                if not (obj.fpt_sparse_vector is None):
                    if obj.fpt_sparse_vector:
                        result['fpt_sparse_vector']=obj.fpt_sparse_vector;

            if hasattr(obj,'merged_annotations'):
                if not (obj.merged_annotations is None):
                    if obj.merged_annotations:
                        result['merged_annotations']=obj.merged_annotations;
            return result;
            
       elif isinstance(obj, MSSpectrum):
            return {"peaks":obj.peaks,"parameters": obj.parameters}

       elif isinstance(obj, MSPeakAnnotation):
           result={'adduct':obj.adduct.definition, 'isotope':obj.isotope, 'isotope_extra_mass':obj.isotope_extra_mass};

           if hasattr(obj,'predicted_fpt'):
               result['predicted_fpt']=encode_to_base64(obj.predicted_fpt);
           if hasattr(obj,'fpt_mask'):
               result['fpt_mask']=encode_to_base64(obj.fpt_mask);
           if not (obj.mol_candidates is None):
               result['mol_candidates']=obj.mol_candidates.mol_list;
           if obj.scores:
               result['scores']=obj.scores;
           if not (obj.formula_scorer is None):
               fscorer=obj.formula_scorer;
               fsc={'unknown_score':fscorer.unknown_score,'formula_scores':[]};
               lst=fsc['formula_scores'];
               for i in range(len(fscorer.formulas)):
                   lst.append({'formula':fscorer.formulas[i],'scores':fscorer.scores[i]});
               result['formula_scorer']=fsc;
           if obj.filters:
               result['filters']=obj.filters;
           
           return result;
           
       elif isinstance(obj, SpectralManager):
           return {"spectra":obj.ms_spectra};
       
       elif isinstance(obj, ElementCompositionFilter):
           return {"elements_filter":{"allowed_elements":element_vector_to_list(obj.element_filter),"compulsory_elements":element_vector_to_list(obj.compulsory_element_filter)}};
           
       elif isinstance(obj, InChIFilter):
           return {"InChI_filter":{"reference_InChI":obj.ref_inchi,'use_short_InChI':obj.use_short_inchi,'match_type':obj.match_type}};

       elif isinstance(obj, FormulasFilter):
           return {"formulas_filter": obj.filter};
           
       elif isinstance(obj, np.ndarray):    
           return encode_to_base64(obj);
           
       elif isinstance(obj, bytearray):    
           return encode_to_base64(obj);
           
       elif isinstance(obj, np.float32):    
           return float(obj);
           
       elif isinstance(obj, np.uint16):    
           return int(obj);

       elif isinstance(obj, np.uint8):    
           return int(obj);
           
       
       return json.JSONEncoder.default(self, obj)


def spectra_to_json(fname, spectral_list):
    dirname, filename=os.path.split(fname);
    if not os.path.exists(dirname):
            os.makedirs(dirname);
    with open(fname,'w') as fout:
        json.dump(spectral_list, fout, cls=MS_Data_Encoder, sort_keys=True, indent=4, separators=(',', ': '))
        

if __name__=='__main__':

    spectra_to_json(os.path.join(chemdistiller_path,'testlogs/t1.jsn'),[parse_formula('C2H5OH'),25]);

    