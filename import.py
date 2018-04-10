# -*- coding: utf-8 -*-
"""
==============================================================================
        ChemDistiller Data Import Script
==============================================================================
    This script is designed for conversion of different formats of MS2 spectra 
into internal data format for ChemDistiller.


run python.exe import.py -h for full help about supported options


Copyright:                      Imperial College London
Project:                        MetaSpace
Funding:                        Horizon 2020
License:                        BSD
Chief project investigator:     Dr. Kirill Veselkov
Lead developer:                 Dr. Ivan Laponogov

"""


import sys;
import os;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.dirname(os.path.realpath(__file__));
    sys.path.append(chemdistiller_path);

from chemdistiller.msspectra.manager import SpectralManager;

from chemdistiller.utils.cmdline import OptionsHolder;
from chemdistiller.procconfig import ImportSet_options;

from chemdistiller.utils.filelist import get_files_list;

from chemdistiller.io.mgf import import_mgf;
from chemdistiller.io.massbank import import_massbank;
from chemdistiller.io.mzml import import_mzML;

from chemdistiller.msspectra.preprocessing import purge_ms1, merge_ms1;


'''
Options to implement:
input file name
include formats - mgf, mzML, txt*
input folder - bulk import
recursive?
input metaspace annotation
output folder if input file or folder
ms1 mz_range
ms2 mz_range
ms1 ppm
ms2 ppm
ms1 background level
ms2 background level
spec ranges:
*-10,20,35-45,80-*
remove ms1 with no ms2?
keep largest peak only?
merge ms1-s?
print contents

'''

def postprocess_spectra(spectral_manager, parameters):
    if parameters['purge_ms1'] != 'no':
        purge_ms1(spectral_manager);
    if parameters['merge_ms1'] != 'no':
        merge_ms1(spectral_manager, parameters);


def import_spectra_from_file(fname, spectral_manager, parameters):
    fname = os.path.abspath(fname);
    fn, fex = os.path.splitext(fname);
    formats = parameters['import_formats'].split(',');
    use_formats = set();
    for f in formats:
        if '*' in f:
            use_formats.add('.mgf');
            use_formats.add('.mzML');
            use_formats.add('.txt');
        elif 'mgf' in f:
            use_formats.add('.mgf');
        elif ('mzML' in f) or ('mzml' in f):
            use_formats.add('.mzML');
            use_formats.add('.mzml');
        elif ('txt' in f):
            use_formats.add('.txt');
    fex = fex.lower()
    if fex in use_formats:
        print('Importing %s....'%fname);
        if fex == '.mgf':
            import_mgf(fname, spectral_manager, parameters)
        elif fex == '.txt':
            import_massbank(fname, spectral_manager, parameters)
        else:
            import_mzML(fname, spectral_manager, parameters)


if __name__=='__main__':
    
    settings=OptionsHolder(__doc__, ImportSet_options);
    settings.description='Tandem MS import';
    print(settings.program_description);
    settings.parse_command_line_args();
    print(settings.format_parameters()); #Print current parameters
    print('\nStarting.....')

    if settings.parameters['no_output']=='no':
        if not os.path.isdir(settings.parameters['output']):
            os.makedirs(settings.parameters['output'])
        else:
            print('Error! Output path already exists! Please provide a new path to avoid overwriting of existing files!')
            #sys.exit(-1);


    with SpectralManager() as spectral_manager:
        if os.path.isfile(settings.parameters['inputdatapath']):
            import_spectra_from_file(settings.parameters['inputdatapath'], spectral_manager, settings.parameters);
        else:
            flist = get_files_list(settings.parameters['inputdatapath'], settings.parameters['recursive'] =='yes');
            for fname in flist:
                import_spectra_from_file(fname, spectral_manager, settings.parameters);
        
        if settings.parameters['no_output']=='no':
            postprocess_spectra(spectral_manager, settings.parameters);
            spectral_manager.export_textfile_spectra_to_folder(settings.parameters['output']);

    
    print('\nFinished.');
    print(settings.description_epilog);
    
    
    

