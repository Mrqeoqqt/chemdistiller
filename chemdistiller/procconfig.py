# -*- coding: utf-8 -*-
"""
ChemDistiller

Command line arguments definition and defaults

Dr. Ivan Laponogov
 
"""


import os

from chemdistiller.utils.cmdline import Option, Value, AllValues;


ImportSet_options=[Option('-u,--units', help = 'Units to process data in.', values=[\
                            Value('ppm', help = 'Parts per million, ppm.', parameters=[\
                                Option('--ms1_mzrange', help = 'm/z limits of imported MS1 data, e.g. 20.0-50.0 or * for all.', values=['*', None], type=str, targets=['/ms1']),\
                                Option('--ms1_mzprecision', help = 'm/z tolerance for MS1 data.', values=[20, None], type=float, targets=['/ms1']),\
                                Option('--ms2_mzrange', help = 'm/z limits of imported MS2 data, e.g. 20.0-50.0 or * for all.', values=['*', None], type=str, targets=['/ms2']),\
                                Option('--ms2_mzprecision', help = 'm/z tolerance for MS2 data.', values=[20, None], type=float, targets=['/ms2'])\
                                ]),
                                
                            Value('Da', help = 'Daltons, Da', parameters=[\
                                Option('--ms1_mzrange', help = 'm/z limits of imported MS1 data, e.g. 20.0-50.0 or * for all.', values=['*', None], type=str, targets=['/ms1']),\
                                Option('--ms1_mzprecision', help = 'm/z tolerance for MS1 data.', values=[0.01, None], type=float, targets=['/ms1']),\
                                Option('--ms2_mzrange', help = 'm/z limits of imported MS2 data, e.g. 20.0-50.0 or * for all.', values=['*', None], type=str, targets=['/ms2']),\
                                Option('--ms2_mzprecision', help = 'm/z tolerance for MS2 data.', values=[0.01, None], type=float, targets=['/ms2'])\
                                ]),\
                            ]),\
                   Option('-b1,--ms1_background', help = 'Noise level for MS1.', values=[0.0, None], type=float, targets=['/ms1']),\
                   Option('-b2,--ms2_background', help = 'Noise level for MS2.', values=[0.0, None], type=float, targets=['/ms2']),\
                   Option('-r,--recursive', help = 'Recursively import files from folder.', values=['no', 'yes'], type=str),\
                   Option('-p,--purge_ms1', help = 'Remove MS1 peaks with no associated MS2 spectra.', values=['no', 'yes'], type=str),\
                   Option('-m,--merge_ms1', help = 'Merge MS1 spectra within individual input files (internally) or between all input files (externally).', values=['no', 'internally', 'externally'], type=str),\
                   Option('-l,--largest_only', help = 'If merging, keep only the largest MS1 peak with corresponding MS2.', values=['no', 'yes'], type=str),\
                   Option('-i,--import_formats', help = 'Import files from folder of the listed formats. * - All. txt - for MassBank, mgf - for GNPS, mzML - for general mzML format.', values=['*', 'txt', 'mzML', 'mgf'], type=str, is_list = True),\
                   Option('-P,--print', help = 'Print info from selected files', values=['no', 'yes'], type=str),\
                   Option('-n,--no_output', help = 'Suppress file conversion. Combine with -P option to only display file info-s.', values=['no', 'yes'], type=str),\
                   Option('-M,--metaspace', help = 'Try to import Metaspace annotation csv files', values=['no', 'yes'], type=str),\
                   Option('-s,--spectra_range', help = 'Limit imported spectra from mzML to this range. E.g. *-20, 36, 48-50, 90-*', values=['*', None], type=str, is_list = True),\
                   
                   Option('inputdatapath', help = 'Sets input data directory path or file.', values=[os.getcwd(), None], type=str, conditions=[('Path or file must exist!', lambda x: os.path.isdir(x.inputdatapath) or os.path.isfile(x.inputdatapath))]),\
                   Option('output', help = 'Output folder for converted spectra. This path should not exist to avoid overwriting other files and will be created by the import.py', values=[os.getcwd()+'/ImportedSpectra', None], type=str)\
             
                   ];


