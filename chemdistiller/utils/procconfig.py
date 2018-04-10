# -*- coding: utf-8 -*-
"""

***********************************************
Defaults for different modules to be used for 
command line options generation
***********************************************

 
"""


import os

import sys;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);


from chemdistiller.utils.cmdline import Option, Value, AllValues;

ImportMZML_options=[Option('-t,--filetype', help = 'Input File Format.', values=[\
                        Value('HDI', help = 'HDI File Format', parameters=[\
                            Option('--fileext', help = 'Set input file extension.', values=['.txt', None], type=str, targets=['/filereadinfo']),\
                            Option('--headerlines', help = 'Set input file header line count.', values=[5, None], type=int, conditions=[('headerlines>=0', lambda x: x.header_lines>=0)], targets=['/filereadinfo']),\
                            Option('--mzline', help = 'Set input file mzline value.', values=[3, None], type=int, conditions=[('mzline>=0', lambda x: x.mzline>=0)], targets=['/filereadinfo']),\
                            Option('--delimiter', help = 'Set input file delimiter.', values=['\t', None], type=str, targets=['/filereadinfo'])
                            ])#,
                        #Value('MSI', help = 'MSI File Format', parameters=[\
                        #    Option('--fileext', help = 'Set input file extension', values=['.msi', None], type=str),\
                        #    Option('--headerlines', help = 'Set input file header line count', values=[15, None], type=int),\
                        #    Option('--mzline', help = 'Set input file mzline value', values=[1, None], type=int),\
                        #    Option('--delimiter', help = 'Set input file delimiter', values=[',', None], type=str)
                        #    ]),
                        ], type=str),\
                   
                   #Option('--units', help = 'Input M/Z units', values=['Da'], type=str),\
                   Option('datapath', help = 'Sets data directory path.', values=[os.getcwd(), None], type=str, conditions=[('Path must exist!', lambda x: os.path.isdir(x.datapath))]),\
                   Option('h5dbname', help = 'Output HDF5 file name.', values=['', None], type=str)\
                   ];


'''

PeakAlign_options=[\

                   Option('--peakalignmethod', help = 'The type of a peak alignment method.', values=[\
                        Value('NN', help = 'Nearest Neightbour', parameters=[\
                            Option('--units', help = 'M/Z units', values=[\
                                Value('Da', help = 'Dalton units', parameters=[
                                    Option('--maxpeakshift', help = 'Maximum allowed peak shift.', values=[0.1, None], type=float, conditions=[('maxpeakshift>=0.0', lambda x: x.maxpeakshift>=0.0)], targets=['/palignparams'])\
                                    ]),
                                Value('ppm', help = 'ppm units', parameters=[
                                    Option('--maxpeakshift', help = 'Maximum allowed peak shift.', values=[100.0, None], type=float, conditions=[('maxpeakshift>=0.0', lambda x: x.maxpeakshift>=0.0)], targets=['/palignparams'])\
                                    ]),
                                Value('CCS', help = 'drift time', parameters=[])
                                ], type=str, targets=['/mzinfo']),\
                            ])
                        ], type=str),\

                   Option('--cmzmethod', help = 'The type of an estimation method for reference mz or drift feature vector.', values=[\
                        Value('density', help = 'Density method.', parameters=[\
                            Option('--units', help = 'M/Z units.', values=[\
                                Value('Da', help = 'Dalton units', parameters=[
                                    Option('--resolution', help = 'Resolution.', values=[0.01, None], type=float, conditions=[('resolution>=0.0', lambda x: x.resolution>=0.0)], targets=['/cmzparams'])\
                                    ]),
                                Value('ppm', help = 'ppm units', parameters=[
                                    Option('--resolution', help = 'Resolution.', values=[10.0, None], type=float, conditions=[('resolution>=0.0', lambda x: x.resolution>=0.0)], targets=['/cmzparams'])\
                                    ]),
                                Value('CCS', help = 'drift time', parameters=[])                                    
                                ], type=str, targets=['/mzinfo']),\
                            ])
                        ], type=str),\
                   #Option('--min_mz', help = 'Lower limit for M/Z values', values=[0.0, None], type=float, conditions=[('min_mz>=0.0', lambda x: x.min_mz>=0.0)]),\
                   #Option('--max_mz', help = 'Upper limit for M/Z values', values=[50000.0, None],  type=float, conditions=[('max_mz>=0.0', lambda x: x.max_mz>=0.0),
                   #                                                                                                         ('max_mz>=min_mz', lambda x: x.max_mz>=x.min_mz)]),\
        
                   Option('dbfile', help =  'HDF5-based database file with multiple deposited msi datasets '+\
                                            'for peak alignment. Each MSI dataset assumes to contain peak picked '+\
                                            'MSI data with their spatial coordinates and mz (and optionally drift time) '+\
                                            'feature vectors.', \
                                            values=['', None], type=str, optional=False),
                   Option('dbprocfile', help = 'Path to a HDF5-based msi database for storage and organization of pre-processed data. '+\
                                               'if this db file exists, all pre-processing parameters will be extracted from it, '+\
                                               'to make sure that the newly imported data are compatible with the ones stored in '+\
                                               'the processed db. The pre-processing workflow can be customized for newly '+\
                                               'created db instance.',\
                                               values=['', None], type=str)
                   #Option('--datasetids', help = 'The paths to datasets in a database file', values=['', None], type=str),\
                   ]; #end of value PeakAlign



                        
IntraNorm_options=[\
                   Option('--method', help = 'Intranormalization method.', values=[\
                        Value('mfc', help = 'Median fold change', parameters=[\
                            Option('--reference', help = 'Refence dataset with respect to which the fold intensity changes of other datasets are calculated.', \
                                values=['mean'], type=str, targets=['/params']),\
                            Option('--offset', help = 'Disregard peak intensity smaller that this value.', \
                                values=[0, None], type=float, conditions=[('offset>=0.0', lambda x: x.offset>=0.0)], targets=['/params']),\
                            Option('--outliers', help = 'Outliers.', \
                                values=['yes', 'no'], type=str, targets=['/params']),\
                            ]),\

                        Value('mean', help = 'Mean', parameters=[\
                            Option('--offset', help = 'Disregard peak intensity smaller that this value.', \
                                values=[0, None], type=float, conditions=[('offset>=0.0', lambda x: x.offset>=0.0)], targets=['/params']),\
                            Option('--outliers', help = 'Outliers.', \
                                values=['yes', 'no'], type=str, targets=['/params']),\
                            ]),\
                        
                        Value('median', help = 'Median', parameters=[\
                            Option('--offset', help = 'Disregard peak intensity smaller that this value.', \
                                values=[0, None], type=float, conditions=[('offset>=0.0', lambda x: x.offset>=0.0)], targets=['/params']),\
                            Option('--outliers', help = 'Outliers.', \
                                values=['yes', 'no'], type=str, targets=['/params']),\
                            ])\
                        ], type=str),\
                   Option('--min_mz', help = 'Lower limit for M/Z values.', values=[0.0, None], type=float, conditions=[('min_mz>=0.0', lambda x: x.min_mz>=0.0)]),\
                   Option('--max_mz', help = 'Upper limit for M/Z values.', values=[50000.0, None],  type=float, conditions=[('max_mz>=0.0', lambda x: x.max_mz>=0.0)]),\
                   
                   Option('h5dbname', help = 'Input HDF5 file name(s).', is_list=False, values=['', None], type=str, optional=False)                   
                   ]; #end of value IntraNorm

InterNorm_options=[\
                   Option('--method', help = 'Inter-normalization method.', values=[\
                        Value('mfc', help = 'Median fold change', parameters=[\
                            Option('--reference', help = 'Refence dataset with respect to which the fold intensity changes of other datasets are calculated.', \
                                values=['mean'], type=str, targets=['/params']),\
                            Option('--offset', help = 'Disregard peak intensity smaller that this value.', \
                                values=[0.0, None], type=float, conditions=[('offset>=0.0', lambda x: x.offset>=0.0)], targets=['/params']),\
                            Option('--outliers', help = 'Outliers.', \
                                values=['yes', 'no'], type=str, targets=['/params']),\
                            ]),\

                        Value('mean', help = 'Mean', parameters=[\
                            Option('--offset', help = 'Disregard peak intensity smaller that this value.', \
                                values=[0.0, None], type=float, conditions=[('offset>=0.0', lambda x: x.offset>=0.0)], targets=['/params']),\
                            Option('--outliers', help = 'Outliers.', \
                                values=['yes', 'no'], type=str, targets=['/params']),\
                            ]),\
                        
                        Value('median', help = 'Median', parameters=[\
                            Option('--offset', help = 'Disregard peak intensity smaller that this value.', \
                                values=[0.0, None], type=float, conditions=[('offset>=0.0', lambda x: x.offset>=0.0)], targets=['/params']),\
                            Option('--outliers', help = 'Outliers.', \
                                values=['yes', 'no'], type=str, targets=['/params']),\
                            ])\
                        ], type=str),\
                   Option('--min_mz', help = 'Lower limit for M/Z values.', values=[0.0, None], type=float, conditions=[('min_mz>=0.0', lambda x: x.min_mz>=0.0)]),\
                   Option('--max_mz', help = 'Upper limit for M/Z values.', values=[50000.0, None],  type=float, conditions=[('max_mz>=0.0', lambda x: x.max_mz>=0.0)]),\
                   
                   Option('h5dbname', help = 'Input HDF5 file name(s).', is_list=False, values=['', None], type=str, optional=False)
                   ]; #end of value InterNorm

VST_options=[\
                   Option('-m,--method', help = 'Variance stabilizing transformation.', values=[\
                        Value('started-log', help = 'StartedLog', parameters=[\
                            Option('--offset', help = 'Offset.', values=[50.0, None], type=float, targets=['/params']),\
                            ])\
                        ], type=str),\
                   Option('h5dbname', help = 'Input HDF5 file name(s).', is_list=False, values=['', None], type=str, optional=False)
                   ]; #end of value VST
    

#    Option('--merge', help = 'Merge HDF5 files into the first HDF5 file', values=['yes', 'no'], type=str),\
#    Option('--rawdatagroup', help = 'Group to assign to raw HDF5 input data', values=['train', None], type=str),\

#    Option('h5dbname', help = 'Input HDF5 file name(s)', is_list=True, values=['', None], type=str, optional=False),\


'''

