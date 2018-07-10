# -*- coding: utf-8 -*-
"""
==============================================================================
        ChemDistiller File Parser Script
==============================================================================
    This is an external script for the purpose of converting of different
    formats of MS2 spectra into internal data format for ChemDistiller.


python MSP2ChemDistiller.py <input> mode (<output_folder>) (<sample_file>)

author: Mingyi Xue
"""

import os
import sys
import pandas as pd
import re
import datetime
import MSP2DF
import time


def getName(df, index):
    name = ''
    if ';' in str(df.loc[index, 'NAME']):
        name = df.loc[index, 'NAME'].split(';')[0].strip()
    else:
        name = df.loc[index, 'NAME']
    return name


def getFileName(filename_dict, name):
    ## usually use molecule name as filename
    ## add underscore to molecule name if some molecule
    ## appears more than once
    fname = ''
    if name in filename_dict:
        filename_dict[name] += 1
        fname = name + '_' + str(filename_dict[name])
    else:
        filename_dict[name] = 0
        fname = name
    return fname


def write2file(df, index, mode, outputfile, sample):
    with open(outputfile, 'w') as f:
        notation = 0
        # notation for peak information
        for line in sample:
            if '[' not in line:
                # directly write a line that need not be filled
                f.write(line)
                if line.startswith('inchi'):
                    if not pd.isnull(df.loc[index, 'INCHIKEY']):
                        f.write('inchikey=' + str(df.loc[index, 'INCHIKEY'].strip()) + '\n')
                if line.lstrip().startswith('peaks'):
                    notation += 1
                    # notation = 1, write precursor peak
                    # notation = 2, write all peaks
                    # notation = -1, all peaks have been written
            elif re.match("^[A-Za-z]", line.lstrip()):
                # '[' in line and start with characters
                # means this line is key=[value]
                s = line.rstrip().split('=')
                head = s[0].lstrip()
                if head == 'db_molecule_name':
                    if mode == 0:
                        if pd.isnull(df.loc[index, 'INCHIKEY']):
                            f.write(s[0] + '=\n')
                        else:
                            f.write(s[0] + '=' + str(df.loc[index, 'INCHIKEY']).strip() + '\n')
                    elif mode == 1:
                        if pd.isnull(df.loc[index, 'NAME']):
                            f.write(s[0] + '=\n')
                        else:
                            f.write(s[0] + '=' + str(df.loc[index, 'NAME']).strip() + '\n')
                elif head == 'precursor_ion':
                    f.write(s[0] + '=' + str(df.loc[index, 'PRECURSORTYPE']).strip() + '\n')
                elif head == 'precursor_mz':
                    f.write(s[0] + '=' + str(df.loc[index, 'PRECURSORMZ']).strip() + '\n')
                elif head == 'ion_type':
                    f.write(s[0] + '=' + str(df.loc[index, 'PRECURSORTYPE']).strip() + '\n')
            else:
                # '[' in line and start with numbers
                # means this line is a peak
                if notation < 0:
                    # only write peaks once
                    # if peaks have been written, skip
                    continue
                if notation == 1:
                    # write precursor peak
                    f.write('\t1,' + str(df.loc[index, 'PRECURSORMZ']).strip() + ',100.0\n')
                else:
                    # write peaks and mark notation = 1
                    notation = -1
                    indice = 0
                    peaks = str(df.loc[index, 'peaks'])
                    # try to store peaks in list
                    peaks = peaks.split('\n')
                    for peak in peaks:
                        p = re.split(r"\s+", peak.strip())
                        if len(peak) >= 2:
                            mass, intensity = p[0], p[1]
                            indice += 1
                            f.write('\t\t\t\t' + str(indice) + ',' + str(mass) +
                                    ',' + str(intensity) + '\n')


def process_msp_to_cd(input, mode, output_folder='', sample_file=''):
    """
    :param input: input directory or file
    :param mode: 0/1
    :param output_folder: output directory or by default
    :return:
    """
    assert int(mode) in {0, 1}
    print("Input:"+input)
    files = []
    # store all MSP files' path
    input_folder = ''
    if os.path.isfile(input):
        # single file input
        files.append(input)
        input_folder = os.path.split(input)[0]
    elif os.path.isdir(input):
        # input is a folder
        files = [os.path.join(input, f) for f in os.listdir(input)
                 if os.path.isfile(os.path.join(input, f))]
        input_folder = input
    print("Input folder:"+input_folder)

    if len(output_folder) == 0:
        # using default output directory
        print('Using default output directory...')
        output_folder = os.path.join(input_folder, 'output')
    if os.path.exists(output_folder) == False:
        # create output directory if not exists
        print("output folder does not exist, creating...")
        os.mkdir(output_folder)
    print("output folder "+output_folder)

    if len(sample_file) == 0:
        chemdistiller_path = os.path.dirname(os.path.realpath(__file__))
        sample_file = os.path.join(chemdistiller_path, 'Input_Sample')
        sample_file = os.path.join(sample_file, '0_empty.txt')
        # find the default empty sample template under chemdistiller/Input_Sample/ folder
    sample = []
    # save contents in 0_empty.txt
    with open(sample_file, 'r') as f:
        sample = [line for line in f]
        print("Template loaded from " + sample_file)
    # copy sample file to memory, save io time

    print("loading files...")
    start = time.time()
    filename_dict = {}
    # store names and counts to avoid molecules with the same name
    for file in files:
        # process molecules in each file under input folder
        # to a single chemdistiller input
        df = MSP2DF.read_from_file(file)
        # read molecules to a whole dataframe
        for index in df.index:
            # each row in df represents a molecule
            # write every row to a single chemdistiller input file
            print('\r %d/%d molecules have been processed' % (index, len(df.index)), end='')
            name = getName(df, index)
            fname = getFileName(filename_dict, name)
            output_file = os.path.join(output_folder, fname+'.txt')
            write2file(df, index, int(mode), output_file, sample)
        print('\n')
    end = time.time()
    print("Completed. Parse time cost in total: %lf secs" % (end - start))


def generate_cd_file(content, output_file):
    pass


if __name__ == "__main__":
    # copy from import.py
    if sys.byteorder != 'little':
        print('Only little endian machines currently supported! bye bye ....')
        quit()
    # deal with parameters
    try:
        if len(sys.argv) <= 2:
            print("python MSP2ChemDistiller.py <input> <mode> (<output_folder>) (<sample_file>)")
            quit()
        elif len(sys.argv) == 3:
            process_msp_to_cd(sys.argv[1], sys.argv[2])
        elif len(sys.argv) == 4:
            process_msp_to_cd(sys.argv[1], sys.argv[2], sys.argv[3])
        elif len(sys.argv) == 5:
            process_msp_to_cd(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
        else:
            print("parameter error...\n")
            print("python MSP2ChemDistiller.py <input> <mode> (<output_folder>) (<sample_file>)")
            quit()
    except AssertionError as a:
        print('''parameter [mode] accept values 0, 1.
                0 indicates using INCHIKEY as db_molecule_name,
                1 indicates using NAME as db_molecule_name.
                db_molecule_name will be used as keys in html 
                report and output json file in chemdistiller.''')
