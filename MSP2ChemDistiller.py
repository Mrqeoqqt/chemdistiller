# -*- coding: utf-8 -*-
"""
==============================================================================
        ChemDistiller File Parser Script
==============================================================================
    This is an external script for the purpose of converting of different
    formats of MS2 spectra into internal data format for ChemDistiller.


python MSP2ChemDistiller.py <input_file> (<output_folder>) (<sample_file>)

author: Mingyi Xue
"""

import os
import sys
import pandas as pd
import re
import datetime
import MSP2DF
from chemdistiller.utils.sysutils import print_in_the_same_line



def process_msp_to_cd(input, output_folder='', sample_file=''):
    """

    :param input: input directory or file
    :param output_folder: output directory or by default
    :return:
    """
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

    for file in files:
        # process molecules in each file under input folder
        # to a single chemdistiller input
        df = MSP2DF.read_from_file(file)
        # read molecules to a whole dataframe
        for index in df.index:
            # each row in df represents a molecule
            # write every row to a single chemdistiller input file
            name = df.loc[index, 'NAME'].split(';')[0]
            # use 'NAME' of molecule as chemdistiller input filename
            output_file = os.path.join(output_folder, name+'.txt')
            print("processing molecule : "+name)
            with open(output_file, 'w') as f:
                notation = 0
                # notation for peak information
                for line in sample:
                    if '[' not in line:
                        # directly write a line that need not be filled
                        f.write(line)
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
                            if pd.isnull(df.loc[index, 'InChIKey']):
                                f.write(s[0] + '=\n')
                            else:
                                f.write(s[0] + '=' + str(df.loc[index, 'InChIKey']).strip() + '\n')
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
                            f.write('\t1,'+str(df.loc[index, 'PRECURSORMZ']).strip()+',100.0\n')
                        else:
                            # write peaks and mark notation = 1
                            notation = -1
                            indice = 0
                            peaks = str(df.loc[index, 'peaks'])
                            peaks = peaks.split('*')
                            for peak in peaks:
                                if len(peak) != 0:
                                    p = peak.split()
                                    if len(p) >= 2:
                                        p1 = p[0].split()
                                        p2 = p[1].split()
                                        indice += 1
                                        f.write('\t\t\t\t' + str(indice) + ',' + str(p[0].strip()) +
                                                ',' + str(p[1].strip()) + '\n')
    print("Completed.")


def generate_cd_file(content, output_file):
    pass


if __name__ == "__main__":
    # copy from import.py
    start = datetime.datetime.now()
    if sys.byteorder != 'little':
        print('Only little endian machines currently supported! bye bye ....')
        quit()
    # deal with parameters
    if len(sys.argv) == 1:
        print("python MSP2ChemDistiller.py <input_file> (<output_folder>) (<sample_file>)")
        quit()
    elif len(sys.argv) == 2:
        process_msp_to_cd(sys.argv[1])
    elif len(sys.argv) == 3:
        process_msp_to_cd(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        process_msp_to_cd(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("parameter error...\n")
        print("python MSP2ChemDistiller.py <input_file> (<output_folder>) (<sample_file>)")
        quit()
    finish = datetime.datetime.now()
    print("Time cost:%f s" % (finish-start).seconds)