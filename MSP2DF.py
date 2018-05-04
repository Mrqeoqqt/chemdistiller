# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 06:09:16 2018
This module aims at read MSP files to a pandas DataFrame,
and write the DataFrame back to MSP files.
@author:  Mingyi Xue
"""

import pandas as pd
import numpy as np
import os
import sys
import datetime
import re
from chemdistiller.utils.sysutils import print_in_the_same_line

def read_from_file(filename):
    """
    read all lines from an MSP file,
    store them in pandas DataFrame and return.
    :param filename: an input file
    :return: pandas.DataFrame
    """
    df = pd.DataFrame(data=[], index=[], columns=[])
    with open(filename, 'r') as f:
        i = 0
        print_in_the_same_line("processing file %s to dataframe" % filename)
        # initialize a dataframe and insert a column named 'peaks'
        for line in f:
            line = str(line.lstrip().rstrip())
            if len(line) == 0:
                continue
            lower = line.lower()
            if lower.startswith('name'):
                if i < 0:
                    # a new row/molecule begins
                    i = -i
                    i += 1
                else:
                    # in case mass spectrum has no peaks
                    i += 1
            if re.match(r"^[A-Za-z]", line):
                # if start with characters, not peaks
                if ':' in line:
                    # other info
                    s = line.split(':', 1)
                    df.loc[i, s[0]] = s[1]
            else:
                if i > 0:
                    df.loc[i, 'peaks'] = ''
                    i = -i
                    # a notation to recognize if df.loc[i, 'peaks'] is empty
                df.loc[-i, 'peaks'] = str(df.loc[-i, 'peaks']) + str(line) + '*'

                # peaks info, seperated by '*'
        print_in_the_same_line("loaded data from %s." % filename)
    return df


def write_to_file(df, filename):
    """
    write a DataFrame to textfile
    :param df: pandas.DataFrame
    :param filename: an output file
    :return: None
    """
    with open(filename, 'w') as f:
        print_in_the_same_line("writing dataframe to %s" % filename)
        # for (k, v) in df.iteritems():
        #     if k == 'peaks':
        #         continue
        #     else:
        #         f.write("%s: %s\r\n" % (k, v))
        columns = df.columns.tolist()
        # get column name of a dateframe
        for index in df.index:
            for column in columns:
                if column != 'peaks':
                    s = df.loc[index, column]
                    if (s != '') & (pd.isnull(s) == False):
                        # if character is null, do not write it to files
                        f.write(column+':' + str(s) + '\n')
                    else:
                        continue
                else:
                    continue
            peaks = str(df.loc[index, 'peaks'])
            peaks = peaks.split('*')
            for peak in peaks:
                if len(peak) != 0:
                    f.write(str(peak) + '\n')
            f.write('\n')

        print_in_the_same_line("loaded data to %s" % filename)
    return

def MSP_run(input,output_folder=''):
    """
    test read_from_file and write_to_file in this module
    :param input: an input MSP file or folder
    :param output_folder: an output folder,
    default in the same directory as the input file,
    or under the input directory.
    :return:
    """
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

    for file in files:
        # process each MSP file
        df = read_from_file(file)
        filename = os.path.split(file)[1]
        if len(output_folder) > 0:
            s = os.path.join(output_folder, filename)
            if os.path.exists(output_folder) == False:
                print("output folder does not exist, creating...")
                os.mkdir(output_folder)
            print("Output folder:%s" % output_folder)
            write_to_file(df, s)
        else:
            output_folder = os.path.join(input_folder, 'output')
            s = os.path.join(output_folder, filename)
            if os.path.exists(output_folder) == False:
                print("output folder does not exist, creating...")
                os.mkdir(output_folder)
            print("Output folder:%s" % output_folder)
            write_to_file(df, s)


if __name__ == '__main__':
    start = datetime.datetime.now()
    if len(sys.argv) == 1:
        print("Missing input file or folder.")
        print("Usage: python MSP2DF.py <input> <output_folder>")
        exit(1)
    elif len(sys.argv) == 2:
        MSP_run(sys.argv[1])
        finish = datetime.datetime.now()
        print("Finished.")
        print("time cost:%d s" % (finish-start).seconds)
    elif len(sys.argv) == 3:
        MSP_run(sys.argv[1], sys.argv[2])
        finish = datetime.datetime.now()
        print("time cost:%d s" % (finish - start).seconds)
        print("Finished.")
    else:
        print("Too many parameters.")
        print("Usage: python MSP2DF.py <input> <output_folder>")
        exit(1)
