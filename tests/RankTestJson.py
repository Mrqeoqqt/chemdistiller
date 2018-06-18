# -*- coding: utf-8 -*-
"""
==============================================================================
        ChemDistiller File Parser Script
==============================================================================
    This is an external script that generates the rank for each candidate
    that was processed.

python MSP2ChemDistiller.py <input_json_file> <outputfile>
author: Mingyi Xue
Date: June/16 2018
"""

from pandas.io.json import json_normalize
import pandas as pd
import json
import os
import sys
import re
import datetime
from chemdistiller.utils.inchi import inchikey_from_inchi;


def getDF(jsonname):
    """
    load info from json to pd.df
    :param jsonname: input should be a jason file
    :return: pandas.DataFrame
    """
    data_str = open(jsonname).read()
    data_list = json.loads(data_str)
    df = json_normalize(data_list)
    return df


def seekRank(df):
    # inchi.inchikey_from_inchi
    # usage : the same way in html_report.generate_candidate_html
    lst = []
    for index, row in df.iterrows():
        target = row['parameters.db_molecule_name']
        content = row['peaks'][0]
        annotations = content['annotations'][0]
        candidates = annotations['mol_candidates']
        indice = 0
        i = 0
        for c in candidates:
            inchi = 'InChI=1S/%s' % c['InChI']
            inchikey = inchikey_from_inchi(inchi)
            i +=1
            if target == inchikey:
                indice = i
                break
        lst.append((target, indice))
    return lst


def makeFile(lst, output_file):
    with open(output_file, 'w') as f:
        f.write('No.' + '\t' + 'Target InChIKey' + '\t' + 'Rank' + '\n')
        for i in range(len(lst)):
            target = lst[i][0]
            indice = lst[i][1]
            if indice == 0:
                indice = 'None'
            f.write(str(i+1) + '\t' + str(target) + '\t' + str(indice) + '\n')
    return


def run(jsonname, output_file):
    if os.path.splitext(os.path.basename(jsonname))[1] != '.json':
        return
    df = getDF(jsonname)
    lst = seekRank(df)
    makeFile(lst, output_file)
    return


if __name__ == "__main__":
    # copy from import.py
    start = datetime.datetime.now()
    if sys.byteorder != 'little':
        print('Only little endian machines currently supported! bye bye ....')
        quit()
    # deal with parameters
    if len(sys.argv) <= 2:
        print("Files needed!")
        print("python MSP2ChemDistiller.py <input_json_file> <outputfile>")
        quit()
    elif len(sys.argv) == 3:
        run(sys.argv[1], sys.argv[2])
    else:
        print("parameter error...\n")
        print("python MSP2ChemDistiller.py <input_json_file> <outputfile>")
        quit()
    finish = datetime.datetime.now()
    print("Time cost:%f s" % (finish - start).seconds)
