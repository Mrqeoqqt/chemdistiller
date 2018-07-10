# -*- coding: utf-8 -*-
"""
==============================================================================
        ChemDistiller File Parser Script
==============================================================================
    This is an external script that generates the rank for each candidate
    that was processed.
    Dependency: chemdistiller.utils.inchi.inchikey_from_inchi

python RankTestJson.py <input_json_file> mode <outputfile>
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


def seekRank(df, mode):
    # inchi.inchikey_from_inchi
    # usage : the same way in html_report.generate_candidate_html
    assert mode in {0, 1}
    lst = []
    for index, row in df.iterrows():
        target_inchikey = ''
        if pd.isnull(row['parameters.inchi']) | \
                len(str(row['parameters.inchi']).strip()) == 0:
            if 'parameters.inchikey' not in row.index:
                lst.append((row['parameters.db_molecule_name'], 'No inchi/inchikey'))
                continue
            else:
                target_inchikey = row['parameters.inchikey']
        else:
            inchi = 'InChI=1S/%s' % row['parameters.inchi']
            target_inchikey = inchikey_from_inchi(inchi)
        content = row['peaks'][0]
        annotations = []
        if 'annotations' in content.keys():
            annotations = content['annotations'][0]
        else:
            if mode == 1:
                lst.append((row['parameters.db_molecule_name'], 'No annotations'))
            elif mode == 0:
                lst.append((target_inchikey, 'No annotations'))
            continue
        candidates = annotations['mol_candidates']
        indice = 0
        i = 0
        for c in candidates:
            inchi = 'InChI=1S/%s' % c['InChI']
            inchikey = inchikey_from_inchi(inchi)
            i += 1
            if target_inchikey == inchikey:
                indice = i
                break
        if mode == 1:
            lst.append((row['parameters.db_molecule_name'], indice))
        elif mode == 0:
            lst.append((target_inchikey, indice))
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


def run(jsonname, mode, output_file):
    """

    :param jsonname: input for rank
    :param output_file: result <No.> <InChIKey> <Rank>
    :return: None
    """
    if os.path.splitext(os.path.basename(jsonname))[1] != '.json':
        return
    df = getDF(jsonname)
    lst = seekRank(df, mode)
    makeFile(lst, output_file)
    return


if __name__ == "__main__":
    # copy from import.py
    try:
        start = datetime.datetime.now()
        if sys.byteorder != 'little':
            print('Only little endian machines currently supported! bye bye ....')
            quit()
        # deal with parameters
        if len(sys.argv) <= 3:
            print("Not enough parameters!")
            print("python MSP2ChemDistiller.py <input_json_file> mode <outputfile>")
            quit()
        elif len(sys.argv) == 4:
            run(sys.argv[1], int(sys.argv[2]), sys.argv[3])
        else:
            print("parameter error...\n")
            print("python MSP2ChemDistiller.py <input_json_file> mode <outputfile>")
            quit()
        finish = datetime.datetime.now()
        print("Time cost:%f s" % (finish - start).seconds)
    except AssertionError as a:
        print('''parameter [mode] accept values 0, 1.
                0 indicates using INCHIKEY in output,
                1 indicates using db_molecule_name in output.''')
    except ValueError as v:
        print("Value Error!")
