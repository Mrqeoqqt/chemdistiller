# -*- coding: utf-8 -*-
"""
==============================================================================
        ChemDistiller File Parser Script
==============================================================================
    This is an external script for the purpose of converting of different
    formats of MS2 spectra into internal data format for ChemDistiller.

run python.exe MSP2ChemDistiller.py -h for full help about supported options

author: Mingyi Xue
"""

import io
import sys
import pandas as pd
import numpy as np
import scipy as sp


def import_spectra_from_msp():
    pass


def process_msp_to_cd():
    pass


def generate_cd_file():
    pass


if __name__ == "__main__":
    # copy from import.py
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....')
        quit()

    # deal with parameters

    #