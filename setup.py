"""
Distiller setup file.

@author: Dr. Ivan Laponogov

"""

import sys;
if sys.byteorder!='little':
    print('Only little endian machines currently supported! bye bye ....');
    quit();
import os;

print('Source path: %s'%os.getcwd());
sys.path.append(os.getcwd());


