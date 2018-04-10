# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:40:24 2017

@author: ilaponog
"""

#def print_in_the_same_line(s):
#    pass

import sys;
import os;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);




#print(sys.version_info.major)
if sys.version_info.major==3:
    from chemdistiller.utils.sysutils_p3 import print_in_the_same_line_p3;
    print_in_the_same_line=print_in_the_same_line_p3;
else:
    from chemdistiller.utils.sysutils_p2 import print_in_the_same_line_p2;
    print_in_the_same_line=print_in_the_same_line_p2;

python_version=sys.version_info.major;    
        
if __name__=='__main__':
    print_in_the_same_line(sys.version_info)
        
        