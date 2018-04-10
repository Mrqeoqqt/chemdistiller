# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 18:46:36 2016

@author: ilaponog
"""

import timeit;


print(timeit.timeit('pff.DoDefaultProcessing();','import ProcessFragmentFile_V4 as pff',number=1));
