# -*- coding: utf-8 -*-
"""
Created on Thu Sep 08 22:24:46 2016

@author: ilaponog
"""

import FilterOutMissinginDB_v1 as DBF;


inpath='f:/DB/HMDB';
outpath='e:/DB/HMDB';
deletedpath='e:/DB/HMDB/Deleted';
DBF.PurgeDBRaw(inpath,outpath,deletedpath,'HMDB','Clean, not 100% complete',4);
