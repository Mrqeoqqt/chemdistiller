# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 20:43:45 2016

@author: ilaponog
"""

import sys;
if sys.byteorder!='little':
    print('Only little endian machines supported! bye bye ....');
    quit();

from glob import glob;
import os;
from operator import itemgetter;   
#subpaths=['ZincEMolecules',\
#'FooDB',\
#'PubChem',\
#'UNPD'];
   
subpaths=['PubChem'\
];



currentsubf='';

dbfile=[];

def closecurrentsubf():
    global currentsubf;
    global dbfile;
    global outpath;
    if currentsubf!='':
        fout=open(outpath+'/'+currentsubf+'.st2','w');
        if 'Neutral' in currentsubf:
            for s in dbfile:
                fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(s[0],s[1],s[2],s[3],s[4],s[5],s[6]));
        else:
            for s in dbfile:
                fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9]));
        fout.close();
        currentsubf='';
        dbfile=[];
        



def loadcurrentsubf():
    global currentsubf;
    global dbfile;
    global indbpath;
    global sourcepath;
    global missingcount;
    if currentsubf!='':
        dbfile=[];
        if os.path.isfile(indbpath+'/'+currentsubf+'.st2'):
            finp=open(indbpath+'/'+currentsubf+'.st2','r');
            for s in finp:
                ss=s.replace('\n','').split('\t');
                dbfile.append(ss);
            finp.close();
        if os.path.isfile(sourcepath+'/'+currentsubf+'.st2'):
            finp=open(sourcepath+'/'+currentsubf+'.st2','r');
            for s in finp:
                ss=s.replace('\n','').split('\t');
                dbfile.append(ss);
            finp.close();
        if 'Neutral' in currentsubf:
            dbfile.sort(key=itemgetter(0,1,2,4));
        else:
            dbfile.sort(key=itemgetter(0,1,2,3,4,6));


for subpath in subpaths:
    sourcepath='e:/Imperial/ready/'+subpath;
    indbpath='e:/DB/'+subpath;
    outpath='e:/Imperial/ready/patch/'+subpath;
    
    result = [y for x in os.walk(sourcepath) for y in glob(os.path.join(x[0], '*.st2'))];
    for fname in result:
        fname=fname.replace('\\','/');
        dirname,ff=os.path.split(fname);
        dirname=dirname.replace(sourcepath,outpath);
        if not os.path.exists(dirname):
            os.makedirs(dirname);

        currentsubf=os.path.splitext(fname.replace(sourcepath+'/',''))[0];
        print(currentsubf);
        loadcurrentsubf();
        #Write duplicate search here
        closecurrentsubf();
        
        
    
    
    
    
    

