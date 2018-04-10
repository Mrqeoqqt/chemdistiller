# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 01:14:07 2016

@author: ilaponog
"""

subpaths=['PubChem'];

import os;

import sys;

if sys.byteorder!='little':
    print('Only little endian machines supported! bye bye ....');
    quit();
from glob import glob;
from operator import itemgetter;   
currentsubf='';

dbfile=[];

changed=False;

def closecurrentsubf(dirname):
    global currentsubf;
    global dbfile;
    global outpath;
    global changed;
    if currentsubf!='':
        if changed:
            if not os.path.exists(dirname):
                os.makedirs(dirname);

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
    global changed;
    if currentsubf!='':
        dbfile=[];
        if os.path.isfile(indbpath+'/'+currentsubf+'.st2'):
            finp=open(indbpath+'/'+currentsubf+'.st2','r');
            for s in finp:
                ss=s.replace('\n','').split('\t');
                dbfile.append(ss);
            finp.close();
        changed=False;
        if 'Neutral' in currentsubf:
            dbfile.sort(key=itemgetter(0,1,2,4));
            offs=0;
        else:
            dbfile.sort(key=itemgetter(0,1,2,3,4,6));
            offs=2;
        for i in reversed(range(1,len(dbfile))):
            if dbfile[i][2+offs]==dbfile[i-1][2+offs]:
                dbfile[i-1][4+offs]+=','+dbfile[i][4+offs];
                dbfile.pop(i);
                changed=True;
            


for subpath in subpaths:
    indbpath='e:/DB/'+subpath;
    outpath='e:/Imperial/ready2/patch/'+subpath;
    
    result = [y for x in os.walk(indbpath) for y in glob(os.path.join(x[0], '*.st2'))];
    for fname in result:
        fname=fname.replace('\\','/');
        dirname,ff=os.path.split(fname);
        dirname=dirname.replace(indbpath,outpath);
        
        currentsubf=os.path.splitext(fname.replace(indbpath+'/',''))[0];
        print(currentsubf);
        loadcurrentsubf();
        #Write duplicate search here
        closecurrentsubf(dirname);




