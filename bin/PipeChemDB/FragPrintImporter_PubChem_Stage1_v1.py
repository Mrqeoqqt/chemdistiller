# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 20:40:40 2016

@author: ilaponog
"""

import os;

inpath='e:/Imperial/DB/PubChem';
outpath='f:/DB/PubChem';
fragpath='e:/Imperial/DB/PubChem/FragsnoH';

import sys;
if sys.byteorder!='little':
    print('Only little endian machines supported! bye bye ....');
    quit();

print('Hello! System is little endian, so continuing ;-) ');

currentsubf='';

dbfile=[];
idxs=[];
missingcount=0;

def closecurrentsubf():
    global currentsubf;
    global dbfile;
    global idxs;
    global outpath;
    global missingcount;
    if currentsubf!='':
        fout=open(outpath+'/'+currentsubf,'w');
        if 'Neutral' in currentsubf:
            for s in dbfile:
                fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(s[0],s[1],s[2],s[3],s[4],s[5],s[6]));
                if s[6]=='':
                    missingcount+=1;
        else:
            for s in dbfile:
                fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9]));
                if s[8]=='':
                    missingcount+=1;
            
            
        
        fout.close();
        
        
        
        currentsubf='';
        dbfile=[];
        idxs=[];



def loadcurrentsubf():
    global currentsubf;
    global dbfile;
    global inpath;
    global idxs;
    global missingcount;
    if currentsubf!='':
        dbfile=[];
        if os.path.isfile(outpath+'/'+currentsubf):
            finp=open(outpath+'/'+currentsubf,'r');
            for s in finp:
                ss=s.replace('\n','').split('\t');
                if 'Neutral' in currentsubf:
                    idxs.append(ss[4]);
                    if ss[6]=='':
                        missingcount-=1;

                else:
                    idxs.append(ss[6]);
                    if ss[8]=='':
                        missingcount-=1;
                dbfile.append(ss);
            finp.close();
        else:
            finp=open(inpath+'/'+currentsubf,'r');
            for s in finp:
                ss=s.replace('\n','').split('\t');
                if 'Neutral' in currentsubf:
                    ss.append('');
                    idxs.append(ss[4]);
                else:
                    ss.append('');
                    ss.append('');
                    idxs.append(ss[6]);
                dbfile.append(ss);
            finp.close();
    
    

def addtosubf(cid,ss):
    global currentsubf;
    global dbfile;
    global idxs;
    ii=idxs.index(cid);
    if 'Neutral' in currentsubf:
        dbfile[ii][6]=ss[0];
    else:
        dbfile[ii][8]=ss[0];
        dbfile[ii][9]=ss[1];
        
    


from glob import glob
result = [y for x in os.walk(fragpath) for y in glob(os.path.join(x[0], '*.frg'))];

for fname in result:
    fname=fname.replace('\\','/');
    finp=open(fname,'r');
    for s in finp:
        ss=s.replace('\n','');
        if ('Neutral' in ss) or ('Negative' in ss) or ('Positive' in ss):
            ss=ss.split('_');
            Cfolder=ss[0];
            i=ss[1];j=ss[2];k=ss[3];l=ss[4];cid=ss[5];
        else:
            ss=ss.replace(' ','').split('\t');
            ss.pop(0);
            subf=Cfolder+'/'+i+'/'+j+'/'+k+'/'+l+'.st2';
            if subf<>currentsubf:
                closecurrentsubf();
                currentsubf=subf;
                if not os.path.exists(outpath+'/%s/%s/%s/%s'%(Cfolder,i,j,k)):
                        os.makedirs(outpath+'/%s/%s/%s/%s'%(Cfolder,i,j,k));
                        print('Folder created: %s'%(outpath+'/%s/%s/%s/%s'%(Cfolder,i,j,k)));
                loadcurrentsubf();
            addtosubf(cid,ss);
    finp.close();
    
print('Missing: %s'%missingcount);
