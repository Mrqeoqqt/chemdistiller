# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 14:08:15 2016

@author: ilaponog
"""
#%%
inpath='e:/Imperial/tempDB/TestDB';
outpath='e:/Imperial/tempDB/TestDB';

from glob import glob;
import os;

import sys;

if sys.byteorder!='little':
    print('Only little endian machines supported! bye bye ....');
    quit();

fileslist = [y for x in os.walk(inpath) for y in glob(os.path.join(x[0], '*.set'))];

for i in range(len(fileslist)):
    fileslist[i]=fileslist[i].replace('\\','/');
print(fileslist);

subpaths=[];
for infile in fileslist:
    dirname,fname=os.path.split(infile);
    dirname=dirname.replace(inpath,'');
    print(dirname);
    if dirname not in subpaths:
        subpaths.append(dirname);

from operator import itemgetter;    

for dirname in subpaths:
    suboutpath=outpath+dirname;

    outerror=suboutpath+'/totalentries.error2';




    if not os.path.exists(suboutpath):
        os.makedirs(suboutpath);

#from glob import glob;

    fileslist = [y for x in os.walk(inpath+dirname) for y in glob(os.path.join(x[0], '*.set'))];

    dblist=[];
    idx=[];
    #totcount=-1;
    for fname in fileslist:
        print(fname);
        print('Reading...');
        finp=open(fname,'r');
        for s in finp:
            s=s.replace('\n','').split('\t');
            s[1]=float(s[1]);
            dblist.append(s);
        finp.close();
        print('Completed');
    
    print('Finished reading. Total records: %s'%len(dblist));
    print('Sorting now...');

    
    #%%
    dblist.sort(key=itemgetter(0,1,2,3));
    print('Sorting Finished.');
    #%%
    ferr=open(outerror,'w');
    for i in reversed(range(0,len(dblist))):
        db=dblist[i];
        if db[0]=='Neutral':
            if db[1]=='' or db[2]=='' or db[3]=='' or db[4]=='':
                ferr.write('%s\t%s\n'%(db[4],db[5]));
                dblist.pop(i);
        else:        
            if db[3]=='' or db[4]=='' or db[5]=='' or db[6]=='':
                ferr.write('%s\t%s\n'%(db[6],db[7]));
                dblist.pop(i);
            
    
    ferr.close();
    print('Removed bad records.');
    print('New Count: %s'%len(dblist));
    
    
    for i in reversed(range(1,len(dblist))):
        db=dblist[i];
        dbprev=dblist[i-1];
        if db[0]==dbprev[0]:
            if db[0]=='Neutral':
                if db[3]==dbprev[3]:
                    if dbprev[5]==db[5]:
                        dblist.pop(i);
                        #print(dbprev);
            else:        
                if db[5]==dbprev[5]:
                    if dbprev[7]==db[7]:
                        dblist.pop(i);
                        #print(dbprev);
    print('Removed duplicates');
    print('New Count: %s'%len(dblist));
    
    dbsubi=1;
    for i in range(1,len(dblist)):
        db=dblist[i];
        dbprev=dblist[i-1];
        #print(db);
        #print(dbprev);
        if db[0]==dbprev[0]:
            if db[0]=='Neutral':
                if db[5]==dbprev[5]:
                        dbsubi+=1;
                        db[5]=db[5]+'_%s'%dbsubi;
                else:
                    dbsubi=1;
            else:        
                if db[7]==dbprev[7]:
                        dbsubi+=1;
                        db[7]=db[7]+'_%s'%dbsubi;
                else:
                    dbsubi=1;
        else:
            dbsubi=1;
    print('Renamed multiples');
    print('New Count: %s'%len(dblist));
    
    
    
    for i in reversed(range(1,len(dblist))):
        db=dblist[i];
        dbprev=dblist[i-1];
        if db[0]==dbprev[0]:
            if db[0]=='Neutral':
                if db[3]==dbprev[3]:
                    dbprev[5]=dbprev[5]+','+db[5];
                    dblist.pop(i);
            else:        
                if db[5]==dbprev[5]:
                    dbprev[7]=dbprev[7]+','+db[7];
                    dblist.pop(i);
    print('Condensed duplicates.');
    print('New Count: %s'%len(dblist));
    
    print('Exporting for Fingerprint calculation...');
    
    fptbatch=0;
    fptcount=0;
    
    if not os.path.exists(suboutpath+'/FPTInput'):
        os.makedirs(suboutpath+'/FPTInput');
    
    fptout=open(suboutpath+'/FPTInput/%s.smi'%fptbatch,'w');
    
    for db in dblist:
        mz=db[1];
        #print(mz,db[0]);
        intmass=int(mz*100);
        i=intmass//1000;
        j=intmass%1000//100;
        k=intmass%100//10;
        l=intmass%10;
        fptcount+=1;
        if fptcount>10000:
            fptcount=0;
            fptbatch+=1;
            fptout.close();
            fptout=open(suboutpath+'/FPTInput/%s.smi'%fptbatch,'w');
        if db[0]=='Neutral':
            fptout.write('%s %s_%s_%s_%s_%s_%s\n'%(db[4],db[0],i,j,k,l,db[5][:50]));
        else:
            fptout.write('%s %s_%s_%s_%s_%s_%s\n'%(db[6],db[0],i,j,k,l,db[7][:50]));
    fptout.close();
    
    print('Exporting Preprocessed DB...');
    
    currentfile='';
    
    for db in dblist:
        mz=db[1];
        #print(mz,db[0]);
        intmass=int(mz*100);
        i=intmass//1000;
        j=intmass%1000//100;
        k=intmass%100//10;
        l=intmass%10;
        fname=suboutpath+'/%s/%s/%s/%s/%s.st'%(db[0],i,j,k,l);
        fpath=suboutpath+'/%s/%s/%s/%s'%(db[0],i,j,k);
        if fname!=currentfile:
            if currentfile!='':
                fout.close();
            if not os.path.exists(fpath):
                os.makedirs(fpath);
            fout=open(fname,'w');
            currentfile=fname;
        if db[0]=='Neutral':
            fout.write('%s\t%s\t%s\t%s\t%s\n'%(db[1],db[2],db[3],db[4],db[5]));
        else:
            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(db[1],db[2],db[3],db[4],db[5],db[6],db[7]));
    
    if currentfile!='':
        fout.close();
    print('Finished All');
    
    
