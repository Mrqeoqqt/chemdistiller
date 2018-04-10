# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 14:08:15 2016

@author: ilaponog
"""
#%%
import os;

import sys;

if sys.byteorder!='little':
    print('Only little endian machines supported! bye bye ....');
    quit();
from glob import glob;


subpaths=['PubChem'];

#%%
for subpath in subpaths:
    inpath='e:/Imperial/tempDB/'+subpath;
    outpath=inpath;
    
    outerror=outpath+'/totalentries.error2';
    
    
    if not os.path.exists(outpath):
        os.makedirs(outpath);
    
    
    fileslist = [y for x in os.walk(inpath) for y in glob(os.path.join(x[0], '*.set'))];
    
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
    from operator import itemgetter;
    
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
    
    #dblist.sort(key=itemgetter(0,1,2,3));
    dbsubi=1;
    for i in range(0,len(dblist)-1):
        db=dblist[i];cc=0;
        for j in range(i+1,len(dblist)):
            if db[0]==dblist[j][0]:
                if db[0]=='Neutral':
                    offs=0;
                else:
                    offs=2;
                if db[5+offs]==dblist[j][5+offs]:
                    cc+=1;
                    dblist[j][5+offs]=dblist[j][5+offs].replace(',','#%s,'%(cc+1))+'#%s'%(cc+1);
            else:
                break;
        if cc>0:
            db[5+offs]=db[5+offs].replace(',','#%s,'%(1))+'#%s'%(1);
            
                   

    print('Renamed multiples');
    print('New Count: %s'%len(dblist));
    #dblist.sort(key=itemgetter(0,1,2,3));
    
    
    
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
    
    if not os.path.exists(outpath+'/FPTInput'):
        os.makedirs(outpath+'/FPTInput');
    
    fptout=open(outpath+'/FPTInput/%s.smi'%fptbatch,'w');
    
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
            fptout=open(outpath+'/FPTInput/%s.smi'%fptbatch,'w');
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
        fname=outpath+'/%s/%s/%s/%s/%s.st'%(db[0],i,j,k,l);
        fpath=outpath+'/%s/%s/%s/%s'%(db[0],i,j,k);
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
    
    
