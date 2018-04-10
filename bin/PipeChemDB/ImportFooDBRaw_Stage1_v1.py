# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 14:08:15 2016

@author: ilaponog
"""

infile='e:/RawData/FooDB/compounds.txt';
outpath='e:/Imperial/tempDB/FooDB';
outfile=outpath+'/totalentries.set';
outerror=outpath+'/totalentries.error';
outnames=outpath+'/totalentries.names';
nproc=1;

import os;

import sys;

if sys.byteorder!='little':
    print('Only little endian machines supported! bye bye ....');
    quit();

if not os.path.exists(outpath):
    os.makedirs(outpath);

import pybel;

def initgloballock(l):
    global lock;
    lock = l;

#from multiprocessing import Pool,Lock;

#lock = Lock();
#pool = Pool(initializer=initgloballock, initargs=(lock,),processes=nproc);


def ProcessDBList(inputrecord):
    try:
        mol=pybel.readstring('smi',inputrecord[0]);
        mol.addh();
        mass=mol.exactmass;
        charge=mol.charge;
        if charge!=0:
            mz=abs(mass/charge);
            charged=True;
            if charge>0:
                chargestate='Positive';
            else:
                chargestate='Negative';
        else:
            mz=mass;
            charged=False;
            chargestate='Neutral';
            
        inchi=mol.write('inchi').replace('\n','').replace('InChI=1S/','');
        sinchi=inchi.split('/');
        shortinchi=sinchi[0];
        for j in range(1,len(sinchi)):
                    if sinchi[j][0]=='c' or sinchi[j][0]=='h':
                        shortinchi+='/'+sinchi[j];
        if charged:        
            outstring='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(chargestate,mz,mass,charge,shortinchi,inchi,inputrecord[0],inputrecord[1]);
        else:
            outstring='%s\t%s\t%s\t%s\t%s\t%s'%(chargestate,mz,shortinchi,inchi,inputrecord[0],inputrecord[1]);
        return [0,outstring];
    except:
        return [-1,inputrecord[0],inputrecord[1]];

cc=0;
finp=open(infile,'r');
fout=open(outfile,'w');
ferr=open(outerror,'w');
fnames=open(outnames,'w');

dblist=[];
totcount=0;
for s in finp:
    print(s);
    totcount+=1;
    if totcount>1:
        s=s.replace('\n','');
        s=s.split('\t');
        names=s[2];
        smi=s[1];
        if smi!='NULL':
            idx=int(s[0]);
            fnames.write('%s\t%s\n'%(idx,names));
            smi=smi.split('.');
            for ss in smi:
                dblist.append([ss,idx]);
            cc+=1;
            if cc>10000:
                print('Processing....');
                cc=0;
                for db in dblist:
                    i=ProcessDBList(db);
                    if i[0]==0:
                        fout.write('%s\n'%i[1]);
                        #print(i[1]);
                    else:
                        ferr.write('%s\t%s\n'%(i[1],i[2]));
                dblist=[];
if cc>0:
    print('Processing....');
    cc=0;
    for db in dblist:
       i=ProcessDBList(db);
       if i[0]==0:
                fout.write('%s\n'%i[1]);
                #print(i[1]);
       else:
                ferr.write('%s\t%s\n'%(i[1],i[2]));
    dblist=[];


finp.close();
fout.close();
ferr.close();
fnames.close();
#pool.close();
#pool.join();

print('Finished. Total records: %s'%totcount);


