# -*- coding: utf-8 -*-
"""
Created on Mon Aug 08 12:24:51 2016

@author: ilaponog
"""
#%%
inpath='e:/Imperial/tempDB/STITCH/2';

outpath=inpath+'/FragsnoH';


#
#

missingcount=0;
totalcount=0;
#batch=0;
subcount=0;
subindex=0;
from os.path import isfile;
import os;
#import FingerPrintUtils_v7 as fp;
#import numpy as np;
import pybel;
import PybelTest_v3 as pt;

#%%
    
#if not os.path.exists(outpath+'/%s'%batch):
#      os.makedirs(outpath+'/%s'%batch);

if not os.path.exists(outpath):
      os.makedirs(outpath);

#fout=open(outpath+'/%s/fraginput_%s.txt'%(batch,subindex),'w');
fout=open(outpath+'/fraginput_%s.txt'%(subindex),'w');

for i in range(0,3000):
      #if i%8==batch:
      for j in range(0,10):
        for k in range(0,10):
            for l in range(0,10):
                if isfile(inpath+'/Neutral'+'/%s/%s/%s/%s.st'%(i,j,k,l)):
                    print(inpath+'/Neutral'+'/%s/%s/%s/%s.st'%(i,j,k,l));
                    dbfile=[];
                    fin=open(inpath+'/Neutral'+'/%s/%s/%s/%s.st'%(i,j,k,l),'r');
                    for s in fin:
                        ss=s.replace('\n','').split('\t');
                        dbfile.append(ss);
                    fin.close();
                    totalcount+=len(dbfile);
                    print('Reading complete');
                    for s in dbfile:
                        smiles=s[3];
                        idx=s[4];
                        fullidx='Neutral_%s_%s_%s_%s_%s'%(i,j,k,l,idx);
                        mol=pybel.readstring('smi',smiles);
                        mol.addh();
                        model=pt.PybelModel_To_Fragmenter(mol);
                        subcount+=1;
                        if subcount>10000:
                            subcount=0;
                            fout.close();
                            subindex+=1;
                            fout=open(outpath+'/fraginput_%s.txt'%(subindex),'w');
                            #fout=open(outpath+'/%s/fraginput_%s.txt'%(batch,subindex),'w');
                        fout.write('%s\n'%fullidx);
                        print(fullidx);
                        fout.write('%s,%s,%s,%s,%s'%(model[0][0][0],model[0][0][1],model[0][0][2],model[0][0][3],model[0][0][4]));
                        for ii in range(1,len(model[0])):
                            fout.write('\t%s,%s,%s,%s,%s'%(model[0][ii][0],model[0][ii][1],model[0][ii][2],model[0][ii][3],model[0][ii][4]));
                        fout.write('\n');
                        if len(model[1])>0:
                            fout.write('%s,%s,%s,%s'%(model[1][0][0],model[1][0][1],model[1][0][2],model[1][0][3]));
                            for ii in range(1,len(model[1])):
                                fout.write('\t%s,%s,%s,%s'%(model[1][ii][0],model[1][ii][1],model[1][ii][2],model[1][ii][3]));
                            fout.write('\n');
                        
                        
                        
                        
                        
                        
                    
for i in range(0,3000):
      #if i%8==batch:
      for j in range(0,10):
        for k in range(0,10):
            for l in range(0,10):
                if isfile(inpath+'/Negative'+'/%s/%s/%s/%s.st'%(i,j,k,l)):
                    print(inpath+'/Negative'+'/%s/%s/%s/%s.st'%(i,j,k,l));
                    dbfile=[];
                    fin=open(inpath+'/Negative'+'/%s/%s/%s/%s.st'%(i,j,k,l),'r');
                    for s in fin:
                        ss=s.replace('\n','').split('\t');
                        dbfile.append(ss);
                    fin.close();
                    totalcount+=len(dbfile);
                    print('Reading complete');
                    for s in dbfile:
                        smiles=s[5];
                        idx=s[6];
                        fullidx='Negative_%s_%s_%s_%s_%s'%(i,j,k,l,idx);
                        mol=pybel.readstring('smi',smiles);
                        mol.addh();
                        model=pt.PybelModel_To_Fragmenter(mol);
                        subcount+=1;
                        if subcount>10000:
                            subcount=0;
                            fout.close();
                            subindex+=1;
                            fout=open(outpath+'/fraginput_%s.txt'%(subindex),'w');
                            #fout=open(outpath+'/%s/fraginput_%s.txt'%(batch,subindex),'w');
                        fout.write('%s\n'%fullidx);
                        print(fullidx);
                        fout.write('%s,%s,%s,%s,%s'%(model[0][0][0],model[0][0][1],model[0][0][2],model[0][0][3],model[0][0][4]));
                        for ii in range(1,len(model[0])):
                            fout.write('\t%s,%s,%s,%s,%s'%(model[0][ii][0],model[0][ii][1],model[0][ii][2],model[0][ii][3],model[0][ii][4]));
                        fout.write('\n');
                        if len(model[1])>0:
                            fout.write('%s,%s,%s,%s'%(model[1][0][0],model[1][0][1],model[1][0][2],model[1][0][3]));
                            for ii in range(1,len(model[1])):
                                fout.write('\t%s,%s,%s,%s'%(model[1][ii][0],model[1][ii][1],model[1][ii][2],model[1][ii][3]));
                            fout.write('\n');
                            

for i in range(0,3000):
      #if i%8==batch:
      for j in range(0,10):
        for k in range(0,10):
            for l in range(0,10):
                if isfile(inpath+'/Positive'+'/%s/%s/%s/%s.st'%(i,j,k,l)):
                    print(inpath+'/Positive'+'/%s/%s/%s/%s.st'%(i,j,k,l));
                    dbfile=[];
                    fin=open(inpath+'/Positive'+'/%s/%s/%s/%s.st'%(i,j,k,l),'r');
                    for s in fin:
                        ss=s.replace('\n','').split('\t');
                        dbfile.append(ss);
                    fin.close();
                    totalcount+=len(dbfile);
                    print('Reading complete');
                    for s in dbfile:
                        smiles=s[5];
                        idx=s[6];
                        fullidx='Positive_%s_%s_%s_%s_%s'%(i,j,k,l,idx);
                        mol=pybel.readstring('smi',smiles);
                        mol.addh();
                        model=pt.PybelModel_To_Fragmenter(mol);
                        subcount+=1;
                        if subcount>10000:
                            subcount=0;
                            fout.close();
                            subindex+=1;
                            fout=open(outpath+'/fraginput_%s.txt'%(subindex),'w');
                            #fout=open(outpath+'/%s/fraginput_%s.txt'%(batch,subindex),'w');
                        fout.write('%s\n'%fullidx);
                        print(fullidx);
                        fout.write('%s,%s,%s,%s,%s'%(model[0][0][0],model[0][0][1],model[0][0][2],model[0][0][3],model[0][0][4]));
                        for ii in range(1,len(model[0])):
                            fout.write('\t%s,%s,%s,%s,%s'%(model[0][ii][0],model[0][ii][1],model[0][ii][2],model[0][ii][3],model[0][ii][4]));
                        fout.write('\n');
                        if len(model[1])>0:
                            fout.write('%s,%s,%s,%s'%(model[1][0][0],model[1][0][1],model[1][0][2],model[1][0][3]));
                            for ii in range(1,len(model[1])):
                                fout.write('\t%s,%s,%s,%s'%(model[1][ii][0],model[1][ii][1],model[1][ii][2],model[1][ii][3]));
                            fout.write('\n');
                        


print(totalcount);
fout.close();