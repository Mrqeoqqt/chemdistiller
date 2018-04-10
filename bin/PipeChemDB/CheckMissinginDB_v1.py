# -*- coding: utf-8 -*-
"""
Created on Mon Aug 08 12:24:51 2016

@author: ilaponog
"""
#%%
inpath='f:/DB/PubChem';
outpath='f:/DB/PubChem/Missing';

missingcount=0;
import os;

if not os.path.exists(outpath):
     os.makedirs(outpath);

fout=open(outpath+'/missing.txt','w');
foutsmi=open(outpath+'/missing.smi','w');

for i in range(0,1300):
      for j in range(0,10):
        for k in range(0,10):
            for l in range(0,10):
                if os.path.isfile(inpath+'/Neutral/%s/%s/%s/%s.st2'%(i,j,k,l)):
                    print(inpath+'/Neutral/%s/%s/%s/%s.st2'%(i,j,k,l));
                    
                    
                    fin=open(inpath+'/Neutral/%s/%s/%s/%s.st2'%(i,j,k,l),'r');
                    for s in fin:
                        ss=s.replace('\n','').split('\t');
                        if ss[1]=='' or ss[2]=='' or ss[3]=='' or ss[4]=='' or ss[5]=='' or ss[6]=='':
                            spath='Neutral_%s_%s_%s_%s_%s'%(i,j,k,l,ss[4]);
                            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(spath,ss[0],ss[1],ss[2],ss[3],ss[4],ss[5],ss[6]));
                            if ss[6]=='':
                                foutsmi.write('%s %s\n'%(spath,ss[3]));
                                #foutsmi.write(s);
                                #foutsmi.write('\n');
                                #foutsmi.write(str(ss));
                                #foutsmi.write('\n');
                                
                            missingcount+=1;
                            if ss[1]=='':
                                fout.write('Missing InChi\n');
                            if ss[2]=='':
                                fout.write('Missing InChiFull\n');
                            if ss[3]=='':
                                fout.write('Missing Smiles\n');                                
                            if ss[5]=='':
                                fout.write('Missing Fingerprint\n');                                
                            if ss[6]=='':
                                fout.write('Missing Fragprint\n');                                

                    fin.close();
                    print('Reading complete');




for i in range(0,1300):
      for j in range(0,10):
        for k in range(0,10):
            for l in range(0,10):
                if os.path.isfile(inpath+'/Negative/%s/%s/%s/%s.st2'%(i,j,k,l)):
                    print(inpath+'/Negative/%s/%s/%s/%s.st2'%(i,j,k,l));
                    fin=open(inpath+'/Negative/%s/%s/%s/%s.st2'%(i,j,k,l),'r');
                    for s in fin:
                        ss=s.replace('\n','').split('\t');
                        if ss[1]=='' or ss[2]=='' or ss[3]=='' or ss[4]=='' or ss[5]=='' or ss[6]=='' or ss[7]=='' or ss[8]=='' or ss[9]=='':
                            spath='Negative_%s_%s_%s_%s_%s'%(i,j,k,l,ss[6]);
                            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(spath,ss[0],ss[1],ss[2],ss[3],ss[4],ss[5],ss[6],ss[7],ss[8],ss[9]));
                                
                            if ss[8]=='':
                                foutsmi.write('%s %s\n'%(spath,ss[5]));
                            missingcount+=1;
                            if ss[3]=='':
                                fout.write('Missing InChi\n');
                            if ss[4]=='':
                                fout.write('Missing InChiFull\n');
                            if ss[5]=='':
                                fout.write('Missing Smiles\n');                                
                            if ss[7]=='':
                                fout.write('Missing Fingerprint\n');                                
                            if ss[8]=='':
                                fout.write('Missing Fragprint\n');                                

                                
                            
                            missingcount+=1;
                    fin.close();
                    print('Reading complete');


for i in range(0,1300):
      for j in range(0,10):
        for k in range(0,10):
            for l in range(0,10):
                if os.path.isfile(inpath+'/Positive/%s/%s/%s/%s.st2'%(i,j,k,l)):
                    print(inpath+'/Positive/%s/%s/%s/%s.st2'%(i,j,k,l));
                    fin=open(inpath+'/Positive/%s/%s/%s/%s.st2'%(i,j,k,l),'r');
                    for s in fin:
                        ss=s.replace('\n','').split('\t');
                        if ss[1]=='' or ss[2]=='' or ss[3]=='' or ss[4]=='' or ss[5]=='' or ss[6]=='' or ss[7]=='' or ss[8]=='' or ss[9]=='':
                            spath='Positive_%s_%s_%s_%s_%s'%(i,j,k,l,ss[6]);
                            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(spath,ss[0],ss[1],ss[2],ss[3],ss[4],ss[5],ss[6],ss[7],ss[8],ss[9]));
                            if ss[8]=='':
                                foutsmi.write('%s %s\n'%(spath,ss[5]));
                            missingcount+=1;
                            if ss[3]=='':
                                fout.write('Missing InChi\n');
                            if ss[4]=='':
                                fout.write('Missing InChiFull\n');
                            if ss[5]=='':
                                fout.write('Missing Smiles\n');                                
                            if ss[7]=='':
                                fout.write('Missing Fingerprint\n');                                
                            if ss[8]=='':
                                fout.write('Missing Fragprint\n');                                
                            missingcount+=1;
                    fin.close();
                    print('Reading complete');

fout.close();                    
foutsmi.close();
                    
                    
    
print(missingcount);

                    