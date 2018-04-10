# -*- coding: utf-8 -*-
"""
Created on Thu Sep 01 15:09:44 2016

@author: ilaponog
"""

inpath='e:/Imperial/DB/PubChem/Stats/batches';

fptlength=11416;

import numpy as np;

fcounts=np.zeros((fptlength,),dtype=np.float32);
counts=np.zeros((fptlength,),dtype=np.int32);

finp=open(inpath+'/counts.python');
for s in finp:
    ss=s.replace('\n','').split('\t');
    index=int(ss[0]);
    counts[index]=int(ss[1]);
    fcounts[index]=float(ss[2]);
finp.close();

#%%
countmask=np.zeros((fptlength,),dtype=np.int8);
for i in range(0,fptlength):
    if (fcounts[i]<0.05) or (fcounts[i]>0.95):
        countmask[i]=1;
        
#%%
finp=open(inpath+'/cross.total_xor','rb');
fcorr=np.fromfile(finp,dtype=np.float32);
finp.close();
fcorr=fcorr.reshape((fptlength,fptlength));

corrmask=np.zeros((fptlength,),dtype=np.int8);
for i in range(0,fptlength-1):
    print(i);
    if corrmask[i]==0:
        for j in range(i+1,fptlength):
            if fcorr[i,j]>0.95 or fcorr[i,j]<0.05 :
                corrmask[j]=1;
                
                

#%%
finp=open(inpath+'/cross.total_pearson','rb');
fcorr=np.fromfile(finp,dtype=np.float32);
finp.close();
fcorr=fcorr.reshape((fptlength,fptlength));

corrmaskpearson=np.zeros((fptlength,),dtype=np.int8);
for i in range(0,fptlength-1):
    print(i);
    if corrmaskpearson[i]==0:
        for j in range(i+1,fptlength):
            if abs(fcorr[i,j])>0.5:
                corrmaskpearson[j]=1;

#%%
print(fptlength-np.sum(corrmask));
print(fptlength-np.sum(corrmaskpearson));
print(fptlength-np.sum(countmask));

#%%
print(fptlength-np.sum(np.bitwise_or(corrmask,corrmaskpearson)));
print(fptlength-np.sum(np.bitwise_or(corrmask,countmask)));
print(fptlength-np.sum(np.bitwise_or(corrmaskpearson,countmask)));
print(fptlength-np.sum(np.bitwise_or(np.bitwise_or(corrmask,corrmaskpearson),countmask)));





