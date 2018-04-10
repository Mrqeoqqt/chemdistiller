# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 16:21:30 2016

@author: ilaponog
"""
from os.path import isfile;
import numpy as np;


__base64encoderrecord='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=';

def __decodequad(dquad):
    #print(dquad[0],dquad[1],dquad[2],dquad[3]);
    a=dquad[0];
    b=dquad[1];
    c=dquad[2];
    d=dquad[3];
    #print(a,b,c,d);
    i=(a<<2)|((b & 0b00110000 ) >>4);
    j=((b & 0b00001111) <<4) | ((c & 0b00111100)>>2);
    k=((c & 0b00000011)<<6) | ((d & 0b00111111));
    return bytearray([i,j,k]);
        
    
    
    
    
    

def __convtripl(tripl):
    i=tripl[0];
    j=tripl[1];
    k=tripl[2];
    #print(i,j,k);
    a=(i & 0b11111100)>>2;
    b=((i & 0b00000011)<<4) | ((j & 0b11110000) >>4);
    c=((j & 0b00001111) <<2) | ((k & 0b11000000) >>6);
    d=((k & 0b00111111));
    #print(a,b,c,d);
    return __base64encoderrecord[a]+__base64encoderrecord[b]+\
    __base64encoderrecord[c]+__base64encoderrecord[d];
    

    
    
    

def EncodeToBase64(btarray):
    '''
    EncodeToBase64(btarray: bytearray):string;

    - Converts bytearray input into Base64 encoded string, padding as necessary;
    '''
    bb=bytearray(btarray);
    j=len(bb)//3;
    k=len(bb)%3;
    #print(bb,j,k);
    encstring='';
    for i in range(0,j):
        tripl=bb[i*3:i*3+3];
        ss=__convtripl(tripl);
      #  print(tripl,ss);
        encstring=encstring+ss;
      #  print(ss);
    
    if k==1:
        tripl=bytearray([bb[j*3],0,0]);
      #  print(tripl);
        ss=__convtripl(tripl);
      #  print(ss);
        encstring=encstring+ss[0:2];
        encstring=encstring+'==';
      #  print(encstring);
    elif k==2:
        tripl=bytearray([bb[j*3],bb[j*3+1],0]);
        ss=__convtripl(tripl);
      #  print(ss);
        encstring=encstring+ss[0:3];
        encstring=encstring+'=';        
    
    
    
    return encstring;
    

def DecodeFromBase64(encstring):
    '''
    DecodeFromBase64(encstring:string):bytearray

    - Decodes Base64 encoded string into bytearray. String length must be of 
    multiple of 4 and it should only contain the following symbols:
    
    ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=

    with '=' indicating the padding at the end of the string only.
    '''
    encstring=encstring.encode('ascii');
    if len(encstring)%4!=0:
        raise NameError('String length not divisible by 4!');
    if len(encstring)==0:
        return bytearray();
    bb=bytearray();
    q=len(encstring)//4;
    try:
        for i in range(0,len(encstring)):
            bb.append(__base64encoderrecord.index(encstring[i]));
    except: raise NameError('Non-standard characters in Base64 string!');
    
    #for i in bb:
        #print(i);
        
    decbb=bytearray();
    
    for i in range(0, q):
        dquad=bb[i*4:i*4+4];
        qq=__decodequad(dquad);
        decbb=decbb+qq;
    if encstring[len(encstring)-1]=='=':
        decbb.pop();
    if encstring[len(encstring)-2]=='=':
        decbb.pop();
    
    return decbb;
    




#%%

def importfpt(fname):
    bits=[];
    try:
        if isfile(fname):
            f=open(fname,'r');
            for s in f:
                s=s.replace('\n','');
                s=s.split(',');
                for i in s:
                    bits.append(int(i));
            f.close();
        #print(len(bits));
        bb=np.packbits(bits);
        #print(bb);
        #print(len(bb));
        if len(bb)!=1427:
            raise NameError('Error reading fingerprint');
        else:
            return bb;
        
    except:
        #print('Error reading fingerprint');
        #raise NameError('Error reading fingerprint');
        return [];


import os;
from glob import glob;

inpath= os.getcwd();
print(inpath);



fileslist = [y for x in os.walk(inpath) for y in glob(os.path.join(x[0], '*.fp'))];
if os.path.isfile(inpath+'/TotalFingerprints.txt'):
    fout=open(inpath+'/TotalFingerprints.txt','a');
else: 
    fout=open(inpath+'/TotalFingerprints.txt','w');
for fname in fileslist:
    bb=importfpt(fname);
    if len(bb)>0:
        ss=EncodeToBase64(bb);
    ff=os.path.split(fname)[1];
    ff,ext=os.path.splitext(ff);
    print(ff);
    fout.write('%s\t%s\n'%(ss,ff));
fout.close();


