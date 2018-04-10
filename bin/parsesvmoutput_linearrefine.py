# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 15:30:02 2016

@author: Dr. Ivan Laponogov
"""
import os;
import numpy as np;
fout=open('do_models_linearrefine.bat','w');
fout.write('#!/bin/bash\n');
if os.path.isfile('linear_weighted.out'):
    fout.write('echo linear weighted \n');
    
    finp=open('linear_weighted.out','r');
    lines=finp.readlines();
    finp.close();

    bestresult=0;
    bestc=0;
    datafile='';
    w0=0;
    w1=0;
    currentc=0;
    for s in lines:
        s=s.replace('%','').rstrip('\n').replace('\r','');
        if s.startswith('./svm-train'):
            s=s.split(' ');
            for i in range(len(s)):
                if s[i]=='-w0':
                    w0=s[i+1];
                elif s[i]=='-w1':
                    w1=s[i+1];
                elif s[i]=='-c':
                    currentc=s[i+1];
            datafile=s[len(s)-1];
            
        elif s.startswith('Cross Validation Accuracy'):
            s=s.split(' = ')[1];
            
            s=s.replace('%','');
            
            s=float(s);
            if s>bestresult:
                bestresult=s;
                bestc=currentc;
    fout.write('echo Best Score=%s Best c=%s\n'%(bestresult,bestc));                
    fout.write('echo Linear Classifier Weighted >>linear_weighted.out\n');
    powc=int(np.log2(float(bestc)));
    for i in range(powc-2,powc+3):
        fout.write('echo ./svm-train -t 0 -c %s -v 5 -w0 %s -w1 %s %s >>linear_weighted.out\n'%(2**i,w0,w1,datafile));
        fout.write('./svm-train -t 0 -c %s -v 5 -w0 %s -w1 %s %s >>linear_weighted.out\n'%(2**i,w0,w1,datafile));
            

    #fout.write('./svm-train -t 0 -c %s -w0 %s -w1 %s %s %s.linear_weight\n'%(bestc,w0,w1,datafile,datafile));


if os.path.isfile('linear_notweighted.out'):
    fout.write('echo linear notweighted \n');
    
    finp=open('linear_notweighted.out','r');
    lines=finp.readlines();
    finp.close();

    bestresult=0;
    bestc=0;
    datafile='';
    w0=0;
    w1=0;
    currentc=0;
    for s in lines:
        s=s.replace('%','').rstrip('\n').replace('\r','');
        if s.startswith('./svm-train'):
            s=s.split(' ');
            for i in range(len(s)):
                if s[i]=='-c':
                    currentc=s[i+1];
            
            datafile=s[len(s)-1];
            
            
        elif s.startswith('Cross Validation Accuracy'):
            s=s.split(' = ')[1];
            
            s=s.replace('%','');
            
            s=float(s);
            if s>bestresult:
                bestresult=s;
                bestc=currentc;
    fout.write('echo Best Score=%s Best c=%s\n'%(bestresult,bestc));                
    fout.write('echo Linear Classifier NotWeighted >>linear_notweighted.out\n');
    powc=int(np.log2(float(bestc)));
    for i in range(powc-2,powc+3):
        fout.write('echo ./svm-train -t 0 -c %s -v 5  %s >>linear_notweighted.out\n'%(2**i,datafile));
        fout.write('./svm-train -t 0 -c %s -v 5  %s >>linear_notweighted.out\n'%(2**i,datafile));
    
    #fout.write('./svm-train -t 0 -c %s %s %s.linear_noweight\n'%(bestc,datafile,datafile));


fout.close();
    

#%%