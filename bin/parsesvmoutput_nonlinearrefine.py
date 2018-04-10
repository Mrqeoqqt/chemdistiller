# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 15:30:02 2016

@author: Dr. Ivan Laponogov
"""
import os;
import numpy as np;
fout=open('do_models_nonlinearrefine.bat','w');
fout.write('#!/bin/bash\n');
if os.path.isfile('radial_weighted.out'):
    fout.write('echo radial weighted \n');
    
    finp=open('radial_weighted.out','r');
    lines=finp.readlines();
    finp.close();

    bestresult=0;
    bestc=0;
    bestg=0;
    datafile='';
    w0=0;
    w1=0;
    currentc=0;
    currentg=0;
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
                elif s[i]=='-g':
                    currentg=s[i+1];                    
            datafile=s[len(s)-1];
            
        elif s.startswith('Cross Validation Accuracy'):
            s=s.split(' = ')[1];
            
            s=s.replace('%','');
            
            s=float(s);
            if s>bestresult:
                bestresult=s;
                bestc=currentc;
                bestg=currentg;
    fout.write('echo Radial Basis Classifier Weighted >>radial_weighted.out\n');
    fout.write('echo Best Score=%s Best c=%s Best g=%s\n'%(bestresult,bestc,bestg));                
    powc=int(np.log2(float(bestc)));
    powg=int(np.log2(float(bestg)));
    for i in range(powc-2,powc+3):
        for j in range(powg-2,powg+3):
                fout.write('echo ./svm-train -t 2 -c %s -g %s -v 5 -w0 %s -w1 %s %s >>radial_weighted.out\n'%(2**i,2**j,w0,w1,datafile));
                fout.write('./svm-train -t 2 -c %s -g %s -v 5 -w0 %s -w1 %s %s >>radial_weighted.out\n'%(2**i,2**j,w0,w1,datafile));
    
    #fout.write('./svm-train -t 2 -c %s -g %s -w0 %s -w1 %s %s %s.radial_weight\n'%(bestc,bestg,w0,w1,datafile,datafile));


if os.path.isfile('radial_notweighted.out'):
    fout.write('echo radial not weighted \n');
    
    finp=open('radial_notweighted.out','r');
    lines=finp.readlines();
    finp.close();

    bestresult=0;
    bestc=0;
    bestg=0;
    datafile='';
    w0=0;
    w1=0;
    currentc=0;
    currentg=0;
    for s in lines:
        s=s.replace('%','').rstrip('\n').replace('\r','');
        if s.startswith('./svm-train'):
            s=s.split(' ');
            for i in range(len(s)):
                if  s[i]=='-c':
                    currentc=s[i+1];                    
                elif s[i]=='-g':
                    currentg=s[i+1];                    
            datafile=s[len(s)-1];
            
        elif s.startswith('Cross Validation Accuracy'):
            s=s.split(' = ')[1];
            
            s=s.replace('%','');
            
            s=float(s);
            if s>bestresult:
                bestresult=s;
                bestc=currentc;
                bestg=currentg;
    powc=int(np.log2(float(bestc)));
    powg=int(np.log2(float(bestg)));
    fout.write('echo Radial Basis Classifier NotWeighted >>radial_notweighted.out\n');                
    fout.write('echo Best Score=%s Best c=%s Best g=%s\n'%(bestresult,bestc,bestg));                
    for i in range(powc-2,powc+3):
        for j in range(powg-2,powg+3):
                fout.write('echo ./svm-train -t 2 -c %s -g %s -v 5 %s >>radial_notweighted.out\n'%(2**i,2**j,datafile));
                fout.write('./svm-train -t 2 -c %s -g %s -v 5 %s >>radial_notweighted.out\n'%(2**i,2**j,datafile));
    
    
    #fout.write('./svm-train -t 2 -c %s -g %s %s %s.radial_weight\n'%(bestc,bestg,datafile,datafile));



fout.close();
    

#%%