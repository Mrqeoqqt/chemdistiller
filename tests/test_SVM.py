# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 15:55:57 2016

@author: Dr. Ivan Laponogov
"""

import sys;
import numpy as np;
import scipy.sparse as sp;

if __name__=='__main__':

    
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

import os;
#dir_path = os.path.dirname(os.path.realpath(__file__));
chemdistiller_path = os.path.abspath(os.path.join(os.getcwd(), ".."));
#print('Expected chemdistiller path: %s'%chemdistiller_path);
sys.path.append(chemdistiller_path);

import time;

from chemdistiller.scorers.svm import predict_fingerprints_radial, predict_fingerprints_linear;

kernel=0;

if __name__=='__main__':
        testfname='e:/Imperial/TestDB/SVM_Train/0/Positive/0.test';
        testblock=[];
        reffpt=[];
        with open(testfname,'r') as finp:
            for s in finp:
                s=s.rstrip('\n');
                s=s.split(' ');
                reffpt.append(int(s.pop(0)));
                subs=[];
                for ss in s:
                    if ss!='':
                        ss=ss.split(':');
                        subs.append((int(ss[0]),float(ss[1])));
                testblock.append(subs);
        
        data=[];
        row=[];
        col=[];
                    
        for i in range(len(testblock)):
            subblock=testblock[i];
            for j in range(len(subblock)):
                row.append(i);
                col.append(subblock[j][0]);
                data.append(subblock[j][1]);
    
        row=np.array(row,dtype=np.int32);
        col=np.array(col,dtype=np.int32);
        data=np.array(data,dtype=np.float32);
    
        testvector=sp.csc_matrix((data, (row, col)), shape=(len(testblock), 20001));
        print(testvector.shape);
    
        if kernel==0:
            SVM_fname='e:/Imperial/TestDB/SVM_Train/Ported_csc/linear_noweight/0/1/w.dat';
            print(SVM_fname);
            
            curtime=time.time();
            result=predict_fingerprints_linear(SVM_fname, testvector);
            dtime=time.time()-curtime;

            result=result[:,0];
            print(result);
            print(reffpt);
        
            cc=0;    
            for i in range(len(reffpt)):
                print(reffpt[i],result[i]);
                if reffpt[i]==result[i]:
                    cc+=1;
            print('Precision: %s'%(float(cc)*100/len(reffpt)));
            print('Time: %s:%s:%s'%(int(dtime)//3600,
                                int(dtime)%3600//60,
                                int(dtime)%60));
        
        else:
            SVM_fname='e:/Imperial/TestDB/SVM_Train/Ported_csc/radial_noweight/0/1/0.svm';
            print(SVM_fname);
            curtime=time.time();
            result=predict_fingerprints_radial(SVM_fname, testvector, None);
            dtime=time.time()-curtime;
            print(result);
            print(reffpt);
            
        
            cc=0;    
            for i in range(len(reffpt)):
                print(reffpt[i],result[i]);
                if reffpt[i]==result[i]:
                    cc+=1;
            print('Precision: %s'%(float(cc)*100/len(reffpt)));
            print('Time: %s:%s:%s'%(int(dtime)//3600,
                                int(dtime)%3600//60,
                                int(dtime)%60));
