# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 17:23:38 2016

@author: Dr. Ivan Laponogov
"""


import sys;
import os;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);



import numpy as np;
import scipy.sparse as sp;
from chemdistiller.settings import test_as_single_process;

if not test_as_single_process:
    import multiprocessing;


def predict_fingerprints_linear(SVM_fname, testvector):
    if not os.path.isfile(SVM_fname):
        return None;
    with open(SVM_fname,'rb') as finp:
        kernel=np.fromfile(finp, dtype=np.uint8, count=1)[0];
        if kernel!=0:
            return None;
        
        count=np.fromfile(finp, dtype=np.int32, count=1)[0];
        rho=np.fromfile(finp, dtype=np.float64, count=count);
        datalen=np.fromfile(finp, dtype=np.int32, count=1)[0];
        data=np.fromfile(finp, dtype=np.float32, count=datalen);
        indiceslen=np.fromfile(finp, dtype=np.int32, count=1)[0];
        indices=np.fromfile(finp, dtype=np.int32, count=indiceslen);
        indptrlen=np.fromfile(finp, dtype=np.int32, count=1)[0];
        indptr=np.fromfile(finp, dtype=np.int32, count=indptrlen);
        w=sp.csc_matrix((data, indices, indptr), shape=(count, 20001));
    
    

    fpt_predicted=np.array(\
            (np.sign(\
                np.subtract(\
                    testvector.dot(w.transpose()).toarray()\
                ,rho)\
             )+1)/2,\
        dtype=np.uint8);
    
    
    return fpt_predicted;


def predict_fingerprints_radial(arg):
    SVM_fname, testvector, progress_report_queue=arg;
    if not os.path.isfile(SVM_fname):
        if not (progress_report_queue is None):
            progress_report_queue.put(1);   
        return None;
        
    with open(SVM_fname,'rb') as finp:

        kernel=np.fromfile(finp, dtype=np.uint8, count=1)[0];

        if kernel!=1:
            if not (progress_report_queue is None):
                progress_report_queue.put(1);   
            return None;
            
        rho=np.fromfile(finp, dtype=np.float64, count=1)[0];
        count1=np.fromfile(finp, dtype=np.int32, count=1)[0];
        count2=np.fromfile(finp, dtype=np.int32, count=1)[0];
        label1=np.fromfile(finp, dtype=np.uint8, count=1)[0];
        label2=np.fromfile(finp, dtype=np.uint8, count=1)[0];
        gamma=np.fromfile(finp, dtype=np.float32, count=1)[0];
        sv_coeff=np.fromfile(finp, dtype=np.float32, count=count1+count2);
        datalen=np.fromfile(finp, dtype=np.int32, count=1)[0];
        data=np.fromfile(finp, dtype=np.float32, count=datalen);
        indiceslen=np.fromfile(finp, dtype=np.int32, count=1)[0];
        indices=np.fromfile(finp, dtype=np.int32, count=indiceslen);
        indptrlen=np.fromfile(finp, dtype=np.int32, count=1)[0];
        indptr=np.fromfile(finp, dtype=np.int32, count=indptrlen);
        sparse1=sp.csc_matrix((data, indices, indptr), shape=(count1+count2, 20001));
    
    if label1!=1 or label2!=0:
        if not (progress_report_queue is None):
            progress_report_queue.put(1);
        return None;

    #Calculate SVM here

    ntest=testvector.shape[0];
    
    dotmat=testvector.dot(sparse1.transpose()).toarray()*2*gamma;
    
    testvector=testvector.dot(testvector.transpose()).diagonal()*(-gamma);

    sparse1=sparse1.dot(sparse1.transpose()).diagonal()*(-gamma);

    predicted_fpt=np.array(\
        (np.sign(\
            np.sum(\
                np.multiply(\
                    np.exp(\
                        np.add(\
                            np.add(dotmat, testvector.reshape(ntest,1))\
                        ,sparse1)\
                    )\
                ,sv_coeff)\
            ,1)\
        -rho)+1\
        )/2, dtype=np.uint8);

    if not (progress_report_queue is None):
        progress_report_queue.put(1);

    return predicted_fpt;



class SVM:
    supported_adducts=set(['[M+H]+','[M-H]-']);
    supported_isotopes=set([0]);    

    def __init__(self, fingerprint_model_path, settings, nproc=1):
        self.fingerprint_model_path=fingerprint_model_path;
        self.nproc=nproc;
        self.settings=settings;

    def process_annotations(self, annotations):
        if annotations:
            self.mode=annotations[0].adduct.mode;
            self.adduct=annotations[0].adduct;
            self.isotope=annotations[0].isotope;
            if self.__configure_svm():

                data=[];
                row=[];
                col=[];
                
                for i in range(len(annotations)):
                    annotation=annotations[i];
                    annotation.fpt_mask=self.fpt_mask;
                    annotation.predicted_fpt=np.zeros((11416,), dtype=np.uint8);
                    subblock=annotation.parent_peak.fpt_sparse_vector;
                    for j in range(len(subblock)):
                        row.append(i);
                        col.append(subblock[j][0]);
                        data.append(subblock[j][1]);

                row=np.array(row,dtype=np.int32);
                col=np.array(col,dtype=np.int32);
                data=np.array(data,dtype=np.float32);

                testvector=sp.csc_matrix((data, (row, col)), shape=(len(annotations), 20001));
                count=len(self.fpt_indexes);
                
                if self.svm_type==0:
                    if self.mode==1:                    
                        print('Starting Linear SVM fingerprint prediction for positive mode...');
                    else:
                        print('Starting Linear SVM fingerprint prediction for negative mode...');

                    fname=os.path.join(self.fingerprint_model_path,str(self.mode),'w.dat');
                    fpt_predicted=predict_fingerprints_linear(fname, testvector);

                    if fpt_predicted is None:
                        print('Failed prediction of Fingerprints');
                    else:      
                        for j in range(len(annotations)):
                            pfpt=annotations[j].predicted_fpt;
                            for i in range(count):
                                pfpt[self.fpt_indexes[i]]=fpt_predicted[j,self.pre_indexes[i]];
                    
                elif self.svm_type==1:
                    
                    if test_as_single_process:
                        if self.mode==1:                    
                            print('Starting Radial SVM fingerprint prediction for positive mode...');
                        else:
                            print('Starting Radial SVM fingerprint prediction for negative mode...');
                            
                        
                        for i in range(count):
                            fpt_index=self.fpt_indexes[i];
                            print('\rFingerPrit SVM %s of %s '%(i+1,count)),
                            fname=os.path.join(self.fingerprint_model_path, str(self.mode),'%s.svm'%fpt_index);
                            fpt_predicted=predict_fingerprints_radial((fname, testvector, None));
                            if fpt_predicted is None:
                                print('Failed prediction of Fingerprint %s'%fpt_index);
                            else:
                                for j in range(len(annotations)):
                                    annotations[j].predicted_fpt[fpt_index]=fpt_predicted[j];
                    else:

                        #progress_display=ProgressDisplay('Predicting FingerPrints using Radial SVMs...', self.nproc);
                        pool=multiprocessing.Pool(self.nproc);
                        
                        args=[];
                        m = multiprocessing.Manager();
                        progress_report_queue = m.Queue();
                        
                        for i in range(count):
                            fpt_index=self.fpt_indexes[i];
                            fname=os.path.join(self.fingerprint_model_path,str(self.mode),'%s.svm'%fpt_index);
                            args.append((fname, testvector, progress_report_queue));
                        print('Starting Pool (%s) Processing'%self.nproc);
                        poolresult=pool.map_async(predict_fingerprints_radial, args);
                        progress=0;
                        while not poolresult.ready():
                            if not progress_report_queue.empty():
                                progress+=progress_report_queue.get();
                                print('\rProgress: %s of %s'%(progress, count));
                        
                        result=poolresult.get();
                        pool.close();
                        pool.join();
                        
                        for i in range(count):
                            fpt_index=self.fpt_indexes[i];
                            fpt_predicted=result[i];
                            if fpt_predicted is None:
                                print('Failed prediction of Fingerprint %s'%fpt_index);
                            else:
                                for j in range(len(annotations)):
                                    annotations[j].predicted_fpt[fpt_index]=fpt_predicted[j];
                #Finished prediction                
                for annotation in annotations:
                    annotation.predicted_fpt=np.packbits(annotation.predicted_fpt);

                    
        #Consider Negative/Positive mode, annotation adducts and isotopes                    
        
    def __configure_svm(self):

        fname=os.path.join(self.fingerprint_model_path,str(self.mode),'fpt_mask.default');
        if not os.path.isfile(fname):
            return False;
            
        with open(fname,'rb') as finp:
            self.fpt_mask=np.fromfile(finp, dtype=np.uint8);

        fname=os.path.join(self.fingerprint_model_path,str(self.mode),'fpt_stats.dat');
        if not os.path.isfile(fname):
            return False;
            
        with open(fname,'rb') as finp:
            self.svm_type=np.fromfile(finp, dtype=np.uint8, count=1)[0];            
            self.fpt_stats=np.fromfile(finp, dtype=np.float32);
            self.fpt_stats.reshape((-1,10));  # occupancy, prediction score, w0, w1, c, g, nsvm0, nsvm1, totalsvm ntrainingset

        self.fpt_indexes=[];
        self.pre_indexes=[];
        mask=np.unpackbits(self.fpt_mask);
        #self.raw_mask=mask;
        ii=-1;
        for i in range(len(mask)):
            if mask[i]>0:
                ii+=1;
                passed=True;
                if self.settings:
                    for j in range(10):
                        if self.fpt_stats[ii,j]<self.settings[j][0] or self.fpt_stats[ii,j]>self.settings[j][1]:
                            passed=False;
                            break;
                    if abs(self.fpt_stats[ii,1]-50.0)-abs(self.fpt_stats[ii,0]*100.0-50.0)<self.settings[10][0] or abs(self.fpt_stats[ii,1]-50.0)-abs(self.fpt_stats[ii,0]*100.0-50.0)>self.settings[10][1] :
                        passed=False;
                if passed:
                    self.fpt_indexes.append(i);
                    self.pre_indexes.append(ii);
                else:
                    mask[i]=0;
        self.fpt_mask=np.packbits(mask);
        
        return True;
        
            
if __name__ == '__main__':
    if not test_as_single_process:
        multiprocessing.freeze_support();

        
        