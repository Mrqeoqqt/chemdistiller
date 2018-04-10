# -*- coding: utf-8 -*-
"""
Created on Thu Mar 02 15:21:45 2017

@author: ilaponog
"""

import numpy as np;
import scipy.sparse as sp;
import os;
import time;
from glob import glob;
import h5py;

reftime=0;

def tic():
    global reftime;
    reftime=time.time();
    
def toc():
    global reftime;
    dtime=time.time()-reftime;
    print('%s:%s:%s'%(int(dtime)//3600,
                        int(dtime)%3600//60,
                        int(dtime)%60));


def load_svm_data(infname):
    print('Importing Test Data...');
    tic();
    with open(infname, 'r') as finp:
        v1=[];
        sparse1=[];
        for s in finp:
            s=s.rstrip('\n');
            s=s.split(' ');
            v1.append(float(s[0]));
            sublist=[];
            for i in range(1,len(s)):
                if s[i]!='':
                    ss=s[i].split(':');
                    sublist.append((int(ss[0]),float(ss[1])));
            sparse1.append(sublist);
            
        v1=np.array(v1, dtype=np.uint8);
    
        data=[];
        row=[];
        col=[];
            
        for i in range(len(sparse1)):
            subblock=sparse1[i];
            for j in range(len(subblock)):
                row.append(i);
                col.append(subblock[j][0]);
                data.append(subblock[j][1]);
        row=np.array(row,dtype=np.int32);
        col=np.array(col,dtype=np.int32);
        data=np.array(data,dtype=np.float32);
        sparse1=sp.csc_matrix((data, (row, col)));
        toc();
        return (sparse1, v1)

class SVM:
    def __init__(self, modelname, maskname):
        if os.path.isfile(modelname):
            print('Importing SVM model...')
            self._import_libSVM_model(modelname);
        else:
            raise IOError('File %s not found!'%(modelname));

        if os.path.isfile(maskname):
            print('Importing FPT mask...')
            self._import_mask(maskname);
        else:
            raise IOError('File %s not found!'%(maskname));

        
        print('Finished Importing SVM model.');
        
            
            
    def _import_mask(self, infname):
        with open(infname,'r') as finp:
            mask=[];
            for s in finp:
                mask.append(int(s));
        self.mask=np.array(mask, dtype=np.uint8)>0;
        

    def _import_masked_mean(self, infname):
        with open(infname,'r') as finp:
            masked_mean=[];
            for s in finp:
                masked_mean.append(float(s));
        self.masked_mean=np.array(masked_mean, dtype=np.float32);

    def _import_masked_std(self, infname):
        with open(infname,'r') as finp:
            masked_std=[];
            for s in finp:
                masked_std.append(float(s));
        self.masked_std=np.array(masked_std, dtype=np.float32);
    
    def _import_libSVM_model(self, infname):
        print('Importing libSVM model');
        tic();
        with open(infname,'r') as finp:
            reading_sv=False;
            self.count1=np.int32(0);
            self.count2=np.int32(0);
            v1=[];
            sparse1=[];
            
            self.use_probabilities=False;
            for s in finp:
                s=s.rstrip('\n');
                if reading_sv:
                    s=s.split(' ');
                    v1.append(float(s[0]));
                    sublist=[];
                    for i in range(1,len(s)):
                        if s[i]!='':
                            ss=s[i].split(':');
                            sublist.append((int(ss[0]),float(ss[1])));
                    sparse1.append(sublist);
                else:            
                    if s.startswith('svm_type'):
                        if s.split(' ')[1]!='c_svc':
                            raise TypeError('Wrong SVM type!');
                    elif s.startswith('kernel_type'):
                        s=s.split(' ')[1];
                        if s=='rbf':
                            self.kernel=1;
                        elif s=='linear':
                            self.kernel=0;
                        else:
                            raise TypeError('Wrong SVM kernel!');
                    elif s.startswith('gamma'):
                        self.gamma=np.float32(s.split(' ')[1]);
                    elif s.startswith('nr_class'):
                        s=s.split(' ')[1];
                        if s!='2':
                            raise TypeError('Wrong SVM number of classes!');
                    elif s.startswith('rho'):
                        self.rho=np.float64(s.split(' ')[1]);
                    elif s.startswith('label'):
                        s=s.split(' ');
                        self.label1=np.uint8(s[1]);
                        self.label2=np.uint8(s[2]);
                    elif s.startswith('probA'):
                        s=s.split(' ');
                        self.probA=np.float32(s[1]);
                        self.use_probabilities=True;
                    elif s.startswith('probB'):
                        s=s.split(' ');
                        self.probB=np.float32(s[1]);
                        self.use_probabilities=True;
                    elif s.startswith('nr_sv'):
                        s=s.split(' ');
                        self.count1=np.int32(s[1]);
                        self.count2=np.int32(s[2]);
                    elif s.startswith('SV'):
                        reading_sv=True;
        
        if len(v1)!=self.count1+self.count2:
            raise TypeError('Wrong SVM format!');
        
        v1=np.array(v1,dtype=np.float32);
    
        if self.label1!=1:
            self.label2,self.label1=self.label1,self.label2;
            self.count1,self.count2=self.count2,self.count1;
            v1=v1*-1;
            self.rho=self.rho*-1;
        
        if self.kernel==0:
            data=[];
            row=[];
            col=[];
            
            for i in range(len(sparse1)):
                subblock=sparse1[i];
                for j in range(len(subblock)):
                    row.append(i);
                    col.append(subblock[j][0]);
                    data.append(subblock[j][1]);
            row=np.array(row,dtype=np.int32);
            col=np.array(col,dtype=np.int32);
            data=np.array(data,dtype=np.float32);
            sparse1=sp.csc_matrix((data, (row, col)));
            self.w=v1.dot(sparse1);        

        else:
            data=[];
            row=[];
            col=[];
            self.sv_coef=v1;            
            
            for i in range(len(sparse1)):
                subblock=sparse1[i];
                for j in range(len(subblock)):
                    row.append(i);
                    col.append(subblock[j][0]);
                    data.append(subblock[j][1]);
            row=np.array(row,dtype=np.int32);
            col=np.array(col,dtype=np.int32);
            data=np.array(data,dtype=np.float32);
            sparse1=sp.csc_matrix((data, (row, col)));
            tic();
            print('Sparse*Sparse');
            self.sp_transp=sparse1.transpose().toarray();
            self.sparse1dotsparse1=sparse1.dot(self.sp_transp).diagonal()*(-self.gamma);
            toc();
            print('pre-computed')
            
        print('done')
        toc();
    
    def predict(self, X_masked): #masked and corrected array expected
            
        if self.kernel==1:

            ntest=X_masked.shape[0];
            #print('dotmat...')
            #tic();
            #dotmat=X_masked.dot(self.sp_transp).toarray()*2*self.gamma;
            dotmat=X_masked.dot(self.sp_transp)*2*self.gamma;
            #print(dotmat.shape)
            #toc();
            #print('computed')
            #print('X*X')
            #tic();
            X_masked=X_masked.dot(X_masked.transpose()).diagonal()*(-self.gamma);
            #print(X_masked.shape)
            #toc();
            #print('computed')

            #print('first difference')
            #tic();
            s1=np.add(dotmat, X_masked.reshape(ntest,1))
            #toc();
            #print('computed')

            #print('second difference')
            #tic();
            s1=np.add(s1, self.sparse1dotsparse1);
            #toc();
            #print('computed')


            #print('exponent')
            #tic();
            s1=np.exp(s1);
            #toc();
            #print('computed')

            #print('multiply')
            #tic();
            s1=np.multiply(s1,self.sv_coef);
            #toc();
            #print('computed')

            #print('sum')
            #tic();
            s1=np.sum(s1,1);
            #toc();
            #print('computed')
            
            predicted=np.subtract(s1, self.rho);
            
        else:
            #print('Linear prediction');
            #tic();
            predicted=np.subtract(X_masked.dot(self.w.transpose()).toarray(), self.rho);
            #toc();
            #print('computed');
            
        if self.use_probabilities:
            #print('probability')
            #tic();
            predicted=predicted*self.probA+self.probB;
            
            mask = predicted<0;
            
            mdata = np.ma.masked_array(predicted, mask);
            
            mdata1=np.exp(-mdata);
            mdata1=np.divide(mdata1,(1.0+mdata1));
            
            mask = predicted>=0;
            mdata = np.ma.masked_array(predicted, mask);
            
            mdata2=np.exp(mdata);
            mdata2=np.divide(1.0,(1+mdata2));

            predicted=np.add(np.ma.filled(mdata1,0.0), np.ma.filled(mdata2, 0.0))
            #toc();
            #print('computed')

        else:
            #print('sign');
            #tic();
            predicted=np.array((np.sign(predicted)+1)/2, dtype=np.uint8);
            #toc();
            #print('computed');
            
        return predicted;

#%%
#smX=X[0:100,:]

#print(smX);

#mean=svm.masked_mean;

#std=svm.masked_std;

#print(mean)
#print(std)

#X_masked=X[:, self.mask];
#X_masked=X;
#X_masked=np.divide(X_masked-self.masked_mean), self.masked_std);
        

#%% 

dbpath='g:/DB_HDF5_MetaboliteLikeness_twoclass_best';
outfile='g:/DB_HDF5_MetaboliteLikeness_twoclass_best/likeness.log';
   
    
if __name__=='__main__':
    modelname='g:/MetaboliteLikeness/results3/best/model_4.model'
    maskname='g:/MetaboliteLikeness/results3_home3/UniqueFPTsFullInChi_fixed_subsampled_notcorrected_randsel.mask'
    svm = SVM(modelname, maskname);
    
#%%
if __name__=='__main__':
    onlyfiles = [ff for ff in os.listdir(dbpath) if os.path.isfile(os.path.join(dbpath, ff))];
    
    
    for dbfile in onlyfiles:
        if dbfile.lower().endswith('.h5'):
            print(dbfile);
            ccc=0;
            ttt=0;
            #subf=;
            
            for subf in ['/Neutral', '/Negative', '/Positive']:
                with h5py.File(os.path.join(dbpath,dbfile)) as hdf5_container:
                    #ascii_dataset=hdf5_container['%s%s'%(subf,'/ChemInfo/ASCII')];
                    #inchi_dataset=hdf5_container['%s%s'%(subf,'/ChemInfo/InChi')];
                    fpts_dataset=hdf5_container['%s%s'%(subf,'/FingerPrints/FPTArray')];
                    #mz_dataset=hdf5_container['%s%s'%(subf,'/MZ')];
                    total_count=len(fpts_dataset);
                    ttt+=total_count;
                    print('Total: %s'%total_count);
                    
                    if '%s%s'%(subf,'/Likeness/MetaLikeness') in hdf5_container:
                        metalikeness_dataset=hdf5_container['%s%s'%(subf,'/Likeness/MetaLikeness')];
                    else:
                        metalikeness_dataset=hdf5_container.create_dataset('%s%s'%(subf,'/Likeness/MetaLikeness'), (total_count, ), chunks=(10000, ), maxshape=(None, ), compression="gzip", compression_opts=4, dtype=np.int8);
                        metalikeness_dataset[:]=0;
                    ranges=[];
                        
                    for i in range(total_count//10000):
                        ranges.append((i*10000,(i+1)*10000));
                    ranges.append(((total_count//10000)*10000, total_count));
                    
                    for r in ranges:
                        print(r)
                        fpts=np.unpackbits(fpts_dataset[r[0]:r[1],:],axis=1);
                        fpts=fpts[:, svm.mask].astype(np.float32);
                        #fpts=np.divide(np.subtract(fpts,svm.masked_mean), svm.masked_std);
                        #print(fpts.shape)
                        #likeness=svm.predict(sp.csc_matrix(fpts));
                        l=svm.predict(fpts)[:];
                        metalikeness_dataset[r[0]:r[1]]=l
                        ccc+=np.sum(l>0);
                        
            fout=open(outfile,'a');            
            rate=float(ccc)/ttt*100;
            fout.write('%s\t%s\n'%(dbfile, rate));
            print('%s\t%s\n'%(dbfile, rate));
            fout.close();
            