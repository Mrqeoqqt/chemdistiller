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
        sparse1=sparse1.toarray();
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
    
    def _import_libSVM_model(self, infname):
        print('Importing libSVM model');
        tic();
        with open(infname,'r') as finp:
            reading_sv=False;
            self.count1=np.int32(0);
            v1=[];
            sparse1=[];
            
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
                        if s.split(' ')[1]!='svdd':
                            raise TypeError('Wrong SVM type!');
                    elif s.startswith('kernel_type'):
                        s=s.split(' ')[1];
                        if s=='rbf':
                            self.kernel=1;
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
                    elif s.startswith('total_sv'):
                        s=s.split(' ');
                        self.count1=np.int32(s[1]);
                        
                    elif s.startswith('SV'):
                        reading_sv=True;
        
        if len(v1)!=self.count1:
            raise TypeError('Wrong SVM format!');
        else:
            v1=np.array(v1,dtype=np.float32);
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
            #self.a=sparse1.toarray();
            tic();
            print('Sparse*Sparse');
            self.sp_transp=sparse1.transpose().toarray();
            self.sparse1dotsparse1=sparse1.dot(self.sp_transp).diagonal()*(-self.gamma);
            toc();
            print('pre-computed')
            
        print('done')
        toc();
    
    def predict(self, X_masked): #masked and corrected array expected
        ntest=X_masked.shape[0];
        dotmat=X_masked.dot(self.sp_transp)*2*self.gamma;
        X_masked=X_masked.dot(X_masked.transpose()).diagonal()*(-self.gamma);
        s1=np.add(dotmat, X_masked.reshape(ntest,1))
        s1=np.add(s1, self.sparse1dotsparse1);
        s1=np.exp(s1);
        s1=np.multiply(s1,self.sv_coef);
        s1=np.multiply(s1,-2.0);
        s1=np.sum(s1,1);
        s1=np.add(s1,1+2*self.rho);
        s1=np.sign(-s1).astype(dtype=np.int8);    
        return s1;


dbpath='g:/NoMetabolite/';
outfile='g:/NoMetabolite/Total_prediction.txt';   
#%%

worklist=[
["Metabolites",	99	,	0.001	,	0.00012207	,	99.179	],
#["Metabolites",	36	,	0.2	,	0.000244141	,	99.8678	],
#["Metabolites",	4	,	0.5	,	0.000244141	,	99.8517	],
#["Metabolites",	92	,	0.002	,	0.000244141	,	99.5667	],
#["Metabolites",	82	,	0.005	,	6.10E-05	,	99.7548	],
#["Metabolites",	76	,	0.01	,	0.000244141	,	99.8241	],
#["Metabolites",	69	,	0.025	,	0.000488281	,	99.8222	],
#["Metabolites",	60	,	0.05	,	0.000244141	,	99.8395	],
#["Metabolites",	52	,	0.1	,	0.000244141	,	99.8402	],
#["Metabolites",	44	,	0.15	,	0.000244141	,	99.8357	],

#["Metabolites",	28	,	0.25	,	0.000244141	,	99.8312	],
#["Metabolites",	20	,	0.3	,	0.000244141	,	99.853	],
#["Metabolites",	12	,	0.4	,	0.000244141	,	99.8183	],

["BMDB",	81	,	0.005	,	3.05E-05	,	96.2599	],
#["BMDB",	77	,	0.01	,	0.000488281	,	96.5407	],
#["BMDB",	69	,	0.025	,	0.000488281	,	97.0896	],
#["BMDB",	60	,	0.05	,	0.000244141	,	97.179	],
#["BMDB",	53	,	0.1	,	0.000488281	,	96.8343	],
#["BMDB",	44	,	0.15	,	0.000244141	,	96.796	],
#["BMDB",	36	,	0.2	,	0.000244141	,	97.0513	],
#["BMDB",	28	,	0.25	,	0.000244141	,	96.8471	],
#["BMDB",	21	,	0.3	,	0.000488281	,	96.7577	],
#["BMDB",	13	,	0.4	,	0.000488281	,	97.3321	],
#["BMDB",	5	,	0.5	,	0.000488281	,	97.2173	],
["ChEBI",	98	,	0.001	,	6.10E-05	,	98.4807	],
#["ChEBI",	92	,	0.002	,	0.000244141	,	99.1535	],
#["ChEBI",	85	,	0.005	,	0.000488281	,	99.5327	],
#["ChEBI",	76	,	0.01	,	0.000244141	,	99.6318	],
#["ChEBI",	68	,	0.025	,	0.000244141	,	99.6954	],
#["ChEBI",	60	,	0.05	,	0.000244141	,	99.6942	],
#["ChEBI",	52	,	0.1	,	0.000244141	,	99.7284	],
#["ChEBI",	44	,	0.15	,	0.000244141	,	99.6979	],
#["ChEBI",	36	,	0.2	,	0.000244141	,	99.6807	],
#["ChEBI",	28	,	0.25	,	0.000244141	,	99.6905	],
#["ChEBI",	21	,	0.3	,	0.000488281	,	99.6477	],
#["ChEBI",	12	,	0.4	,	0.000244141	,	99.7138	],
#["ChEBI",	4	,	0.5	,	0.000244141	,	99.6991	],
["DrugBank",	73	,	0.01	,	3.05E-05	,	97.2337	],
#["DrugBank",	68	,	0.025	,	0.000244141	,	98.1178	],
#["DrugBank",	59	,	0.05	,	0.00012207	,	98.4457	],
#["DrugBank",	50	,	0.1	,	6.10E-05	,	98.1606	],
#["DrugBank",	43	,	0.15	,	0.00012207	,	98.4885	],
#["DrugBank",	34	,	0.2	,	6.10E-05	,	98.5598	],
#["DrugBank",	27	,	0.25	,	0.00012207	,	98.7452	],
#["DrugBank",	20	,	0.3	,	0.000244141	,	98.3459	],
#["DrugBank",	10	,	0.4	,	6.10E-05	,	98.4743	],
#["DrugBank",	3	,	0.5	,	0.00012207	,	98.4457	],
["ECMDB",	81	,	0.005	,	3.05E-05	,	93.6712	],
#["ECMDB",	73	,	0.01	,	3.05E-05	,	95.1998	],
#["ECMDB",	66	,	0.025	,	6.10E-05	,	95.0657	],
#["ECMDB",	62	,	0.05	,	0.000976563	,	96.1384	],
#["ECMDB",	50	,	0.1	,	6.10E-05	,	93.4567	],
#["ECMDB",	46	,	0.15	,	0.000976563	,	93.3226	],
#["ECMDB",	37	,	0.2	,	0.000488281	,	93.0813	],
#["ECMDB",	28	,	0.25	,	0.000244141	,	94.8243	],
#["ECMDB",	21	,	0.3	,	0.000488281	,	92.384	],
#["ECMDB",	11	,	0.4	,	0.00012207	,	92.7863	],
#["ECMDB",	5	,	0.5	,	0.000488281	,	93.3226	],
["FooDB",	85	,	0.005	,	0.000488281	,	98.269	],
#["FooDB",	77	,	0.01	,	0.000488281	,	98.6995	],
#["FooDB",	69	,	0.025	,	0.000488281	,	98.7786	],
#["FooDB",	60	,	0.05	,	0.000244141	,	98.4535	],
#["FooDB",	53	,	0.1	,	0.000488281	,	98.8533	],
#["FooDB",	45	,	0.15	,	0.000488281	,	98.9544	],
#["FooDB",	37	,	0.2	,	0.000488281	,	98.95	],
#["FooDB",	28	,	0.25	,	0.000244141	,	98.9368	],
#["FooDB",	20	,	0.3	,	0.000244141	,	98.8138	],
#["FooDB",	12	,	0.4	,	0.000244141	,	98.8094	],
#["FooDB",	4	,	0.5	,	0.000244141	,	98.6468	],
["HMDB",	92	,	0.002	,	0.000244141	,	98.5631	],
#["HMDB",	81	,	0.005	,	3.05E-05	,	99.5642	],
#["HMDB",	73	,	0.01	,	3.05E-05	,	99.6504	],
#["HMDB",	69	,	0.025	,	0.000488281	,	99.1427	],
#["HMDB",	61	,	0.05	,	0.000488281	,	99.2768	],
#["HMDB",	54	,	0.1	,	0.000976563	,	99.2337	],
#["HMDB",	45	,	0.15	,	0.000488281	,	99.2696	],
#["HMDB",	38	,	0.2	,	0.000976563	,	99.2504	],
#["HMDB",	29	,	0.25	,	0.000488281	,	99.2385	],
#["HMDB",	21	,	0.3	,	0.000488281	,	99.0804	],
#["HMDB",	14	,	0.4	,	0.000976563	,	99.2002	],
#["HMDB",	6	,	0.5	,	0.000976563	,	99.3271	],
["LipidMaps",	99	,	0.001	,	0.00012207	,	97.1336	],
#["LipidMaps",	92	,	0.002	,	0.000244141	,	98.2996	],
#["LipidMaps",	85	,	0.005	,	0.000488281	,	98.9633	],
#["LipidMaps",	78	,	0.01	,	0.000976563	,	99.2666	],
#["LipidMaps",	69	,	0.025	,	0.000488281	,	99.5699	],
#["LipidMaps",	60	,	0.05	,	0.000244141	,	99.458	],
#["LipidMaps",	51	,	0.1	,	0.00012207	,	99.6047	],
#["LipidMaps",	44	,	0.15	,	0.000244141	,	99.5227	],
#["LipidMaps",	36	,	0.2	,	0.000244141	,	99.6122	],
#["LipidMaps",	28	,	0.25	,	0.000244141	,	99.5749	],
#["LipidMaps",	20	,	0.3	,	0.000244141	,	99.458	],
#["LipidMaps",	12	,	0.4	,	0.000244141	,	99.5152	],
#["LipidMaps",	4	,	0.5	,	0.000244141	,	99.6942	],
["MassBank",	85	,	0.005	,	0.000488281	,	97.4618	],
#["MassBank",	73	,	0.01	,	3.05E-05	,	98.7098	],
#["MassBank",	68	,	0.025	,	0.000244141	,	98.8448	],
#["MassBank",	59	,	0.05	,	0.00012207	,	98.9712	],
#["MassBank",	53	,	0.1	,	0.000488281	,	98.8616	],
#["MassBank",	43	,	0.15	,	0.00012207	,	99.0134	],
#["MassBank",	36	,	0.2	,	0.000244141	,	98.9375	],
#["MassBank",	27	,	0.25	,	0.00012207	,	99.0218	],
#["MassBank",	20	,	0.3	,	0.000244141	,	98.9628	],
#["MassBank",	12	,	0.4	,	0.000244141	,	99.0471	],
#["MassBank",	3	,	0.5	,	0.00012207	,	99.064	],
["PlantCyc",	77	,	0.01	,	0.000488281	,	96.0215	],
#["PlantCyc",	67	,	0.025	,	0.00012207	,	97.472	],
#["PlantCyc",	59	,	0.05	,	0.00012207	,	97.6792	],
#["PlantCyc",	53	,	0.1	,	0.000488281	,	97.2441	],
#["PlantCyc",	45	,	0.15	,	0.000488281	,	97.4513	],
#["PlantCyc",	36	,	0.2	,	0.000244141	,	97.9072	],
#["PlantCyc",	27	,	0.25	,	0.00012207	,	97.6378	],
#["PlantCyc",	20	,	0.3	,	0.000244141	,	97.5964	],
#["PlantCyc",	11	,	0.4	,	0.00012207	,	97.3891	],
#["PlantCyc",	3	,	0.5	,	0.00012207	,	97.7	],
["T3DB",	73	,	0.01	,	3.05E-05	,	96.7046	],
#["T3DB",	69	,	0.025	,	0.000488281	,	96.5548	],
#["T3DB",	60	,	0.05	,	0.000244141	,	97.154	],
#["T3DB",	50	,	0.1	,	6.10E-05	,	97.6633	],
#["T3DB",	44	,	0.15	,	0.000244141	,	96.7346	],
#["T3DB",	36	,	0.2	,	0.000244141	,	96.7046	],
#["T3DB",	28	,	0.25	,	0.000244141	,	97.1839	],
#["T3DB",	19	,	0.3	,	0.00012207	,	97.154	],
#["T3DB",	13	,	0.4	,	0.000488281	,	96.6747	],
#["T3DB",	5	,	0.5	,	0.000488281	,	97.0042	],
["YMDB",	72	,	0.01	,	1.53E-05	,	96.0922	]];
#["YMDB",	68	,	0.025	,	0.000244141	,	95.0902	],
#["YMDB",	61	,	0.05	,	0.000488281	,	95.6413	],
#["YMDB",	53	,	0.1	,	0.000488281	,	95.1403	],
#["YMDB",	45	,	0.15	,	0.000488281	,	95.4409	],
#["YMDB",	36	,	0.2	,	0.000244141	,	95.7415	],
#["YMDB",	29	,	0.25	,	0.000488281	,	95.5411	],
#["YMDB",	19	,	0.3	,	0.00012207	,	95.491	],
#["YMDB",	12	,	0.4	,	0.000244141	,	95.5411	],
#["YMDB",	5	,	0.5	,	0.000488281	,	95.5411	]];

    
#%%
if __name__=='__main__':
    onlyfiles = [ff for ff in os.listdir(dbpath) if os.path.isfile(os.path.join(dbpath, ff))];
    
    
    for opt in worklist:
        fout=open(outfile,'a');
        refDB=opt[0]
        coef=opt[2];
        model_index=opt[1];
        print(refDB,coef,model_index)
        
        inname='g:/MetaboliteLikeness/OneClass/results/%s/UniqueFPTsFullInChi_%s_like'%(refDB, refDB);
        maskname=inname+'.mask';
        
        modelname=inname+'.model_%s'%model_index;
        
        
        svm = SVM(modelname, maskname);
        
        for dbfile in onlyfiles:
            if dbfile.lower().endswith('.h5'):
                print(dbfile);
                #subf=;
                with h5py.File(os.path.join(dbpath,dbfile)) as hdf5_container:
                    rate=0;
                    ccc=0;
                    #nonmetafptdataset=inh5nonmeta['/FingerPrints/FPTArray'];
                    #nonmetamasschargedataset=inh5nonmeta['/FingerPrints/MassCharge'];
                    #nonmetauniqueiddataset=inh5nonmeta['/FingerPrints/UniqueID'];
                    #nonmetainchidataset=inh5nonmeta['/FingerPrints/InChi'];

                    for subf in ['']:
                            #try:
                            fpts_dataset=hdf5_container['%s%s'%(subf,'/FingerPrints/FPTArray')];
                            total_count=len(fpts_dataset);
                            ccc+=total_count;
                            print('Total: %s'%total_count);
                            
                            if '%s%s'%(subf,'/Likeness/%s_%s'%(refDB, coef)) in hdf5_container:
                                metalikeness_dataset=hdf5_container['%s%s'%(subf,'/Likeness/%s_%s'%(refDB, coef))];
                            else:
                                metalikeness_dataset=hdf5_container.create_dataset('%s%s'%(subf,'/Likeness/%s_%s'%(refDB, coef)), (total_count, ), chunks=(10000, ), maxshape=(None, ), compression="gzip", compression_opts=4, dtype=np.int8);
                                metalikeness_dataset[:]=0;
                            ranges=[];
                                
                            for i in range(total_count//10000):
                                ranges.append((i*10000,(i+1)*10000));
                            ranges.append(((total_count//10000)*10000, total_count));
                            
                            for r in ranges:
                                print(r)
                                fpts=fpts_dataset[r[0]:r[1],:];
                                print(fpts.shape)
                                fpts=fpts[:, svm.mask].astype(np.float32);
                                l=svm.predict(fpts)[:];
                                metalikeness_dataset[r[0]:r[1]]=l;
                                rate+=np.sum(l>=0);
                                
                                #except:
                                #    print('Error')
                    rate=float(rate)/ccc*100;
                    fout.write('%s\t%s\t%s\t%s\n'%(refDB, coef, dbfile, rate));
                    print('%s\t%s\t%s\t%s\n'%(refDB, coef, dbfile, rate));
                    
        #break;            
        fout.close();
