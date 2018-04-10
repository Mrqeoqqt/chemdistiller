# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 22:08:15 2016

@author: Dr. Ivan Laponogov
"""
import os;
import numpy as np;
import scipy.sparse as sp;


def get_parameters(inpath, index, modelfile, svm_type, settings, fpt_stats):
    #settings[minallowed,maxallowed]: occupancy, prediction_score, w0, w1, c, g, nsvm0, nsvm1, totalsvm, ntrainingset
    finp=open(modelfile,'r');
    count=0;
    seton=0;
    for s in finp:
        s=int(s.rstrip('\n').split(' ')[0]);
        count+=1;
        if s==1:
            seton+=1;
    finp.close();        
    fpt_stats[-1][0]=float(seton)/count;
    fpt_stats[-1][9]=count;
    if 'linear' in svm_type:
        finp=open(os.path.join(inpath,'%s.model.linear.do.bat'%index),'r');
    elif 'radial' in svm_type:
        finp=open(os.path.join(inpath,'%s.model.radial.do.bat'%index),'r');        
    on=False;        
    for s in finp:
        s=s.rstrip('\n');
        if ('linear weighted' in s) and (svmtype=='linear_weight'):
            on=True;
        elif  ('linear notweighted' in s) and (svmtype=='linear_noweight'):
            on=True;
        elif  ('radial not weighted' in s) and (svmtype=='radial_noweight'):
            on=True;
        elif  ('radial weighted' in s) and (svmtype=='radial_weight'):
            on=True;
        
        if on and ('Best' in s):
            s=s.split(' ');
            for ss in s:
                if 'Score' in ss:
                    ss=ss.split('=')[1];
                    fpt_stats[-1][1]=float(ss);
                elif ss.startswith('c'):
                    ss=ss.split('=')[1];
                    fpt_stats[-1][4]=float(ss);
                elif ss.startswith('g'):
                    ss=ss.split('=')[1];
                    fpt_stats[-1][5]=float(ss);

        elif s.startswith('./svm-train') and on:
            s=s.split(' ');
            for i in range(len(s)):
                ss=s[i];
                if ss.startswith('-w0'):
                    
                    fpt_stats[-1][2]=float(s[i+1]);
                elif ss.startswith('-w1'):
                    
                    fpt_stats[-1][3]=float(s[i+1]);
            
            break;

    finp.close();
    
    for i in range(10):
        if settings[i][0]>fpt_stats[-1][i]:
            settings[i][0]=fpt_stats[-1][i];
    
        if settings[i][1]<fpt_stats[-1][i]:
            settings[i][1]=fpt_stats[-1][i];
    
    
    
    
    
    
    

def convert_model(infname, outfname, settings, fpt_stats):
    finp=open(infname,'r');
    fout=open(outfname,'wb');
    reading_sv=False;
    count1=np.int32(0);
    count2=np.int32(0);
    v1=[];
    sparse1=[];
    cc=0;
    gamma=np.float32(0.0);
    for s in finp:
        s=s.rstrip('\n');
        if reading_sv:
            s=s.split(' ');
            cc+=1;
            if 1==1:
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
                    kernel=1;
                elif s=='linear':
                    kernel=0;
                else:
                    raise TypeError('Wrong SVM kernel!');
            elif s.startswith('gamma'):
                gamma=np.float32(s.split(' ')[1]);
            elif s.startswith('nr_class'):
                s=s.split(' ')[1];
                if s!='2':
                    raise TypeError('Wrong SVM number of classes!');
            elif s.startswith('rho'):
                rho=np.float64(s.split(' ')[1]);
            elif s.startswith('label'):
                s=s.split(' ');
                label1=np.uint8(s[1]);
                label2=np.uint8(s[2]);
            elif s.startswith('nr_sv'):
                s=s.split(' ');
                count1=np.int32(s[1]);
                count2=np.int32(s[2]);
            elif s.startswith('SV'):
                reading_sv=True;

    finp.close();
    if len(v1)!=count1+count2:
        raise TypeError('Wrong SVM format!');
    stats=[0.0]*10;    
    stats[5]=gamma;
    stats[6]=count1;
    stats[7]=count2;
    stats[8]=count1+count2;
    fpt_stats.append(stats);
    
    np.uint8(kernel).tofile(fout);
    v1=np.array(v1,dtype=np.float32);

    if label1!=1:
        label2,label1=label1,label2;
        count1,count2=count2,count1;
        v1=v1*-1;
        rho=rho*-1;
        
    rho.tofile(fout);
    count1.tofile(fout);
    count2.tofile(fout);
        
    label1.tofile(fout);
    label2.tofile(fout);
    
    if kernel==0:

        sv_coef=sp.csc_matrix(v1);

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
        sparse1=sp.csc_matrix((data, (row, col)), shape=(count1+count2, 20001));
        w=sv_coef.dot(sparse1);        
        np.int32(w.data.shape[0]).tofile(fout);
        w.data.tofile(fout);
        np.int32(w.indices.shape[0]).tofile(fout);
        w.indices.tofile(fout);
        np.int32(w.indptr.shape[0]).tofile(fout);
        w.indptr.tofile(fout);
        
    else:
    
        gamma.tofile(fout);
        v1.tofile(fout);
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
        sparse1=sp.csc_matrix((data, (row, col)), shape=(count1+count2, 20001));
    
        np.int32(sparse1.data.shape[0]).tofile(fout);
        sparse1.data.tofile(fout);
        np.int32(sparse1.indices.shape[0]).tofile(fout);
        sparse1.indices.tofile(fout);
        np.int32(sparse1.indptr.shape[0]).tofile(fout);
        sparse1.indptr.tofile(fout);

    fout.close();


def condense_linear_svm_ws(inpath, outpath, batch, mode, svm_type):
    outpath=os.path.join(outpath, svm_type, str(batch),str(mode));
    inpath=os.path.join(inpath, svm_type, str(batch),str(mode));
    
    finp=open(os.path.join(inpath,'fpt_mask.default'),'rb');
    fpt_mask=np.unpackbits(np.fromfile(finp,dtype=np.uint8));
    finp.close();
    if not os.path.exists(outpath):
        os.makedirs(outpath);
    ws=[];rr=[];
    for i in range(11416):
        if fpt_mask[i]==1:
            fname=os.path.join(inpath,'%s.svm'%i);
            with open(fname,'rb') as finp:
                print(fname);
                
                kernel=np.fromfile(finp, dtype=np.uint8, count=1)[0];

                if kernel!=0:
                    raise TypeError('Should be linear SVM!');
            
            
                rho=np.fromfile(finp, dtype=np.float64, count=1)[0];
                rr.append(rho);

                count1=np.fromfile(finp, dtype=np.int32, count=1)[0];
                count2=np.fromfile(finp, dtype=np.int32, count=1)[0];
                label1=np.fromfile(finp, dtype=np.uint8, count=1)[0];
                label2=np.fromfile(finp, dtype=np.uint8, count=1)[0];


                datalen=np.fromfile(finp, dtype=np.int32, count=1)[0];
                data=np.fromfile(finp, dtype=np.float32, count=datalen);
                indiceslen=np.fromfile(finp, dtype=np.int32, count=1)[0];
                indices=np.fromfile(finp, dtype=np.int32, count=indiceslen);
                
                indptrlen=np.fromfile(finp, dtype=np.int32, count=1)[0];
                indptr=np.fromfile(finp, dtype=np.int32, count=indptrlen);
        
                w=sp.csc_matrix((data, indices, indptr), shape=(1, 20001));
                ws.append(w);
    ww=sp.vstack(ws);
    ww=sp.csc_matrix(ww);
    
    with open(os.path.join(outpath,'w.dat'),'wb') as fout:
        kernel.tofile(fout);
        count1=np.int32(len(rr));
        rrho=np.array(rr, dtype=np.float64);
        count1.tofile(fout);
        rrho.tofile(fout);
        
        np.int32(ww.data.shape[0]).tofile(fout);
        ww.data.tofile(fout);
        np.int32(ww.indices.shape[0]).tofile(fout);
        ww.indices.tofile(fout);
        np.int32(ww.indptr.shape[0]).tofile(fout);
        ww.indptr.tofile(fout);
    
    
                    
        
        

def process_svm_type(inpath, outpath, batch, mode, svm_type):
    if mode==-1:
        modename='Negative';
    else:
        modename='Positive';
    inpath=os.path.join(inpath,str(batch), modename);
    outpath=os.path.join(outpath, svm_type, str(batch),str(mode));
    if not os.path.exists(outpath):
        os.makedirs(outpath);
    fpt_mask=np.zeros((11416,), dtype=np.uint8);
    settings=[];
    for i in range(10):
        settings.append([1e30,-1e30]);
    fpt_stats=[];        
    for i in range(11416):
        fname=os.path.join(inpath,'%s.model.%s'%(i,svm_type));
        if os.path.isfile(fname):
            fpt_mask[i]=1;
            outfname=os.path.join(outpath,'%s.svm'%i);
            print(outfname);
            convert_model(fname,outfname,settings,fpt_stats);
            
            get_parameters(inpath,i,os.path.join(inpath,'%s.model'%i), svm_type, settings, fpt_stats);
            
    fpt_mask=np.packbits(fpt_mask);
    fout=open(os.path.join(outpath,'fpt_mask.default'),'wb');
    fpt_mask.tofile(fout);
    fout.close();
    
    fout=open(os.path.join(outpath,'fpt_stats.dat'),'wb');
    if 'linear' in svm_type:
        np.array([0],dtype=np.uint8).tofile(fout);
    else:
        np.array([1],dtype=np.uint8).tofile(fout);
        
    np.array(fpt_stats,dtype=np.float32).tofile(fout);
    fout.close();

    fout=open(os.path.join(outpath,'stats.txt'),'w');
    fout.write('[%s,%s], #occupancy\n'%(settings[0][0],settings[0][1]));
    fout.write('[%s,%s], #prediction_score\n'%(settings[1][0],settings[1][1]));
    fout.write('[%s,%s], #w0\n'%(settings[2][0],settings[2][1]));
    fout.write('[%s,%s], #w1\n'%(settings[3][0],settings[3][1]));
    fout.write('[%s,%s], #c\n'%(settings[4][0],settings[4][1]));
    fout.write('[%s,%s], #g\n'%(settings[5][0],settings[5][1]));
    fout.write('[%s,%s], #nsvm0\n'%(settings[6][0],settings[6][1]));
    fout.write('[%s,%s], #nsvm1\n'%(settings[7][0],settings[7][1]));
    fout.write('[%s,%s], #totalsvm\n'%(settings[8][0],settings[8][1]));
    fout.write('[%s,%s], #ntrainingset\n'%(settings[9][0],settings[9][1]));
    fout.close();
    
    
    

inpath='e:/Imperial/TestDB/SVM_Train';

outpath='e:/Imperial/TestDB/SVM_Train/Ported_csc';
   
batch=-1;

mode=-1;
svmtype='linear_noweight';
process_svm_type(inpath, outpath, batch, mode, svmtype);
svmtype='linear_weight';
process_svm_type(inpath, outpath, batch, mode, svmtype);

mode=1;
svmtype='linear_noweight';
process_svm_type(inpath, outpath, batch, mode, svmtype);
svmtype='linear_weight';
process_svm_type(inpath, outpath, batch, mode, svmtype);

mode=-1;
svmtype='radial_noweight';
process_svm_type(inpath, outpath, batch, mode, svmtype);
svmtype='radial_weight';
process_svm_type(inpath, outpath, batch, mode, svmtype);

mode=1;
svmtype='radial_noweight';
process_svm_type(inpath, outpath, batch, mode, svmtype);
svmtype='radial_weight';
process_svm_type(inpath, outpath, batch, mode, svmtype);

mode=1;

svmtype='linear_noweight';
condense_linear_svm_ws(outpath, outpath, batch, mode, svmtype);
svmtype='linear_weight';
condense_linear_svm_ws(outpath, outpath, batch, mode, svmtype);


mode=-1;
svmtype='linear_noweight';
condense_linear_svm_ws(outpath, outpath, batch, mode, svmtype);
svmtype='linear_weight';
condense_linear_svm_ws(outpath, outpath, batch, mode, svmtype);


