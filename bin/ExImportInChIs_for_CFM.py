# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 11:47:54 2017

@author: ilaponog
"""
import os;
import h5py;
import numpy as np;


dbpath='g:/DB_HDF5_CFM/';


infile='e:/Imperial/TestDB/testmassranges_0.005_batch0.txt';
outpath='g:/CFM/TestInChI_round11';
importpath='g:/CFM/TestInChI/test/results4';

batchcount=100;

if not os.path.exists(outpath):
    os.makedirs(outpath);


from glob import glob;
import tarfile;
#import time;





def process_input_data(fname, finp, cfm_dataset, cfm_datavalues):
    fname=os.path.basename(fname);
    fname,ext=os.path.splitext(fname);
    #print(fname);
    index=int(fname);
    print(index);
    if cfm_dataset[index,0,0]>=0:
        print('Already exists!');
        return
    energy0=[];
    energy1=[];
    energy2=[];
    #currentenergy=None;
    for s in finp:
        if s.startswith('energy0'):
            currentenergy=energy0;
        elif s.startswith('energy1'):
            currentenergy=energy1;
        elif s.startswith('energy2'):
            currentenergy=energy2;
        else:
            s=s.rstrip('\n').split();
            if len(s)==2:
                currentenergy.append((float(s[0]),float(s[1])));
    tot=len(energy0)+len(energy1)+len(energy2);
    c=len(cfm_datavalues);
    cfm_datavalues.resize((c+tot,2));
    i=0;
    cfm_dataset[index,0,0]=c+i;
    cfm_dataset[index,0,1]=c+i+len(energy0);
    for record in energy0:
        cfm_datavalues[c+i,0]=record[0];
        cfm_datavalues[c+i,1]=record[1];
        i+=1;
        
    cfm_dataset[index,1,0]=c+i;
    cfm_dataset[index,1,1]=c+i+len(energy1);
    for record in energy1:
        cfm_datavalues[c+i,0]=record[0];
        cfm_datavalues[c+i,1]=record[1];
        i+=1;

    cfm_dataset[index,2,0]=c+i;
    cfm_dataset[index,2,1]=c+i+len(energy2);
    for record in energy2:
        cfm_datavalues[c+i,0]=record[0];
        cfm_datavalues[c+i,1]=record[1];
        i+=1;
        
 
def import_from_log(fname, cfm_dataset, cfm_datavalues):
    #print(fname);
    with open(fname,'r') as finp:
        process_input_data(fname, finp, cfm_dataset, cfm_datavalues);
    
def import_from_gz(fname, cfm_dataset, cfm_datavalues):
    print(fname);
    #time.sleep(1)
    tar = tarfile.open(fname, "r:gz")
    for member in tar.getmembers():
        #print(member);
        f = tar.extractfile(member);
        if not (f is None):
            process_input_data(member.name, f, cfm_dataset, cfm_datavalues);
            f.close();


def try_import(dbfile, importpath, cfm_dataset_neg, cfm_dataset_pos, cfm_datavalues):
    dbfile=dbfile.split('_')[0];
    dbpath=os.path.join(importpath, dbfile);
    if os.path.isdir(dbpath):
        result = [y for x in os.walk(dbpath) for y in glob(os.path.join(x[0], '*.*'))];
        for fname in result:
            fn,fext=os.path.splitext(fname);
            #print(fname);
            if fext=='.log':
                if 'negative' in fname:
                    import_from_log(fname, cfm_dataset_neg, cfm_datavalues);
                else:
                    import_from_log(fname, cfm_dataset_pos, cfm_datavalues);
                    
            elif fext=='.gz':
                if 'negative' in fname:
                    import_from_gz(fname, cfm_dataset_neg, cfm_datavalues);
                else:
                    import_from_gz(fname, cfm_dataset_pos, cfm_datavalues);
                

        

onlyfiles = [ff for ff in os.listdir(dbpath) if os.path.isfile(os.path.join(dbpath, ff))];


for dbfile in onlyfiles:
    if dbfile.lower().endswith('.h5'):
        print(dbfile);
        subf='/Neutral';
        with h5py.File(os.path.join(dbpath,dbfile)) as hdf5_container:
            massranges=[];
            finp=open(infile,'r');
            for s in finp:
                s=s.rstrip('\n').split('\t');
                massranges.append([float(s[0]),float(s[1])]);
            finp.close();
            
            massrange=massranges.pop(0);            
            index=0;
            subcount=0;
            
            ascii_dataset=hdf5_container['%s%s'%(subf,'/ChemInfo/ASCII')];
            inchi_dataset=hdf5_container['%s%s'%(subf,'/ChemInfo/InChi')];
            mz_dataset=hdf5_container['%s%s'%(subf,'/MZ')];
            total_count=len(mz_dataset);
            print('Total: %s'%total_count);
            
            if '%s%s'%(subf,'/CFM/SE/Negative') in hdf5_container:
                cfm_dataset_neg=hdf5_container['%s%s'%(subf,'/CFM/SE/Negative')];
            else:
                cfm_dataset_neg=hdf5_container.create_dataset('%s%s'%(subf,'/CFM/SE/Negative'), (total_count, 3, 2), chunks=(10000, 3, 2), maxshape=(None, 3, 2), compression="gzip", compression_opts=4, dtype=np.int64);
                cfm_dataset_neg[:,:,:]=-1;

            if '%s%s'%(subf,'/CFM/SE/Positive') in hdf5_container:
                cfm_dataset_pos=hdf5_container['%s%s'%(subf,'/CFM/SE/Positive')];
            else:
                cfm_dataset_pos=hdf5_container.create_dataset('%s%s'%(subf,'/CFM/SE/Positive'), (total_count, 3, 2), chunks=(10000, 3, 2), maxshape=(None, 3, 2), compression="gzip", compression_opts=4, dtype=np.int64);
                cfm_dataset_pos[:,:,:]=-1;
            
            if '%s%s'%(subf,'/CFM/datavector') in hdf5_container:
                cfm_datavalues=hdf5_container['%s%s'%(subf,'/CFM/datavector')];
            else:
                cfm_datavalues=hdf5_container.create_dataset('%s%s'%(subf,'/CFM/datavector'), (1, 2), chunks=(50000, 2), maxshape=(None, 2), compression="gzip", compression_opts=4, dtype=np.float32);
            
            try_import(dbfile, importpath, cfm_dataset_neg, cfm_dataset_pos, cfm_datavalues);            
 
            fout=open('%s/%s.%s.inchi'%(outpath, dbfile, index),'w');
            for i in range(total_count):
                if not (cfm_dataset_neg[i,0,0]>=0 and cfm_dataset_pos[i,0,0]>=0):
                    mz=mz_dataset[i];
                    if mz>=massrange[0] and mz<=massrange[1]:
                        j=inchi_dataset[i,:];
                        
                        i0=str(bytearray(ascii_dataset[j[0,0]:j[0,1]]));
                        i1=str(bytearray(ascii_dataset[j[1,0]:j[1,1]]));
                        i2=str(bytearray(ascii_dataset[j[2,0]:j[2,1]]));
                        i3=str(bytearray(ascii_dataset[j[3,0]:j[3,1]]));
                    
                        inchi='%s/%s%s%s'%(i0,i1,i2,i3);
                        fout.write('%s InChI=1S/%s\n'%(i,inchi));
                        subcount+=1;
                        if subcount>batchcount:
                            print('%s of %s'%(index, total_count//batchcount));
                            subcount=0;
                            fout.close();
                            index+=1;
                            fout=open('%s/%s.%s.inchi'%(outpath,dbfile,index),'w');
                    elif mz>massrange[1]:
                        while (mz>massrange[1]) and (len(massranges)>0):
                            print('MassRanges Left %s'%len(massranges));
                            massrange=massranges.pop(0);
                        if mz>=massrange[0] and mz<=massrange[1]:
                            j=inchi_dataset[i,:];
    
                            i0=str(bytearray(ascii_dataset[j[0,0]:j[0,1]]));
                            i1=str(bytearray(ascii_dataset[j[1,0]:j[1,1]]));
                            i2=str(bytearray(ascii_dataset[j[2,0]:j[2,1]]));
                            i3=str(bytearray(ascii_dataset[j[3,0]:j[3,1]]));
    
                            inchi=='%s/%s%s%s'%(i0,i1,i2,i3);
                            fout.write('%s InChI=1S/%s\n'%(i,inchi));
                            subcount+=1;
                            if subcount>batchcount:
                                print('%s of %s'%(index, total_count//batchcount));
                                subcount=0;
                                fout.close();
                                index+=1;
                                fout=open('%s/%s.%s.inchi'%(outpath,dbfile,index),'w');
                        elif mz>massrange[1]:
                            break;
            fout.close();                                        
         
        
print('Finished')        
        