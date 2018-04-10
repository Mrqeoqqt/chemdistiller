# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 11:47:54 2017

@author: ilaponog
"""
import os;
import h5py;


dbpath='e:/DB_HDF5/';


infile='e:/Imperial/TestDB/testmassranges_0.005_batch0-4.txt';
    

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
            #print(len(massranges));
            
            massrange=massranges.pop(0);            
            index=0;
            subcount=0;
            fout=open('g:/CFM/TestInChI_All/%s.%s.inchi'%(dbfile,index),'w');
            ascii_dataset=hdf5_container['%s%s'%(subf,'/ChemInfo/ASCII')];
            inchi_dataset=hdf5_container['%s%s'%(subf,'/ChemInfo/InChi')];
            mz_dataset=hdf5_container['%s%s'%(subf,'/MZ')];
            total_count=len(mz_dataset);
            print('Total: %s'%total_count);
            for i in range(total_count):
                
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
                    if subcount>1000:
                        print('%s of %s'%(index, total_count//1000));
                        subcount=0;
                        fout.close();
                        index+=1;
                        fout=open('g:/CFM/TestInChI_All/%s.%s.inchi'%(dbfile,index),'w');
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
                        if subcount>1000:
                            print('%s of %s'%(index, total_count//1000));
                            subcount=0;
                            fout.close();
                            index+=1;
                            fout=open('g:/CFM/TestInChI_All/%s.%s.inchi'%(dbfile,index),'w');
                    elif mz>massrange[1]:
                        break;
                        
            fout.close();            
        
print('Finished')        
        