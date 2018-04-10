# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 11:47:54 2017

@author: ilaponog
"""
import os;
import h5py;
import sys;


if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."));
    sys.path.append(chemdistiller_path);



from chemdistiller.utils.inchi import parse_inchi, \
                                      inchikeyvalues_from_inchi, \
                                      inchikey_from_inchi,\
                                      inchikey_from_inchikeyvalues;


dbpath='e:/DB_HDF5/';


infile='e:/Imperial/TestDB/testmassranges_0.005_batch0-4.txt';
    

onlyfiles = [ff for ff in os.listdir(dbpath) if os.path.isfile(os.path.join(dbpath, ff))];
#%%
for i in reversed(range(len(onlyfiles))):
    fn=onlyfiles[i];
    print(fn)
    if not (fn.startswith('TestDB') or fn.startswith('MassBank') or fn.startswith('ChEBI') or fn.startswith('HMDB')):
        del onlyfiles[i];
print(onlyfiles)

#%%
fout=open('g:/MetFrag/smallDB0/db0.psv','w');                
fout.write('MonoisotopicMass|InChI|Identifier|InChIKey2|InChIKey1|MolecularFormula\n');

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
            #index=0;
            #subcount=0;
            
            ascii_dataset=hdf5_container['%s%s'%(subf,'/ChemInfo/ASCII')];
            inchi_dataset=hdf5_container['%s%s'%(subf,'/ChemInfo/InChi')];
            ids_dataset=hdf5_container['%s%s'%(subf,'/ChemInfo/IDs')];
            inchikey_values_dataset=hdf5_container['%s%s'%(subf,'/ChemInfo/InChiKeyValues')];   
            mz_dataset=hdf5_container['%s%s'%(subf,'/MZ')];
            total_count=len(mz_dataset);
            print('Total: %s'%total_count);
            for i in range(total_count):
                
                mz=mz_dataset[i];
                if mz>=massrange[0] and mz<=massrange[1]:


                    
                    fout.write('%.5f'%mz);
                    
                    j=inchi_dataset[i,:];
                    
                    i0=str(bytearray(ascii_dataset[j[0,0]:j[0,1]]));
                    i1=str(bytearray(ascii_dataset[j[1,0]:j[1,1]]));
                    i2=str(bytearray(ascii_dataset[j[2,0]:j[2,1]]));
                    i3=str(bytearray(ascii_dataset[j[3,0]:j[3,1]]));
                
                    inchi='%s/%s%s%s'%(i0,i1,i2,i3);
                    fout.write('|InChI=1S/%s'%(inchi));
                    
                    j=ids_dataset[i,:];
                    ids=str(bytearray(ascii_dataset[j[0]:j[1]]));
        
                    kv=inchikey_values_dataset[i,:];
                    key=inchikey_from_inchikeyvalues(kv);
                    key=key.split('-');
                    
                    fout.write('|%s|%s|%s|%s\n'%(ids, key[1],key[0],i0));
                    
                    
                    
                    
                    

                elif mz>massrange[1]:
                    while (mz>massrange[1]) and (len(massranges)>0):
                        print('MassRanges Left %s'%len(massranges));
                        massrange=massranges.pop(0);
                    if mz>=massrange[0] and mz<=massrange[1]:

                        fout.write('%s'%mz);
                        
                        j=inchi_dataset[i,:];
                        
                        i0=str(bytearray(ascii_dataset[j[0,0]:j[0,1]]));
                        i1=str(bytearray(ascii_dataset[j[1,0]:j[1,1]]));
                        i2=str(bytearray(ascii_dataset[j[2,0]:j[2,1]]));
                        i3=str(bytearray(ascii_dataset[j[3,0]:j[3,1]]));
                    
                        inchi='%s/%s%s%s'%(i0,i1,i2,i3);
                        fout.write('|InChI=1S/%s'%(inchi));
                        
                        j=ids_dataset[i,:];
                        ids=str(bytearray(ascii_dataset[j[0]:j[1]]));
            
                        kv=inchikey_values_dataset[i,:];
                        key=inchikey_from_inchikeyvalues(kv);
                        key=key.split('-');
                        
                        fout.write('|%s|%s|%s|%s\n'%(ids, key[1],key[0],i0));

                    elif mz>massrange[1]:
                        break;
                        
            
fout.close();                    
print('Finished')        
        