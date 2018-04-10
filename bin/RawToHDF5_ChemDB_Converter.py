# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 15:42:17 2016

@author: Dr. Ivan Laponogov
"""
import sys;
import os;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."));
    sys.path.append(chemdistiller_path);
    
import h5py;
import numpy as np;
import time;
#from glob import glob;

chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."));
sys.path.append(chemdistiller_path);

from chemdistiller.utils.base64 import decode_from_base64;
from chemdistiller.utils.inchi import inchikeyvalues_from_inchi, parse_inchi;
from chemdistiller.utils.periodictable import parse_formula, formula_to_element_vector, encode_formula_to_array;
                        
from chemdistiller.chemdb.database import parse_string_fragment_charges, parse_string_fragments;

    
    
    

class HDF5DBContainer:

    def __init__(self, output_filename):
        
        self.folderpath, self.filename = os.path.split(os.path.abspath(output_filename));
        if not os.path.exists(self.folderpath):
            os.makedirs(self.folderpath);
        self.HDF5container=h5py.File(os.path.join(self.folderpath,self.filename));
        self.fragprintpos=0;
        self.asciipos=0;

    def __enter__(self):
        return self;

    def __exit__(self, exc_type, exc_value, traceback):
        if hasattr(self,'HDF5container'):
            self.HDF5container.close();
            
            
    def hdf5_import_from_st2raw(self, inpath, fptmask=np.ones((11416,), np.uint8)):
        if not os.path.isfile(os.path.join(inpath,'dbinfo.dat')):
            raise IOError('Database info file not found: %s'%os.path.join(inpath,'dbinfo.dat'));
        
        self.HDF5container.attrs['HDF5ContainerType']=np.string_('DistilledChemicalDatabase');
        self.HDF5container.attrs['HDF5ContainerVersion']=np.string_('1.0');
        finp=open(os.path.join(inpath,'dbinfo.dat'),'r');
        #fout=open(os.path.join(self.folderpath,'dbinfo.dat'),'w');
        for s in finp:
            s=s.rstrip('\n').lstrip().split('=',1);
            if s[0].upper()=='DBFORMAT':
                s[1]='3';
            if s[0]!='':
                #fout.write('%s=%s\n'%(s[0],s[1]));
                self.HDF5container.attrs[s[0]]=np.string_(s[1]);
                
        #fout.close();
        finp.close();
        
        fptlist=[];
        
        for i in range(11416):
            if fptmask[i]==1:
                fptlist.append(i);
        

        fptlen=len(fptlist);
        fptsubmask=np.packbits(np.ones((fptlen,), np.uint8));
        fptmasklen=len(fptsubmask);
        fptindexes=np.array(fptlist,dtype=np.uint32);
        packedmask=np.packbits(fptmask);
        packedmasklen=len(packedmask);
        
        #hdf5_ascii_string = h5py.special_dtype(vlen=bytes);

        fptgroup=self.HDF5container.create_group('FingerPrints');
        #Original mask, packed
        fptoriginalmask=fptgroup.create_dataset("FPTOriginalMask",(packedmasklen,), maxshape=(packedmasklen,), dtype=np.uint8);
        fptoriginalmask[:]=packedmask[:];
        
        #List of indeces of original FPT bits (11416)
        fptmask=fptgroup.create_dataset("FPTMask",(fptlen,), maxshape=(fptlen,), compression="gzip", compression_opts=4,dtype=np.uint32);
        fptmask[:]=fptindexes[:];
        
        #Mask for working bits (packed)       
        fptsubmask=fptgroup.create_dataset("FPTsubmask",(1, fptmasklen), chunks=(100, fptmasklen), maxshape=(None, fptmasklen), compression="gzip", compression_opts=4,dtype=np.uint8);        
        
        #FPT info: 0 - original bit count=11416, 1 - length of new fpt after masking, 2 - length of the packed fpt, 3 - No of padding bits
        fptinfo=fptgroup.create_dataset("FPTInfo",(4, ), maxshape=(4,), dtype=np.uint32);        
        fptinfo[0]=11416;
        fptinfo[1]=fptlen;
        fptinfo[2]=fptmasklen;
        fptinfo[3]=fptmasklen*8-fptlen;
                

        print('Listing input files');
        
        subpaths=['/Negative','/Positive','/Neutral'];
        #subpaths=['/Positive'];
        for subpath in subpaths:
            print(subpath);
            
            fptgroup=self.HDF5container.create_group(subpath+'/FingerPrints');
            fraggroup=self.HDF5container.create_group(subpath+'/FragPrints');
            chemgroup=self.HDF5container.create_group(subpath+'/ChemInfo');
            chargegroup=self.HDF5container[subpath];
            
            #New FPT Array, packed and trimmed to fptmask
            fptdataset=fptgroup.create_dataset("FPTArray",(1, fptmasklen), chunks=(100, fptmasklen), maxshape=(None, fptmasklen), compression="gzip", compression_opts=4,dtype=np.uint8);
            
            if subpath!='/Neutral':
                masschargedataset=chargegroup.create_dataset("MZMassCharge",(1, 3),chunks=(10000, 3), maxshape=(None, 3), compression="gzip", compression_opts=4, dtype=np.float32);
            else:
                mzdataset=chargegroup.create_dataset("MZ",(1,),chunks=(10000,), maxshape=(None, ), compression="gzip", compression_opts=4, dtype=np.float32);
            
            inchikey_dataset=chemgroup.create_dataset("InChiKeyValues",(1, 15),chunks=(10000, 15), maxshape=(None, 15), compression="gzip", compression_opts=4, dtype=np.uint8);
            
            elementsvector_dataset=chemgroup.create_dataset("ElementsVector",(1, 12),chunks=(10000, 12), maxshape=(None, 12), compression="gzip", compression_opts=4, dtype=np.uint8);
            
            formulavector_dataset=chemgroup.create_dataset("FormulaVector",(1, 96),chunks=(10000, 96), maxshape=(None, 96), compression="gzip", compression_opts=4, dtype=np.uint16);
            
            fragprintindex_dataset=fraggroup.create_dataset("FragPrintIndex",(1,2),chunks=(10000,2), maxshape=(None, 2), compression="gzip", compression_opts=4, dtype=np.int64);
            fragprintvalues_dataset=fraggroup.create_dataset("FragPrintValues",(1,),chunks=(10000,), maxshape=(None, ), compression="gzip", compression_opts=4, dtype=np.float32);
            
            smiles_dataset=chemgroup.create_dataset("SMILES",(1, 2),chunks=(10000, 2), maxshape=(None, 2), compression="gzip", compression_opts=4, dtype=np.int64);

            inchi_dataset=chemgroup.create_dataset("InChi",(1, 4, 2),chunks=(10000, 4, 2), maxshape=(None, 4, 2), compression="gzip", compression_opts=4, dtype=np.int64);

            ids_dataset=chemgroup.create_dataset("IDs",(1, 2),chunks=(10000, 2), maxshape=(None, 2), compression="gzip", compression_opts=4, dtype=np.int64);
            
            ascii_dataset=chemgroup.create_dataset("ASCII",(1,),chunks=(10000,), maxshape=(None, ), compression="gzip", compression_opts=4, dtype=np.uint8);

        
            recordindex=-1;
            
            fileslist=[];
            
            
            for i in range(0,2000):
                if os.path.exists(inpath+subpath+'/%s'%i):
                    print(inpath+subpath+'/%s'%i);
                    for j in range(0,10):
                        if os.path.exists(inpath+subpath+'/%s/%s'%(i,j)):
                            for k in range(0,10):
                                if os.path.exists(inpath+subpath+'/%s/%s/%s'%(i,j,k)):
                                    for l in range(0,10):
                                        if os.path.isfile(inpath+subpath+'/%s/%s/%s/%s.st2'%(i,j,k,l)):
                                            fileslist.append(inpath+subpath+'/%s/%s/%s/%s.st2'%(i,j,k,l));

            print('Total number of input files: %s'%len(fileslist));
        
        
            for filename in fileslist:
                    fpath,fname=os.path.split(filename);
                    subpath=fpath.replace(inpath,'');
                    if 'Neutral' in subpath:
                        charged=False;
                        offs=0;
                    else:
                        charged=True;
                        offs=2;
                    print('Importing: .../%s/%s'%(subpath,fname));
                    dblist=[];
                    with open(filename,'r') as finp:
                        for s in finp:
                            try:
                                s=s.replace('\n','').replace('\r','').split('\t');
                                mz=float(s[0]);
                                if charged:
                                    mass=float(s[1]);
                                    charge=float(s[2]);
                                else:
                                    mass=mz;
                                    charge=0.0;
                                if mass>=12.0:    
                                    #shortinchi=s[1+offs];
                                    inchi=s[2+offs];
                                    smiles=s[3+offs];
                                    ids=s[4+offs];
                                    fpt=s[5+offs];
                                    fpt=decode_from_base64(fpt);
                                    fpt=np.unpackbits(fpt);
                                    frag=s[6+offs];
                                    if charged:
                                        fragcharge=s[9];
                                    else:
                                        fragcharge='';
                                    recordindex+=1;
                                    if recordindex%1000==0:
                                        print('Total: %s'%(recordindex+1));
                                    
                                    dblist.append([recordindex, mz, charged, mass, charge, inchi, fpt, frag, fragcharge, smiles, ids]);
                            except:
                                print('Error! Skipping!');
                                
                    if len(dblist)>0:
                        #expand datasets here
                        fptdataset.resize((recordindex+1, fptmasklen));
    
                        if charged:
                            masschargedataset.resize((recordindex+1, 3));
                        else:
                            mzdataset.resize((recordindex+1,));
                        
                        inchikey_dataset.resize((recordindex+1, 15));
                        
                        elementsvector_dataset.resize((recordindex+1, 12));
                        
                        formulavector_dataset.resize((recordindex+1, 96));
            
                        fragprintindex_dataset.resize((recordindex+1,2));
                        
                        smiles_dataset.resize((recordindex+1, 2));
                        
                        ids_dataset.resize((recordindex+1, 2));
                        
                        inchi_dataset.resize((recordindex+1, 4, 2));
                        
                        
                        for db in dblist:
                            currentindex=db[0];
                            fptdataset[currentindex,:]=np.packbits(db[6][fptindexes])[:];
                            #print(inchi)
                            inchi=parse_inchi(db[5]);
                            #print(inchi)
                            inchikeyvalues=inchikeyvalues_from_inchi(db[5]);
    
                            sformula=inchi[0].split('/',1)[0];
                            
                            #print(sformula);
                            formula=parse_formula(sformula);
                            
                            elementsvector=formula_to_element_vector(formula);
                            encodedformula=encode_formula_to_array(formula);
                            
                            charge=db[4];
                            charged=db[2];

                            if charged:
                                #print(db[7],db[8])
                                frags=parse_string_fragment_charges(charge, db[7], db[8]);
                                #print(frags)
                            else:
                                frags=parse_string_fragments(db[7]);

                            if charged:
                                masschargedataset[currentindex,0]=db[1];
                                masschargedataset[currentindex,1]=db[3];
                                masschargedataset[currentindex,2]=charge;
                            else:
                                mzdataset[currentindex]=db[1];
                            
                            inchikey_dataset[currentindex,:]=inchikeyvalues[:];
                            
                            elementsvector_dataset[currentindex,:]=elementsvector[:];
                            
                            formulavector_dataset[currentindex,:]=encodedformula[:];
                            
                            fragcount=len(frags);
                            frags=np.array(frags,dtype=np.float32);
                            fragprintindex_dataset[currentindex,0]=self.fragprintpos;
                            fragprintindex_dataset[currentindex,1]=self.fragprintpos+fragcount;
                            
                            fragprintvalues_dataset.resize((self.fragprintpos+fragcount,));
                            
                            fragprintvalues_dataset[self.fragprintpos:self.fragprintpos+fragcount]=frags[:];
                            
                            self.fragprintpos+=fragcount;
                            
                            

                            smiles=bytearray(db[9].encode('ascii'));
                            
                            smileslen=len(smiles);
                            
                            smiles=np.array(smiles,dtype=np.uint8);

                            ids=bytearray(db[10].encode('ascii'));
                            
                            idslen=len(ids);
                            
                            ids=np.array(ids,dtype=np.uint8);
                            
                            sinchi=inchi[0].split('/',1);
                            if len(sinchi)>1:
                                sinchi=sinchi[1];
                            else:
                                sinchi='';
                            
                            inchi0=bytearray(sformula.encode('ascii'));
                            inchi1=bytearray(sinchi.encode('ascii'));
                            inchi2=bytearray(inchi[2].encode('ascii'));
                            inchi3=bytearray(inchi[1].encode('ascii'));
                            
                            inchi0len=len(inchi0);
                            inchi1len=len(inchi1);
                            inchi2len=len(inchi2);
                            inchi3len=len(inchi3);

                            inchi0=np.array(inchi0,dtype=np.uint8);
                            inchi1=np.array(inchi1,dtype=np.uint8);
                            inchi2=np.array(inchi2,dtype=np.uint8);
                            inchi3=np.array(inchi3,dtype=np.uint8);

                            ascii_dataset.resize((self.asciipos+smileslen+idslen+inchi0len+inchi1len+inchi2len+inchi3len,));
                            
                            smiles_dataset[currentindex,0]=self.asciipos;
                            smiles_dataset[currentindex,1]=self.asciipos+smileslen;
                            ascii_dataset[self.asciipos:self.asciipos+smileslen]=smiles[:];
                            self.asciipos+=smileslen;

                            ids_dataset[currentindex,0]=self.asciipos;
                            ids_dataset[currentindex,1]=self.asciipos+idslen;
                            ascii_dataset[self.asciipos:self.asciipos+idslen]=ids[:];
                            self.asciipos+=idslen;

                            inchi_dataset[currentindex,0,0]=self.asciipos;
                            inchi_dataset[currentindex,0,1]=self.asciipos+inchi0len;
                            ascii_dataset[self.asciipos:self.asciipos+inchi0len]=inchi0[:];
                            self.asciipos+=inchi0len;
                            
                            inchi_dataset[currentindex,1,0]=self.asciipos;
                            inchi_dataset[currentindex,1,1]=self.asciipos+inchi1len;
                            ascii_dataset[self.asciipos:self.asciipos+inchi1len]=inchi1[:];
                            self.asciipos+=inchi1len;
                            
                            inchi_dataset[currentindex,2,0]=self.asciipos;
                            inchi_dataset[currentindex,2,1]=self.asciipos+inchi2len;
                            ascii_dataset[self.asciipos:self.asciipos+inchi2len]=inchi2[:];
                            self.asciipos+=inchi2len;
                            
                            inchi_dataset[currentindex,3,0]=self.asciipos;
                            inchi_dataset[currentindex,3,1]=self.asciipos+inchi3len;
                            ascii_dataset[self.asciipos:self.asciipos+inchi3len]=inchi3[:];
                            self.asciipos+=inchi3len;
                            

        print('Import Finished!');
             
        
        
        


if __name__=='__main__':
    #subpaths=[];
    
    
    subpaths=['HMDB','ChEBI','MassBank','BMDB','TestDB',\
    'PlantCyc/aracyc',\
    'PlantCyc/barleycyc',\
    'PlantCyc/brachypodiumcyc',\
    'PlantCyc/cassavacyc',\
    'PlantCyc/chinesecabbagecyc',\
    'PlantCyc/chlamycyc',\
    'PlantCyc/corncyc',\
    'PlantCyc/grapecyc',\
    'PlantCyc/mosscyc',\
    'PlantCyc/oryzacyc',\
    'PlantCyc/papayacyc',\
    'PlantCyc/plantcyc',\
    'PlantCyc/poplarcyc',\
    'PlantCyc/potatocyc',\
    'PlantCyc/selaginellacyc',\
    'PlantCyc/setariacyc',\
    'PlantCyc/sorghumbicolorcyc',\
    'PlantCyc/soycyc',\
    'PlantCyc/spirodelacyc',\
    'PlantCyc/switchgrasscyc',\
    'PlantCyc/tomatocyc',\
    'PlantCyc/wheatacyc',\
    'PlantCyc/wheatdcyc',\
    'DrugBank','ECMDB','FooDB','EMolecules','LipidMaps','MINE/EcoCycMINE',\
    'MINE/KEGGMINE','MINE/YMDBMINE','PhenolExplorer/Compounds','PhenolExplorer/Metabolites',\
    'SMPDB','T3DB','UNPD','YMDB','PubChem','Zinc'\
    ];

    for subpath in subpaths:
        with HDF5DBContainer('g:/DB_HDF5_fix/' + subpath + "_" + time.strftime("[%H_%M_%S#%d_%m_%Y]") + ".h5") as HDF5:
            print(HDF5.folderpath);
            print(HDF5.filename);
            #finp=open('e:/Imperial/TestDB/SVM_Train/Ported_csc/linear_weight/0/-1/fpt_mask.default','rb');
            #mask=np.fromfile(finp,dtype=np.uint8);
            #mask=np.unpackbits(mask);
            #finp.close();
            HDF5.hdf5_import_from_st2raw('e:/DB/'+subpath);
            
    
    
    