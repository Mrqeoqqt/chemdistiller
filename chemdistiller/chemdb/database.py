# -*- coding: utf-8 -*-
"""
Created on Wed Oct 05 14:16:14 2016
DBase class


@author: Dr. Ivan Laponogov
"""
import os;
import sys;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);




from chemdistiller.utils.periodictable import parse_formula, \
                                              formula_to_element_vector, \
                                              encode_formula_to_array, \
                                              decode_formula_from_array;

from chemdistiller.utils.base64 import decode_from_base64;

from chemdistiller.utils.inchi import parse_inchi, \
                                      inchikeyvalues_from_inchi, \
                                      inchikey_from_inchi,\
                                      inchikey_from_inchikeyvalues;

                                      
from chemdistiller.chemdb.results import MolecularRecord;                                      
from chemdistiller.utils.sysutils import python_version;                                      

import h5py;
import numpy as np;


def parse_string_fragments(s):
    s=s.split(',');
    return [float(i) for i in s];
    
def parse_string_fragment_charges(charge, frag, fragcharges):
    frag=frag.split(',');        
    lst=[];
    inbracket=0;
    ss='';
    if '[' in fragcharges: 
        for i in fragcharges:
            if i=='[':
                inbracket=1;
                ss='';
            elif i==']':
                inbracket=0;
                lst.append([int(i) for i in ss.split(',')]);
                ss='';
            elif inbracket==1:
                ss+=i;
    else:
        lst.append([int(fragcharges)]);
    
    if len(lst)!=len(frag):
        raise TypeError('Different lengths in Frag (%s) and FragCharges (%s)!'%(len(frag),len(lst)));
    
    
    fragnew=set();
    for i in range(len(frag)):
        fv=float(frag[i]);
        for j in lst[i]:
            if j==0:
                fragnew.add(fv);
            elif (j>0 and charge>0) or (j<0 and charge<0):
                fragnew.add(abs(fv/j));
            
    return sorted(list(fragnew),reverse=True);


class DBase:

    def __init__(self,database):
        self.database_path, self.database_name=os.path.split(database);
        self.db_name='';
        self.db_fullname='';        
        self.db_description='';
        self.db_license='';
        self.db_source_url='';
        self.db_citation='';
        self.db_version='';
        self.db_format=0;
        self.db_subgroupindex=0;
        self.db_groupname='';
        self.db_index=0;
        self.record_index=-1;
        self.currentfile='';
        self._initialize_database();
        
    def close(self):
        return 0;
        
    def _initialize_database(self):
        
        if self.database_name.lower().endswith('.dbinfo.dat'):
            fname=os.path.join(self.database_path, self.database_name);
            if os.path.isfile(fname):
                with open(fname,'r') as f:
                    for s in f:
                        s=s.rstrip('\n');
                        if   s.upper().startswith('DBNAME'):
                            self.db_name=s.split('=',1)[1];
                        elif s.upper().startswith('DBFULLNAME'):
                            self.db_fullname=s.split('=',1)[1];
                        elif s.upper().startswith('DBDESCRIPTION'):
                            self.db_description=s.split('=',1)[1];
                        elif s.upper().startswith('DBLICENSE'):
                            self.db_license=s.split('=',1)[1];
                        elif s.upper().startswith('DBSOURCEURL'):
                            self.db_source_url=s.split('=',1)[1];
                        elif s.upper().startswith('DBCITATION'):
                            self.db_citation=s.split('=',1)[1];
                        elif s.upper().startswith('DBVERSION'):
                            self.db_version=s.split('=',1)[1];
                        elif s.upper().startswith('DBFORMAT'):
                            self.db_format=int(s.split('=',1)[1]);
                        elif s.upper().startswith('DBSUBGROUPINDEX'):
                            self.db_subgroupindex=int(s.split('=',1)[1]);
                        elif s.upper().startswith('DBGROUPNAME'):
                            self.db_groupname=s.split('=',1)[1];
        elif self.database_name.lower().endswith('.h5'):   
            fname=os.path.join(self.database_path, self.database_name);
            if os.path.isfile(fname):
                self.hdf5_container=h5py.File(fname,'r');
                if 'HDF5ContainerType' in self.hdf5_container.attrs:
                    #print(self.hdf5_container.attrs['HDF5ContainerType'])
                    if self.hdf5_container.attrs['HDF5ContainerType'].decode('ascii')=='DistilledChemicalDatabase':
                        self.db_format=3;
                if self.db_format!=3:
                    self.hdf5_container.close();
        
        if self.db_format==1:
            self._init_raw();
        elif self.db_format==3:
            self._init_hdf5();
        else:
            print('Format unsupported!');

        
    def __enter__(self):
        return self;

    def __exit__(self, exc_type, exc_value, traceback):
        self.close();
        
    def _close_raw(self):
        if hasattr(self,'datafile'):
            if not self.datafile.closed:
                self.datafile.close();
    
    def _close_hdf5(self):
        if hasattr(self,'hdf5_container'):
            self.hdf5_container.close();

    def _init_raw(self):
        self.retrieve_records=self._raw_retrieve_records;
        self._next_record=self._raw_next_record;
        self.close=self._close_raw;
        #self.fpt_mask=np.full((11416//8,1),255,np.uint8);
        
    def _init_hdf5(self):        
        self.retrieve_records=self._hdf5_retrieve_records;
        self._next_record=self._hdf5_next_record;
        self.close=self._close_hdf5;
        self.hdf5_version=str(self.hdf5_container.attrs['HDF5ContainerVersion']);
        
        if 'DBName' in self.hdf5_container.attrs:
            self.db_name=self.hdf5_container.attrs['DBName'];
            if python_version==3:
                self.db_name=self.db_name.decode('ascii', "backslashreplace");
                                                
        if 'DBFullName' in self.hdf5_container.attrs:
            self.db_fullname=self.hdf5_container.attrs['DBFullName'];
            if python_version==3:
                self.db_fullname=self.db_fullname.decode('ascii', "backslashreplace");
                                                      
                                                      
        if 'DBDESCRIPTION' in self.hdf5_container.attrs:
            self.db_description=self.hdf5_container.attrs['DBDESCRIPTION'];
            if python_version==3:
                self.db_description=self.db_description.decode('ascii', "backslashreplace");
                                                         
                                                         
        if 'DBLICENSE' in self.hdf5_container.attrs:
            self.db_license=self.hdf5_container.attrs['DBLICENSE'];
            if python_version==3:
                self.db_license=self.db_license.decode('ascii', "backslashreplace");
                                                     
                                                     
        if 'DBSOURCEURL' in self.hdf5_container.attrs:
            self.db_source_url=self.hdf5_container.attrs['DBSOURCEURL'];
            if python_version==3:
                self.db_source_url=self.db_source_url.decode('ascii', "backslashreplace");
                                                        
                                                        
        if 'DBCITATION' in self.hdf5_container.attrs:
            self.db_citation=self.hdf5_container.attrs['DBCITATION'];
            if python_version==3:
                self.db_citation=self.db_citation.decode('ascii', "backslashreplace");
                                                      
        if 'DBVERSION' in self.hdf5_container.attrs:
            self.db_version=self.hdf5_container.attrs['DBVERSION'];
            if python_version==3:
                self.db_version=self.db_version.decode('ascii', "backslashreplace");
                                                     
        if 'DBSUBGROUPINDEX' in self.hdf5_container.attrs:
            self.db_subgroupindex=int(self.hdf5_container.attrs['DBSUBGROUPINDEX']);
                                     
        if 'DBGROUPNAME' in self.hdf5_container.attrs:
            self.db_groupname=self.hdf5_container.attrs['DBGROUPNAME'];
            if python_version==3:
                self.db_groupname=self.db_groupname.decode('ascii', "backslashreplace");
                                                       

    def _raw_record_accepted(self,record):
        mz=record['MZ'];
        #print(mz);
        if mz<self.mzmin:
            return False;
        if mz>self.mzmax:
            raise StopIteration();
        if self.filter_on_charge: 
            if self.charged:
                if abs(self.charge-record['Charge'])>1e-8:
                    return False;
        for f in self.filters:
            if f.rejected(record):
                return False;
        #record['DBPath']=self.database_path;
        if 'Revisitable' in self.required_fields:
            record['DBName']=self.db_name;
            record['DBFormat']=self.db_format;
            record['RFile']=self.currentfile;
            record['RIndex']=self.record_index;
            record['DBIndex']=self.db_index;
        return True;
        
        
    def _hdf5_record_accepted(self,record):
        if self.filter_on_charge: 
            if self.charged:
                if abs(self.charge-record['Charge'])>1e-8:
                    return False;
        for f in self.filters:
            if f.rejected(record):
                return False;
        #record['DBPath']=self.database_path;
        if 'Revisitable' in self.required_fields:
            record['DBName']=self.db_name;
            record['DBFormat']=self.db_format;
            record['RFile']=os.path.join(self.database_path, self.database_name);
            record['RIndex']=self.mzcurrent;
            record['DBIndex']=self.db_index;
        return True;
        
    
    def _get_next_hdf5_record(self):
        self.mzcurrent+=1;
        if self.mzcurrent>self.mzmax_index:
            raise StopIteration();                

        record=MolecularRecord();
        
        if self.charged:
            record['MZ']=self.mz_mass_charge_dataset[self.mzcurrent,0]; 
            record['Mass']=self.mz_mass_charge_dataset[self.mzcurrent,1]; 
            record['Charge']=int(self.mz_mass_charge_dataset[self.mzcurrent,2]); 
        else:
            record['MZ']=self.mz_dataset[self.mzcurrent]; 
            record['Mass']=record['MZ'];
            record['Charge']=0;

        if 'FPT' in self.required_fields:
            record['FPT']=self.fpt_dataset[self.mzcurrent,:];

        if 'Frag' in self.required_fields:
            i0=self.fragprint_index_dataset[self.mzcurrent,0];
            i1=self.fragprint_index_dataset[self.mzcurrent,1];
            
            record['Frag']=list(self.fragprint_values_dataset[i0:i1]); 
            
        
        if 'CFM' in self.required_fields:
            
            i0=self.CFM_negative_index_dataset[self.mzcurrent,0,0];
            j0=self.CFM_negative_index_dataset[self.mzcurrent,0,1];

            i1=self.CFM_negative_index_dataset[self.mzcurrent,1,0];
            j1=self.CFM_negative_index_dataset[self.mzcurrent,1,1];

            i2=self.CFM_negative_index_dataset[self.mzcurrent,2,0];
            j2=self.CFM_negative_index_dataset[self.mzcurrent,2,1];
            
            if i0>=0 and i1>=0 and i2>=0:
                n=[self.CFM_values_dataset[i0:j0,:], self.CFM_values_dataset[i1:j1,:], self.CFM_values_dataset[i2:j2,:]]
            else:
                n=[];

            i0=self.CFM_positive_index_dataset[self.mzcurrent,0,0];
            j0=self.CFM_positive_index_dataset[self.mzcurrent,0,1];

            i1=self.CFM_positive_index_dataset[self.mzcurrent,1,0];
            j1=self.CFM_positive_index_dataset[self.mzcurrent,1,1];

            i2=self.CFM_positive_index_dataset[self.mzcurrent,2,0];
            j2=self.CFM_positive_index_dataset[self.mzcurrent,2,1];

            if i0>=0 and i1>=0 and i2>=0:
                p=[self.CFM_values_dataset[i0:j0,:], self.CFM_values_dataset[i1:j1,:], self.CFM_values_dataset[i2:j2,:]]
            else:
                p=[];
            
            record['CFM']=[n, p];


        if 'CFM2' in self.required_fields:
            
            i0=self.CFM2_negative_index_dataset[self.mzcurrent,0,0];
            j0=self.CFM2_negative_index_dataset[self.mzcurrent,0,1];

            i1=self.CFM2_negative_index_dataset[self.mzcurrent,1,0];
            j1=self.CFM2_negative_index_dataset[self.mzcurrent,1,1];

            i2=self.CFM2_negative_index_dataset[self.mzcurrent,2,0];
            j2=self.CFM2_negative_index_dataset[self.mzcurrent,2,1];
            
            if i0>=0 and i1>=0 and i2>=0:
                n=[self.CFM2_values_dataset[i0:j0,:], self.CFM2_values_dataset[i1:j1,:], self.CFM2_values_dataset[i2:j2,:]]
            else:
                n=[];

            i0=self.CFM2_positive_index_dataset[self.mzcurrent,0,0];
            j0=self.CFM2_positive_index_dataset[self.mzcurrent,0,1];

            i1=self.CFM2_positive_index_dataset[self.mzcurrent,1,0];
            j1=self.CFM2_positive_index_dataset[self.mzcurrent,1,1];

            i2=self.CFM2_positive_index_dataset[self.mzcurrent,2,0];
            j2=self.CFM2_positive_index_dataset[self.mzcurrent,2,1];

            if i0>=0 and i1>=0 and i2>=0:
                p=[self.CFM2_values_dataset[i0:j0,:], self.CFM2_values_dataset[i1:j1,:], self.CFM2_values_dataset[i2:j2,:]]
            else:
                p=[];
            
            record['CFM2']=[n, p];

            
        if 'SMILES' in self.required_fields:    
            i0=self.smiles_dataset[self.mzcurrent,0];
            i1=self.smiles_dataset[self.mzcurrent,1];
            arr=self.ascii_dataset[i0:i1];
            if python_version==2:                                  
                record['SMILES']=str(bytearray(arr));
            else:
                record['SMILES']=bytes(list(arr)).decode('ascii');

        if 'IDs' in self.required_fields:            
            i0=self.ids_dataset[self.mzcurrent,0];
            i1=self.ids_dataset[self.mzcurrent,1];
            arr=self.ascii_dataset[i0:i1];
            if python_version==2:
                record['IDs']=str(bytearray(arr));
            else:
                record['IDs']=bytes(list(arr)).decode('ascii');


        if ('InChIKeyValues' in self.required_fields) or ('InChIKey' in self.required_fields):        
            kv=self.inchikey_values_dataset[self.mzcurrent,:];
            if 'InChIKeyValues' in self.required_fields:        
                record['InChIKeyValues']=kv;
                      
                      
            
            if 'InChIKey' in self.required_fields:        
                record['InChIKey']=inchikey_from_inchikeyvalues(kv);

        if ('ShortInChI' in self.required_fields) or ('InChI' in self.required_fields):# or ('Formula' in self.required_fields):            
            i=self.inchi_dataset[self.mzcurrent,0,0];
            j=self.inchi_dataset[self.mzcurrent,0,1];
            arr=self.ascii_dataset[i:j];
            if python_version==2:
                i0=str(bytearray(arr));
            else:
                i0=bytes(list(arr)).decode('ascii');

            i=self.inchi_dataset[self.mzcurrent,1,0];
            j=self.inchi_dataset[self.mzcurrent,1,1];
            arr=self.ascii_dataset[i:j];
            if python_version==2:
                i1=str(bytearray(arr));
            else:
                i1=bytes(list(arr)).decode('ascii');

            i=self.inchi_dataset[self.mzcurrent,2,0];
            j=self.inchi_dataset[self.mzcurrent,2,1];
            arr=self.ascii_dataset[i:j];
            if python_version==2:
                i2=str(bytearray(arr));
            else:
                i2=bytes(list(arr)).decode('ascii');

            i=self.inchi_dataset[self.mzcurrent,3,0];
            j=self.inchi_dataset[self.mzcurrent,3,1];
            arr=self.ascii_dataset[i:j];
            if python_version==2:
                i3=str(bytearray(arr));
            else:
                i3=bytes(list(arr)).decode('ascii');
            
            #if 'Formula' in self.required_fields:    #check which is faster
            #    record['Formula']=parse_formula(i0);
            if 'ShortInChI' in self.required_fields:
                record['ShortInChI']='%s/%s'%(i0,i1);
            if 'InChI' in self.required_fields:
                record['InChI']='%s/%s%s%s'%(i0,i1,i2,i3);
        
        if ('FormulaVector' in self.required_fields) or ('Formula' in self.required_fields):       
            fv=self.formulavector_dataset[self.mzcurrent,:];
            if ('FormulaVector' in self.required_fields):       
                record['FormulaVector']=fv;
            if ('Formula' in self.required_fields):       
                record['Formula']=decode_formula_from_array(fv);
                
        if 'ElementVector' in self.required_fields:
            record['ElementVector']=self.elementvector_dataset[self.mzcurrent,:];
            
        #Likeness
        if ('MetaLike' in self.required_fields):
            record['MetaLike']=self.metalikeness_dataset[self.mzcurrent];
            
        #End of likeness    
            
        return record;
            

    
    def _get_next_raw_record(self):
        if self.currentfile=='':
            self.currentfile=os.path.join(self.database_path, self.db_name, self.subf,
                                          str(self.mzcurrent//1000),
                                          str(self.mzcurrent%1000//100),
                                          str(self.mzcurrent%100//10),
                                          '%s.st2'%str(self.mzcurrent%10));
            while (not os.path.isfile(self.currentfile)) and (self.mzcurrent<=self.mzmax_int):
                self.mzcurrent+=1;
                self.currentfile=os.path.join(self.database_path, self.db_name, self.subf,
                                          str(self.mzcurrent//1000),
                                          str(self.mzcurrent%1000//100),
                                          str(self.mzcurrent%100//10),
                                          '%s.st2'%str(self.mzcurrent%10));
            if self.mzcurrent<=self.mzmax_int:
                self.datafile=open(self.currentfile,'r');
                self.record_index=-1;
            else:
                raise StopIteration();                
        
        s=self.datafile.readline();
        self.record_index+=1;
        while s=='':
            self.datafile.close();
            self.mzcurrent+=1;
            self.currentfile=os.path.join(self.database_path, self.db_name, self.subf,
                                          str(self.mzcurrent//1000),
                                          str(self.mzcurrent%1000//100),
                                          str(self.mzcurrent%100//10),
                                          '%s.st2'%str(self.mzcurrent%10));
            while (not os.path.isfile(self.currentfile)) and (self.mzcurrent<=self.mzmax_int):
                self.mzcurrent+=1;
                self.currentfile=os.path.join(self.database_path, self.db_name, self.subf,
                                          str(self.mzcurrent//1000),
                                          str(self.mzcurrent%1000//100),
                                          str(self.mzcurrent%100//10),
                                          '%s.st2'%str(self.mzcurrent%10));
            if self.mzcurrent<=self.mzmax_int:
                self.datafile=open(self.currentfile,'r');
                self.record_index=-1;
            else:
                raise StopIteration();                
            s=self.datafile.readline();
            self.record_index+=1;
            
        s=s.rstrip('\n').split('\t');
        record=MolecularRecord();
        record['MZ']=float(s[0]);
        if self.charged:
            record['Mass']=float(s[1]);
            record['Charge']=float(s[2]);
        else:
            record['Mass']=record['MZ'];
            record['Charge']=0;
        if 'ShortInChI' in self.required_fields:            
            record['ShortInChI']=parse_inchi(s[2+self.offs])[0];
            
        if 'InChI' in self.required_fields:
            record['InChI']=s[2+self.offs];

        if 'SMILES' in self.required_fields:    
            record['SMILES']=s[3+self.offs];

        if 'IDs' in self.required_fields:            
            record['IDs']=s[4+self.offs];

        if 'FPT' in self.required_fields:
            record['FPT']=decode_from_base64(s[5+self.offs]);
            # Mask FPT here !

        if 'Frag' in self.required_fields:
            record['Frag']=s[6+self.offs];
            if self.charged:
                record['FragCharge']=s[9];

        if 'InChIKeyValues' in self.required_fields:        
            record['InChIKeyValues']=inchikeyvalues_from_inchi(s[2+self.offs]);
            
        if 'InChIKey' in self.required_fields:        
            record['InChIKey']=inchikey_from_inchi(s[2+self.offs]);
        
        if ('Formula' in self.required_fields) or ('ElementVector' in self.required_fields) or ('FormulaVector' in self.required_fields):       
            fla=parse_formula(s[1+self.offs].split('/')[0]);
            
            if 'Formula' in self.required_fields:
                record['Formula']=fla;
        
            if 'ElementVector' in self.required_fields:
                record['ElementVector']=formula_to_element_vector(fla);
            
            if 'FormulaVector' in self.required_fields:
                record['FormulaVector']=encode_formula_to_array(fla);
        

            
        return record;
        
        
    def _raw_next_record(self):
        record=self._get_next_raw_record();

        while not self._raw_record_accepted(record):
            record=self._get_next_raw_record();
        if 'Frag' in self.required_fields:
            if self.charged:
                record['Frag']=parse_string_fragment_charges(record['Charge'], record['Frag'], record['FragCharge']);
                del record['FragCharge'];
            else:
                record['Frag']=parse_string_fragments(record['Frag']);
                
        record['Scores']={};
        for scorer in self.scorers:
            scorer.process_molecular_candidate_record(self,record);
        return record;
        

    def _hdf5_next_record(self):
        record=self._get_next_hdf5_record();

        while not self._hdf5_record_accepted(record):
            record=self._get_next_hdf5_record();
                
        record['Scores']={};
        for scorer in self.scorers:
            scorer.process_molecular_candidate_record(self,record);
        return record;
        


    def _raw_retrieve_records(self, mz, delta_ppm, charge, filters, scorers, required_fields, filter_on_charge=True):
        self.charge=charge;
        self.filter_on_charge=filter_on_charge;
        if charge<0:
            self.subf='Negative';
            self.charged=True;
            self.offs=2;
        elif charge>0:
            self.subf='Positive';
            self.charged=True;
            self.offs=2;
        else:
            self.subf='Neutral';
            self.charged=False;
            self.offs=0;
        self.mzmin=mz*(1-1/1000000.0*delta_ppm);
        self.mzmax=mz*(1+1/1000000.0*delta_ppm);
        self.mzmin_int=int(self.mzmin*100);
        self.mzmax_int=int(self.mzmax*100);
        self.mzcurrent=self.mzmin_int;
        self.filters=filters;
        self.scorers=scorers;
        self.required_fields=required_fields;
        
        if self.currentfile!='':
            if hasattr(self,'datafile'):
                if not self.datafile.closed:
                    self.datafile.close();

            self.currentfile='';
        return DBRecordIterator(self);


    def __hdf5_get_top_index(self):
        if self.charged:
            lo = 0;
            hi = self.mz_mass_charge_dataset.shape[0];
            while lo < hi:
                mid = (lo+hi)//2;
                if self.mz_mass_charge_dataset[mid,0] < self.mzmin: 
                    lo = mid+1;
                else: 
                    hi = mid;
            return lo;
        else:
            lo = 0;
            hi = self.mz_dataset.shape[0];
            while lo < hi:
                mid = (lo+hi)//2;
                if self.mz_dataset[mid] < self.mzmin: 
                    lo = mid+1;
                else: 
                    hi = mid;
            return lo;
        
        
        
    def __hdf5_get_bottom_index(self):
        if self.charged:
            lo = 0;
            hi = self.mz_mass_charge_dataset.shape[0];
            while lo < hi:
                mid = (lo+hi)//2;
                if self.mzmax < self.mz_mass_charge_dataset[mid,0]: 
                    hi = mid;
                else: 
                    lo = mid+1;
            return lo-1;
        else:
            lo = 0;
            hi = self.mz_dataset.shape[0];
            while lo < hi:
                mid = (lo+hi)//2;
                if self.mzmax < self.mz_dataset[mid]: 
                    hi = mid;
                else: 
                    lo = mid+1;
            return lo-1;
    
        
        
        
        
    def _hdf5_retrieve_records(self, mz, delta_ppm, charge, filters, scorers, required_fields, filter_on_charge=True):        
        self.charge=charge;
        self.filter_on_charge=filter_on_charge;
        self.required_fields=required_fields;
        self.filters=filters;
        self.scorers=scorers;
        
        if charge<0:
            self.subf='/Negative';
            self.charged=True;
        elif charge>0:
            self.subf='/Positive';
            self.charged=True;
        else:
            self.subf='/Neutral';
            self.charged=False;

        self.mzmin=mz*(1-1/1000000.0*delta_ppm);
        self.mzmax=mz*(1+1/1000000.0*delta_ppm);
        use_ascii=False;

        if self.charged:
            self.mz_mass_charge_dataset=self.hdf5_container['%s%s'%(self.subf,'/MZMassCharge')];
        else:
            self.mz_dataset=self.hdf5_container['%s%s'%(self.subf,'/MZ')];

        if 'FPT' in self.required_fields:
            self.fpt_dataset=self.hdf5_container['%s%s'%(self.subf,'/FingerPrints/FPTArray')];
            
        if 'Frag' in self.required_fields:
            self.fragprint_index_dataset=self.hdf5_container['%s%s'%(self.subf,'/FragPrints/FragPrintIndex')];
            self.fragprint_values_dataset=self.hdf5_container['%s%s'%(self.subf,'/FragPrints/FragPrintValues')];
            
        if 'CFM' in self.required_fields:
            self.CFM_negative_index_dataset=self.hdf5_container['%s%s'%(self.subf, '/CFM/SE/Negative')];
            self.CFM_positive_index_dataset=self.hdf5_container['%s%s'%(self.subf, '/CFM/SE/Positive')];
            self.CFM_values_dataset=self.hdf5_container['%s%s'%(self.subf, '/CFM/datavector')];


        if 'CFM2' in self.required_fields:
            self.CFM2_negative_index_dataset=self.hdf5_container['%s%s'%(self.subf, '/CFM/SE/Negative')];
            self.CFM2_positive_index_dataset=self.hdf5_container['%s%s'%(self.subf, '/CFM/SE/Positive')];
            self.CFM2_values_dataset=self.hdf5_container['%s%s'%(self.subf, '/CFM/datavector')];

            
        if ('ShortInChI' in self.required_fields) or ('InChI' in self.required_fields):# or ('Formula' in self.required_fields):            
            self.inchi_dataset=self.hdf5_container['%s%s'%(self.subf,'/ChemInfo/InChi')];
            use_ascii=True;
        
        if 'SMILES' in self.required_fields:    
            self.smiles_dataset=self.hdf5_container['%s%s'%(self.subf,'/ChemInfo/SMILES')];
            use_ascii=True;
        
        if 'IDs' in self.required_fields:            
            self.ids_dataset=self.hdf5_container['%s%s'%(self.subf,'/ChemInfo/IDs')];
            use_ascii=True;
        
        if ('InChIKeyValues' in self.required_fields) or ('InChIKey' in self.required_fields):        
            self.inchikey_values_dataset=self.hdf5_container['%s%s'%(self.subf,'/ChemInfo/InChiKeyValues')];
            
        if use_ascii:
            self.ascii_dataset=self.hdf5_container['%s%s'%(self.subf,'/ChemInfo/ASCII')];
        
        if 'ElementVector' in self.required_fields:
            self.elementvector_dataset=self.hdf5_container['%s%s'%(self.subf,'/ChemInfo/ElementsVector')];
                            
        if ('FormulaVector' in self.required_fields) or ('Formula' in self.required_fields):
            self.formulavector_dataset=self.hdf5_container['%s%s'%(self.subf,'/ChemInfo/FormulaVector')];
        
        
        #Likeness section
        
        if ('MetaLike' in self.required_fields):
            self.metalikeness_dataset=self.hdf5_container['%s%s'%(self.subf,'/Likeness/MetaLikeness')];
        
        
        #Likeness section end
        
        
        self.mzmin_index=self.__hdf5_get_top_index();
        self.mzmax_index=self.__hdf5_get_bottom_index();
        self.mzcurrent=self.mzmin_index-1;
        
        return DBRecordIterator(self);
        


class DBRecordIterator:
    def __init__(self,database):
        self.database=database;

    def __iter__(self):
        return self;

    def next(self):
        return self.database._next_record();
    
    def __next__(self):
        return self.database._next_record();
    


        