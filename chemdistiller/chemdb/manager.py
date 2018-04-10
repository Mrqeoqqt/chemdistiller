# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 22:26:27 2016

Chemical Database Manager.

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

from chemdistiller.chemdb.database import DBase;
from chemdistiller.chemdb.results import DBQueryResult;

from chemdistiller.scorers.totalscore import total_multiplicative_score;
from chemdistiller.filters.element import CHNOPS_filter, CHNOPS_halogens_filter, ElementCompositionFilter;
from chemdistiller.filters.formula import FormulasFilter;
from chemdistiller.scorers.formula import FormulaScorer;
from chemdistiller.scorers.element import ElementScorer;


class DBManager:
    
    def __init__(self, database_list_file=''):
        if database_list_file=='':
            subpath=os.getcwd();
            #print(subpath);
            if os.path.isfile(os.path.join(subpath,'databases.list')):
                database_list_file=os.path.join(subpath,'databases.list');
            elif os.path.isfile(os.path.abspath(os.path.join(subpath,'data/databases.list'))):
                database_list_file=os.path.abspath(os.path.join(subpath,'data/databases.list'));
            else:
                raise IOError('databases.list not found in %s or %s!'%(subpath,os.path.abspath(os.path.join(subpath,'Data'))));
                
        self.data_bases=[];
        self.initialize_databases(database_list_file);
        
    def __enter__(self):
        return self;

    def close(self):
        if self.data_bases:
            for db in self.data_bases:
                db.close();
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.close();
        
    def initialize_databases(self, database_list_file):
        with open(database_list_file,'r') as f:
            for s in f:
                subpath=s.rstrip('\n');
                if os.path.isabs(subpath):
                    dbpath=subpath;
                else:
                    dbpath=os.path.abspath(os.path.join(os.path.dirname(database_list_file),subpath));
                onlyfiles = [ff for ff in os.listdir(dbpath) if os.path.isfile(os.path.join(dbpath, ff))];
                
                for dbfile in onlyfiles:
                    if dbfile.lower().endswith('.h5') or dbfile.lower().endswith('.dbinfo.dat'):
                        dbase=DBase(os.path.join(dbpath, dbfile));
                        if dbase.db_format>0 and dbase.db_format<4:
                            self.data_bases.append(dbase);
                            self.data_bases[len(self.data_bases)-1].db_index=len(self.data_bases)-1;
                        else:
                            print('%s - Wrong DB format. Ignored.'%dbpath);
                
    def get_database_names(self, dbindexes=[]):
        result=[];
        if dbindexes:
            for index in dbindexes:
                result.append(self.data_bases[index].db_name);
        else:
            for db in self.data_bases:
                result.append(db.db_name);
        return result;
    
    def get_database_fullnames(self, dbindexes=[]):
        result=[];
        if dbindexes:
            for index in dbindexes:
                result.append(self.data_bases[index].db_fullname);
        else:
            for db in self.data_bases:
                result.append(db.db_fullname);
        return result;
        
    def get_database_descriptions(self, dbindexes=[]):
        result=[];
        if dbindexes:
            for index in dbindexes:
                result.append(self.data_bases[index].db_description);
        else:
            for db in self.data_bases:
                result.append(db.db_description);
        return result;

        
    def get_database_licenses(self, dbindexes=[]):
        result=[];
        if dbindexes:
            for index in dbindexes:
                result.append(self.data_bases[index].db_license);
        else:
            for db in self.data_bases:
                result.append(db.db_license);
        return result;
    
    def get_database_source_urls(self, dbindexes=[]):
        result=[];
        if dbindexes:
            for index in dbindexes:
                result.append(self.data_bases[index].db_source_url);
        else:
            for db in self.data_bases:
                result.append(db.db_source_url);
        return result;
        
    def get_database_citations(self, dbindexes=[]):        
        result=[];
        if dbindexes:
            for index in dbindexes:
                result.append(self.data_bases[index].db_citation);
        else:
            for db in self.data_bases:
                result.append(db.db_citation);
        return result;
                
    def get_database_versions(self, dbindexes=[]):        
        result=[];
        if dbindexes:
            for index in dbindexes:
                result.append(self.data_bases[index].db_version);
        else:
            for db in self.data_bases:
                result.append(db.db_version);
        return result;
        
    def get_database_formats(self, dbindexes=[]):        
        result=[];
        if dbindexes:
            for index in dbindexes:
                result.append(self.data_bases[index].db_format);
        else:
            for db in self.data_bases:
                result.append(db.db_format);
        return result;
        

    def get_database_subgroups(self, dbindexes=[]):        
        result=[];
        if dbindexes:
            for index in dbindexes:
                result.append(self.data_bases[index].db_subgroupindex);
        else:
            for db in self.data_bases:
                result.append(db.db_subgroupindex);
        return result;
    
    def get_database_groupnames(self, dbindexes=[]):        
        result=[];
        if dbindexes:
            for index in dbindexes:
                result.append(self.data_bases[index].db_groupname);
        else:
            for db in self.data_bases:
                result.append(db.db_groupname);
        return result;    
        
    def db_indexes_from_db_names(self, db_names_list, case_sensitive=False):
        result=[];
        if case_sensitive:
            for db_name in db_names_list:
                found=False;
                for db in self.data_bases:
                    if db.db_name==db_name:
                        result.append(db.db_index);
                        found=True;
                        break;
                if not found:
                    raise IndexError('Database name %s not found in the list of databases!'%db_name);
        else:
            for db_name in db_names_list:
                found=False;
                for db in self.data_bases:
                    if db.db_name.upper()==db_name.upper():
                        result.append(db.db_index);
                        found=True;
                        break;
                if not found:
                    raise IndexError('Database name %s not found in the list of databases!'%db_name);
        return result;
                

    def query_by_mz_scored(self, mz, deltappm, charge=None, db_indexes=[], filters=[], scorers=[], total_score=total_multiplicative_score, required_fields=None, results_limit=0, save_memory=True, revisitable=True):

        if required_fields is None:
            required_fields=set();
        else:
            required_fields=required_fields.copy();
        for s in scorers:
            required_fields.update(s.required_fields);
        for s in filters:
            required_fields.update(s.required_fields);

        result=DBQueryResult(results_limit, revisitable);
        

        if db_indexes:
            #print('Indexes')
            for index in db_indexes:
                if save_memory:
                    result.add(self.data_bases[index], mz, deltappm, charge, filters, scorers, required_fields, total_score);
                else:
                    result.addbulk(self.data_bases[index], mz, deltappm, charge, filters, scorers, required_fields, total_score);                    
        else:
            for db in self.data_bases:
                #print(db.db_name)
                if save_memory:                
                    result.add(db, mz, deltappm, charge, filters, scorers, required_fields, total_score);
                else:
                    result.addbulk(db, mz, deltappm, charge, filters, scorers, required_fields, total_score);

        result.cleanup_results();

        return result;
    
if __name__=='__main__':   
     
     dbfile=os.path.join(chemdistiller_path,'tests/databases_DEMO.list');

     
     print(dbfile);
     if os.path.isfile(dbfile):
         dbmanager=DBManager(dbfile);
         import time;
         t=time.time();
         
         print(dbmanager.db_indexes_from_db_names(['HMDB']));
         
         result=dbmanager.query_by_mz_scored(174.125594, 20, charge=0, \
             db_indexes=dbmanager.db_indexes_from_db_names(['HMDB'],case_sensitive=True),\
             filters=[], scorers=[],\
             required_fields=set(['Formula','SMILES','InChI','Frag']), results_limit=20, save_memory=False, revisitable=True);
         t=time.time()-t;
         print(len(result.mol_list));
         for i in result.mol_list:
             print(i);
         print(t);
         
         dbmanager.close();
     

