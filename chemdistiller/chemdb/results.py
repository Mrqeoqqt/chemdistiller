# -*- coding: utf-8 -*-
"""
Created on Wed Oct 05 15:01:29 2016

@author: ilaponog
"""

import os;
import sys;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);

import bisect;
from operator import itemgetter;
from chemdistiller.utils.utils import numpy_array_to_string, string_to_numpy_byte_array,\
                                       string_to_numpy_uint16_array, string_to_float_list;
from chemdistiller.utils.base64 import encode_to_base64, decode_from_base64;
from chemdistiller.utils.periodictable import parse_formula;
from chemdistiller.msspectra.adducts import get_adduct_by_name;
          

class MolecularRecord(dict):
    def _pipe_from_textfile(self, finp):
        while True:
            s=finp.readline();
            if s=='':
                 return;
            s=s.rstrip('\n').lstrip();
            if '##' in s:
                s=s[:s.index('##')];
            
            if '=' in s:
                s=s.split('=',1);
                if s[0].lower().startswith('totalscore'):
                    self['TotalScore']=float(s[1]);

                elif s[0].lower().startswith('adduct'):
                    self['Adduct']=s[1];
                elif s[0].lower().startswith('isotopeextramass'):
                    self['IsotopeExtraMass']=float(s[1]);
                elif s[0].lower().startswith('isotope'):
                    self['Isotope']=int(s[1]);


                elif s[0].lower().startswith('mz'):
                    self['MZ']=float(s[1]);
                elif s[0].lower().startswith('mass'):
                    self['Mass']=float(s[1]);
                elif s[0].lower().startswith('charge'):
                    self['Charge']=int(s[1]);
                elif s[0].lower().startswith('dbformat'):
                    self['DBFormat']=int(s[1]);
                elif s[0].lower().startswith('dbindex'):
                    self['DBIndex']=int(s[1]);
                elif s[0].lower().startswith('rindex'):
                    self['RIndex']=int(s[1]);
                elif s[0].lower().startswith('dbname'):
                    self['DBName']=s[1];
                elif s[0].lower().startswith('rfile'):
                    self['RFile']=s[1];
                elif s[0].lower().startswith('smiles'):
                    self['SMILES']=s[1];
                elif s[0].lower().startswith('ids'):
                    self['IDs']=s[1];
                elif s[0].lower().startswith('annotation'):
                    self['Annotation']=s[1];
                elif s[0].lower().startswith('shortinchi'):
                    self['ShortInChI']=s[1];
                elif s[0].lower().startswith('inchikeyvalues'):
                    self['InChIKeyValues']=string_to_numpy_byte_array(s[1]);
                elif s[0].lower().startswith('inchikey'):
                    self['InChiKey']=s[1];
                elif s[0].lower().startswith('inchi'):
                    self['InChI']=s[1];
                elif s[0].lower().startswith('formulavector'):
                    self['FormulaVector']=string_to_numpy_uint16_array(s[1]);
                elif s[0].lower().startswith('elementvector'):
                    self['ElementVector']=string_to_numpy_byte_array(s[1]);
                elif s[0].lower().startswith('frag'):
                    self['Frag']=string_to_float_list(s[1]);
                elif s[0].lower().startswith('formula'):
                    self['Formula']=parse_formula(s[1]);
                elif s[0].lower().startswith('fpt'):
                    self['FPT']=decode_from_base64(s[1]);
                elif s[0].lower().startswith('scores'):
                    if not ('Scores' in self):
                        self['Scores']={};
                    score=s[1].split(':',1);
                    self['Scores'][score[0]]=float(score[1]);
                    
                    
                    
                    
                    
                    

            elif s.lower().startswith('end'):
                return;
    

    def _pipe_to_textfile(self, fout, indent=''):                
        for key in sorted(self.keys()):
            if key in ('TotalScore', 'MZ', 'Mass', 'Charge', 'DBFormat', 'DBIndex', 'RIndex',\
            'DBName', 'RFile', 'SMILES', 'IDs', 'Annotation', 'ShortInChI', 'InChI', 'InChiKey',\
            'Adduct', 'Isotope', 'IsotopeExtraMass'):
                fout.write('%s%s=%s\n'%(indent,key.lower(), self[key]));
            elif key in ('FormulaVector', 'ElementVector', 'Frag', 'InChIKeyValues'):
                fout.write('%s%s=%s\n'%(indent,key.lower(), numpy_array_to_string(self[key])));
            elif key=='Formula':
                fout.write('%s%s=%s\n'%(indent,key.lower(), self[key].formula_to_string()));
            elif key=='FPT':
                fout.write('%s%s=%s\n'%(indent,key.lower(), encode_to_base64(self[key])));
            elif key=='Scores':
                for score in self[key].keys():
                    fout.write('%s%s=%s:%s\n'%(indent,key.lower(), score, self[key][score]));

                
    def __repr__(self):
        s='MolecularRecord(\n';
        for key in sorted(self.keys()):
            s='%s%s = %s\n'%(s,key,self[key]);
        return '%s)'%s;
        

class DBQueryResult:
    
    def __init__(self, results_limit=0, revisitable=True):
        self.results_limit=results_limit;
        self.revisitable=revisitable;
        self.mol_list=[];
        self._mol_scores=[];
        self._minscore=0;
        self.addition_type=0;
        self.total_candidate_count=0;

    def _pipe_from_textfile(self, finp):
        while True:
            s=finp.readline();
            if s=='':
                 return;
            s=s.rstrip('\n').lstrip();
            if '#' in s:
                s=s[:s.index('#')];
            
            if '=' in s:
                s=s.split('=',1);
                if s[0].lower().startswith('results_limit'):
                    self.results_limit=int(s[1]);
                elif s[0].lower().startswith('revisitable'):
                    self.revisitable=s[1].startswith('True');
                elif s[0].lower().startswith('addition_type'):
                    self.addition_type=int(s[1]);
                elif s[0].lower().startswith('total_candidate_count'):
                    self.total_candidate_count=int(s[1]);
                
            elif s.lower().startswith('candidate'):
                candidate=MolecularRecord();
                candidate._pipe_from_textfile(finp);
                self.mol_list.append(candidate);

            elif s.lower().startswith('end'):
                return;

        
    def _pipe_to_textfile(self, fout, indent=''):                

        fout.write('%s%s=%s\n'%(indent,'results_limit',self.results_limit));
        fout.write('%s%s=%s\n'%(indent,'revisitable',self.revisitable));
        fout.write('%s%s=%s\n'%(indent,'addition_type',self.addition_type));
        fout.write('%s%s=%s\n'%(indent,'total_candidate_count',self.total_candidate_count));
                
        if self.mol_list:    
            for candidate in self.mol_list:
               fout.write('%s%s\n'%(indent,'candidate'));
               candidate._pipe_to_textfile(fout,'%s\t'%indent);
               fout.write('%send\n'%indent);

    def __enter__(self):
        return self;

    def __exit__(self, exc_type, exc_value, traceback):
        return 0;        

    '''def rescore(self, filters, scorers, required_fields, total_score):    
        for record in self.mol_list:
            tot_score_value=total_score(record['Scores']);
            record['TotalScore']=tot_score_value;
    '''        
    
    def add(self, db, mz, delta_ppm, charge, filters, scorers, required_fields, total_score):
        if self.revisitable:
            required_fields.add('Revisitable');
            
        self.addition_type=0;
        if charge is None:
            #Process Neutrals
            db_iterator=db.retrieve_records(mz, delta_ppm, 0, filters, scorers, required_fields, filter_on_charge=False);
            for record in db_iterator:
                tot_score_value=total_score(record['Scores']);
                self.total_candidate_count+=1;
                
                if tot_score_value>=self._minscore:            
                    insertindex=bisect.bisect_left(self._mol_scores,tot_score_value);
                    
                    record['TotalScore']=tot_score_value;
                    
                    self.mol_list.insert(insertindex,record);
                    
                    self._mol_scores.insert(insertindex,tot_score_value);
                    if self.results_limit>0:
                        if len(self.mol_list)>self.results_limit:
                            del self.mol_list[0];
                            del self._mol_scores[0];
                            self._minscore=self._mol_scores[0];
                            
            #Process Negatives
            db_iterator=db.retrieve_records(mz, delta_ppm, -1, filters, scorers, required_fields, filter_on_charge=False);
            for record in db_iterator:
                tot_score_value=total_score(record['Scores']);
                self.total_candidate_count+=1;
                if tot_score_value>=self._minscore:            
                    insertindex=bisect.bisect_left(self._mol_scores,tot_score_value);
                    
                    record['TotalScore']=tot_score_value;
                    
                    self.mol_list.insert(insertindex,record);
                    
                    self._mol_scores.insert(insertindex,tot_score_value);
                    if self.results_limit>0:
                        if len(self.mol_list)>self.results_limit:
                            del self.mol_list[0];
                            del self._mol_scores[0];
                            self._minscore=self._mol_scores[0];
            

            #Process Positives
            db_iterator=db.retrieve_records(mz, delta_ppm, 1, filters, scorers, required_fields, filter_on_charge=False);
            for record in db_iterator:
                tot_score_value=total_score(record['Scores']);
                self.total_candidate_count+=1;
                if tot_score_value>=self._minscore:            
                    insertindex=bisect.bisect_left(self._mol_scores,tot_score_value);
                    
                    record['TotalScore']=tot_score_value;
                    
                    self.mol_list.insert(insertindex,record);
                    
                    self._mol_scores.insert(insertindex,tot_score_value);
                    if self.results_limit>0:
                        if len(self.mol_list)>self.results_limit:
                            del self.mol_list[0];
                            del self._mol_scores[0];
                            self._minscore=self._mol_scores[0];
        
        else:
            #Defined charge only
            db_iterator=db.retrieve_records(mz, delta_ppm, charge, filters, scorers, required_fields);
            for record in db_iterator:
                tot_score_value=total_score(record['Scores']);
                self.total_candidate_count+=1;
                if tot_score_value>=self._minscore:            
                    insertindex=bisect.bisect_left(self._mol_scores,tot_score_value);
                    
                    record['TotalScore']=tot_score_value;
                    
                    self.mol_list.insert(insertindex,record);
                    
                    self._mol_scores.insert(insertindex,tot_score_value);
                    if self.results_limit>0:
                        if len(self.mol_list)>self.results_limit:
                            del self.mol_list[0];
                            del self._mol_scores[0];
                            self._minscore=self._mol_scores[0];
            
    def addbulk(self, db, mz, delta_ppm, charge, filters, scorers, required_fields, total_score):
        if self.revisitable:
            required_fields.add('Revisitable');

        self.addition_type=1;
        if charge is None:
            #Process Neutrals
            db_iterator=db.retrieve_records(mz, delta_ppm, 0, filters, scorers, required_fields, filter_on_charge=False);
            for record in db_iterator:
                tot_score_value=total_score(record['Scores']);
                record['TotalScore']=tot_score_value;
                self.total_candidate_count+=1;
                
                self.mol_list.append(record);
            #Process Negatives
            db_iterator=db.retrieve_records(mz, delta_ppm, -1, filters, scorers, required_fields, filter_on_charge=False);
            for record in db_iterator:
                tot_score_value=total_score(record['Scores']);
                record['TotalScore']=tot_score_value;
                self.total_candidate_count+=1;
                
                self.mol_list.append(record);
            #Process Positives
            db_iterator=db.retrieve_records(mz, delta_ppm, 1, filters, scorers, required_fields, filter_on_charge=False);
            for record in db_iterator:
                tot_score_value=total_score(record['Scores']);
                record['TotalScore']=tot_score_value;
                self.total_candidate_count+=1;
                
                self.mol_list.append(record);
        else:    
            #Defined charge only
            db_iterator=db.retrieve_records(mz, delta_ppm, charge, filters, scorers, required_fields);
            for record in db_iterator:
                tot_score_value=total_score(record['Scores']);
                record['TotalScore']=tot_score_value;
                self.total_candidate_count+=1;
                
                self.mol_list.append(record);
                
                
    def cleanup_results(self):
        if self.addition_type==0:
            self.mol_list=list(reversed(self.mol_list));
        else:
            if self.results_limit>0:
                self.mol_list=sorted(self.mol_list, key=itemgetter('TotalScore'),reverse=True)[:self.results_limit];
            else:
                self.mol_list=sorted(self.mol_list, key=itemgetter('TotalScore'),reverse=True);                
            