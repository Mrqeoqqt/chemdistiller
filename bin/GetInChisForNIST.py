# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 00:49:48 2016

@author: ilaponog
"""


import sys;
import os;
import pybel;
if sys.byteorder!='little':
    print('Only little endian machines currently supported! bye bye ....');
    quit();

sys.path.append("E:/Imperial/Metaspace_WP4_Source");


#%%

from spectra.spectrum import MSSpectrum;
from spectra.mspeak import MSPeak;
from spectralmanager.spectralmanager import SpectralManager;
from chemdbmanager.chemdatabasemanager import DBManager;
from filters.formulasfilter import FormulasFilter;
#%%

dbfile=os.path.abspath('e:/Imperial/Metaspace_WP4_Source/databases.list');

inchikeys=[];
formulas=[];
exactmasses=[];
finp=open('e:/imperial/temp.txt','r');
for s in finp:
    s=s.rstrip('\n').split('\t');
    inchikeys.append(s[0]);
    formulas.append(s[1]);
    exactmasses.append(float(s[2]));
finp.close();

inchi=['']*len(inchikeys);
shortinchi=['']*len(inchikeys);
smileslist=['']*len(inchikeys);
#%%
charges=['']*len(inchikeys);
     
#%%     
if os.path.isfile(dbfile):
         dbmanager=DBManager(dbfile);
         for index in range(len(inchikeys)):
             key=inchikeys[index];
             #if key=='OZLPUNFFCJDMJD-UHFFFAOYSA-N':
             if inchi[index]=='':
                 #test_formulas_filter=FormulasFilter(formulas[index]);
                 mass=exactmasses[index]/1;
                 print(mass,formulas[index]);
                 
                 result=dbmanager.query_by_mz_scored(mass, 20, charge=1, db_indexes=[4,38], filters=[], scorers=[], required_fields=set(['Formula']), results_limit=-1, save_memory=False);
                 
                 #print(result);
                 for i in result.mol_list:
                     
                     smiles=i['SMILES'];
                     print(smiles);
                     mol=pybel.readstring('smi',smiles);
                     mol.addh();
                     inchikey=mol.write('inchikey').rstrip('\n');
                     print(inchikey);
                     keysp=key.split('-');
                     inchikeysp=inchikey.split('-');
                     if inchikeysp[0]==keysp[0]:
                         inchi[index]=i['InChi'];
                         shortinchi[index]=i['ShortInChi'];
                         smileslist[index]=smiles;
                         charges[index]=mol.charge;
                         formulas[index]=mol.formula;
                         print(smiles,inchi[index],shortinchi[index]);
                         break;
             #elif charges[index]=='':
             #    charges[index]=0;
                 
                 
         
         
         
         dbmanager.close();

#%%
fout=open('e:/imperial/temp2.txt','w');
for i in range(len(inchikeys)):
    fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(inchikeys[i],formulas[i],exactmasses[i],inchi[i],shortinchi[i],smileslist[i],charges[i]));
fout.close();
#%%
cc=0;
fout=open('e:/imperial/temp3.txt','w');
for i in range(len(inchikeys)):
    if inchi[i]=='':
        fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(inchikeys[i],formulas[i],exactmasses[i],inchi[i],shortinchi[i],smileslist[i],charges[i]));
        cc+=1;
fout.close();
print(cc);