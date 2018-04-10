# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 21:00:31 2016

@author: ilaponog
"""

infile='e:/RawData/FooDB/compounds11.csv';
outfile='e:/RawData/FooDB/compounds.txt';

#%%
dblist=[];
s='';
inqu=False;
ss=[];

finp=open(infile,'rb');
for sc in finp:
  #print(sc);
  #sc+='\n';
  #break;  
  sc=sc.replace('\r','').replace('\\"','');
  for c in sc:
    if c=='"':
        inqu=not inqu;
    
    if c=='\n' and not inqu:
        ss.append(s);
        dblist.append(ss);
        print(ss[0]);
        if ss[0]!='id':
            int(ss[0])
        ss=[];
        s='';
        #break;
    elif c==',' and not inqu:
        ss.append(s);
        s='';
    else:
        s+=c;
finp.close();

#%%

#infile='e:/RawData/FooDB/compounds11.csv';
#outfile='e:/RawData/FooDB/compounds_nor.csv';

#finp=open(infile,'r');
fout=open(outfile,'w');
for db in dblist:
    fout.write('%s\t%s\t%s\n'%(db[0].replace('"','').replace('"',''),db[len(db)-19].replace('"',''),db[4].replace('"','')));
#for s in finp:
#    s=s.replace('\r','');
#    fout.write(s);

#finp.close();    
fout.close();    
    
    
    
    
    
    