# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 14:08:15 2016

@author: ilaponog
"""

inpath='d:/RawData/PhenolExplorer';
import os;
import sys;
import pybel;


if sys.byteorder!='little':
    print('Only little endian machines supported! bye bye ....');
    quit();

def ProcessDBList(inputrecord):
    try:
        mol=pybel.readstring('smi',inputrecord[0]);
        mol.addh();
        mass=mol.exactmass;
        charge=mol.charge;
        if charge!=0:
            mz=abs(mass/charge);
            charged=True;
            if charge>0:
                chargestate='Positive';
            else:
                chargestate='Negative';
        else:
            mz=mass;
            charged=False;
            chargestate='Neutral';
            
        inchi=mol.write('inchi').replace('\n','').replace('InChI=1S/','');
        sinchi=inchi.split('/');
        shortinchi=sinchi[0];
        for j in range(1,len(sinchi)):
                    if sinchi[j][0]=='c' or sinchi[j][0]=='h':
                        shortinchi+='/'+sinchi[j];
        if charged:        
            outstring='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(chargestate,mz,mass,charge,shortinchi,inchi,inputrecord[0],inputrecord[1]);
        else:
            outstring='%s\t%s\t%s\t%s\t%s\t%s'%(chargestate,mz,shortinchi,inchi,inputrecord[0],inputrecord[1]);
        return [0,outstring];
    except:
        return [-1,inputrecord[0],inputrecord[1]];



from glob import glob;

fileslist = [y for x in os.walk(inpath) for y in glob(os.path.join(x[0], '*.csv'))];

for i in range(len(fileslist)):
    fileslist[i]=fileslist[i].replace('\\','/');
print(fileslist);


#%%

for infile in fileslist:
    dirname,fname=os.path.split(infile);
    dirname=dirname.replace(inpath,'');
    print(dirname);
    outpath='d:/Imperial/tempDB/PhenolExplorer'+dirname;
    outfile=outpath+'/totalentries.set';
    outerror=outpath+'/totalentries.error';
    outnames=outpath+'/totalentries.names';

    if not os.path.exists(outpath):
        os.makedirs(outpath);

    cc=0;
    finp=open(infile,'r');
    fout=open(outfile,'w');
    ferr=open(outerror,'w');
    fnames=open(outnames,'w');
    
    dblist=[];
    totcount=0;
    for s in finp:
            s=s.replace('\n','').split(',');
            idx=s[0];
            names=s[4];
            smi=s[1];

        
            totcount+=1;
            if smi!='' and idx!='':
                print(smi,idx,names);
                fnames.write('%s\t%s\n'%(idx,names));
                smi=smi.split('.');
                for ss in smi:
                    dblist.append([ss,idx]);
                cc+=1;
            smi='';
            idx='';
            names='';
            if cc>10000:
                print('Processing....');
                cc=0;
                for db in dblist:
                    i=ProcessDBList(db);
                    if i[0]==0:
                        fout.write('%s\n'%i[1]);
                        #print(i[1]);
                    else:
                        ferr.write('%s\t%s\n'%(i[1],i[2]));
                dblist=[];
        
    if cc>0:
        print('Processing....');
        cc=0;
        for db in dblist:
           i=ProcessDBList(db);
           if i[0]==0:
                    fout.write('%s\n'%i[1]);
                    #print(i[1]);
           else:
                    ferr.write('%s\t%s\n'%(i[1],i[2]));
        dblist=[];
    
    
    finp.close();
    fout.close();
    ferr.close();
    fnames.close();
    #pool.close();
    #pool.join();
    
    print('Finished. Total records: %s'%totcount);


