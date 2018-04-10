# -*- coding: utf-8 -*-
"""
Created on Mon Aug 08 12:24:51 2016

@author: ilaponog
"""
#%%

deletedcount=0;
import os;
from glob import glob;
  
    
def ProcessDBEntryFile(fnames):
    #global lock;
    
    
    dblist=[];
    fin=open(fnames[0],'r');
    for s in fin:
        dblist.append(s.replace('\n','').split('\t'));
    fin.close();
    idx=[];
    if fnames[3]:
        offs=2;
    else :
        offs=0;
    for i in dblist:
        idx.append(i[4+offs]);
    deleted=[];
    deletedidx=[];
    for ss in dblist:
        if ss[1+offs]=='' or ss[2+offs]=='' or ss[3+offs]=='' or ss[4+offs]=='' or ss[5+offs]=='' or ss[6+offs]=='':
            deletedidx.append(ss[4+offs]);
        
    for midx in deletedidx:
        while midx in idx:
            i=idx.index(midx);
            deleted.append(dblist.pop(i));
            idx.pop(i);
    if len(dblist)>0:
        dirname, filename = os.path.split(fnames[1]);
        lock.acquire();
        if not os.path.exists(dirname):
           os.makedirs(dirname);            
        lock.release();
        fout=open(fnames[1],'w');
        for ss in dblist:
            fout.write(ss[0]);
            for i in range(1,len(ss)):
                fout.write('\t%s'%ss[i]);
            fout.write('\n');
        fout.close();
    if len(deleted)>0:
        fout=open(fnames[2],'w');
        for ss in deleted:
            fout.write(ss[0]);
            for i in range(1,len(ss)):
                fout.write('\t%s'%ss[i]);
            fout.write('\n');
            if ss[1+offs]=='':
                fout.write('Missing Inchi\n');
            if ss[2+offs]=='':
                fout.write('Missing InchiStereo\n');
            if ss[3+offs]=='':
                fout.write('Missing SMILES\n');
            if ss[4+offs]=='':
                fout.write('Missing ID\n');
            if ss[5+offs]=='':
                fout.write('Missing Fingerprint\n');
            if ss[6+offs]=='':
                fout.write('Missing Fragments\n');
                
        fout.close();
        filename,fileext=os.path.splitext(fnames[2]);
        fout=open(filename+'.smi','w');
        cid='';
        cid_count=1;
        for ss in deleted:
            if ss[4+offs]!=cid:
                cid=ss[4+offs];
                cid_count=1;
                outcid=cid;
            else:
                cid_count+=1;
                outcid='%s_%s'%(cid,cid_count);
            fout.write('%s %s_%s\n'%(ss[3+offs],fnames[4],outcid))
        fout.close();
                               
    return [fnames[0],len(deleted)];
     

def ProcessSubRaw(inpath,outpath,deletedpath,subpath,pool,charged=False):
    fileslist = [y for x in os.walk(inpath+'/'+subpath) for y in glob(os.path.join(x[0], '*.st2'))];
    flist=[];
    removed=0;
    for fname in fileslist:
        fname=fname.replace('\\','/');
        dirname, filename = os.path.split(fname);
        outfname=outpath+'/'+fname.replace(inpath+'/','');
        missfname=deletedpath+'/'+fname.replace(inpath+'/','').replace('/','_');
        masspath=fname.replace(inpath+'/','');
        masspath,ext=os.path.splitext(masspath);
        masspath=masspath.replace('/','_');
        
        flist.append([fname,outfname,missfname,charged,masspath]);
    for i in pool.imap_unordered(ProcessDBEntryFile, flist):
        print('In %s deleted %s entries.'%(i[0],i[1]));
        removed+=i[1];
    return removed;
    
    
from multiprocessing import Pool,Lock;
import time;

def LoadDBRecord(inpath):
    dbrecord={'FORMAT':'RAW','TIMESTAMP':'','DBNAME':'','DBPARENTNAME':'','DBPATH':'','DBPARENTPATH':'','DBVERSION':'0','DBPARENTVERSION':'','DBPARENTTIMESTAMP':'','DBLASTOPERATION':'','DBREMARK':'','DBPARENTREMARK':''};
    if os.path.isfile(inpath+'/dbinfo.txt'):
        finp=open(inpath+'/dbinfo.txt','r');
        for s in finp:
            s=s.replace('\n','');
            i=s.index('=');
            s1=s[:i].upper();
            s2=s[i+1:];
            if s1 in dbrecord:
                dbrecord[s1]=s2;
        finp.close();
    return dbrecord;

def GetNewDBRecord(newdbname,outpath,operationdesc,remarkrec):
    newdbrecord={'FORMAT':'RAW','DBPARENTNAME':'','DBPARENTPATH':'','DBPARENTVERSION':'','DBPARENTTIMESTAMP':'','DBLASTPARENTOPERATION':'','DBPARENTREMARK':''};
    newdbrecord['TIMESTAMP']=time.strftime("%H:%M:%S %d/%m/%Y");
    newdbrecord['DBNAME']=newdbname;
    newdbrecord['DBPATH']=outpath;
    newdbrecord['DBVERSION']=1;
    newdbrecord['DBLASTOPERATION']=operationdesc;
    newdbrecord['DBREMARK']=remarkrec;
    return newdbrecord;
    

def GetNextDBRecord(dbrecord,newdbname,outpath,operationdesc,remarkrec):
    newdbrecord={'FORMAT':'RAW','TIMESTAMP':'','DBNAME':'','DBPARENTNAME':'','DBPATH':'','DBPARENTPATH':'','DBVERSION':'0','DBPARENTVERSION':'','DBPARENTTIMESTAMP':'','DBLASTOPERATION':'','DBREMARK':'','DBPARENTREMARK':''};
    newdbrecord['DBPARENTTIMESTAMP']=dbrecord['TIMESTAMP'];
    newdbrecord['TIMESTAMP']=time.strftime("%H:%M:%S %d/%m/%Y");
    
    newdbrecord['DBPARENTNAME']=dbrecord['DBNAME'];
    newdbrecord['DBNAME']=newdbname;
    
    newdbrecord['DBPARENTPATH']=dbrecord['DBPATH'];
    newdbrecord['DBPATH']=outpath;
    
    newdbrecord['DBPARENTVERSION']=dbrecord['DBVERSION'];
    newdbrecord['DBVERSION']=int(dbrecord['DBVERSION'])+1;
    
    newdbrecord['DBLASTPARENTOPERATION']=dbrecord['DBLASTOPERATION'];
    newdbrecord['DBLASTOPERATION']=operationdesc;
    
    newdbrecord['DBPARENTREMARK']=dbrecord['DBREMARK'];
    newdbrecord['DBREMARK']=remarkrec;
    
    return newdbrecord;

def SaveDBRecord(outpath,dbrecord):
    fout=open(outpath+'/dbinfo.txt','w');
    for key in dbrecord.keys():
        fout.write('%s=%s\n'%(key,dbrecord[key]));
    fout.close();
    
def initgloballock(l):
    global lock;
    lock = l;

def PurgeDBRaw(inpath,outpath,deletedpath,newdbname,dbremark,ncpu=1):
    if not os.path.exists(inpath):
        raise NameError('No source !');
    dbrecord=LoadDBRecord(inpath);
    if dbrecord['FORMAT']!='RAW':
        raise NameError('DB must be in RAW format!');
    if not os.path.exists(outpath):
        os.makedirs(outpath);
    if not os.path.exists(deletedpath):
        os.makedirs(deletedpath);
    totalremoved=0;    
    lock = Lock();
    pool = Pool(initializer=initgloballock, initargs=(lock,),processes=ncpu);
    #pool=Pool(ncpu);
    
    print('Processing Neutral');
    totalremoved+=ProcessSubRaw(inpath,outpath,deletedpath,'Neutral',pool,False);
    print('Processing Negative');
    totalremoved+=ProcessSubRaw(inpath,outpath,deletedpath,'Negative',pool,True);
    print('Processing Positive');
    totalremoved+=ProcessSubRaw(inpath,outpath,deletedpath,'Positive',pool,True);
    print('Deleted in total: %s'%totalremoved);
    fileslist = [y for x in os.walk(deletedpath) for y in glob(os.path.join(x[0], '*.smi'))];
    if len(fileslist)>0:
        
        fout=open(deletedpath+'/TotalDeleted.smi','w');
        for s in fileslist:
            s=s.replace('\\','/');
            if s!=deletedpath+'/TotalDeleted.smi':
                finp=open(s,'r');
                for ss in finp:
                    fout.write(ss);
                finp.close();
        fout.close();
    
    newdbrecord=GetNextDBRecord(dbrecord,newdbname,outpath,'Deleted compounds with missing records','Deleted records can be found in subfolder: %s'%deletedpath);
    
    SaveDBRecord(outpath,newdbrecord);
    pool.close();
    pool.join();

    
if __name__=='__main__':
    subpath='PubChem';
    inpath='e:/Imperial/tempDB2/'+subpath;
    outpath='e:/Imperial/tempDB3/'+subpath;
    deletedpath='e:/Imperial/tempDB3/Deleted/'+subpath;
    PurgeDBRaw(inpath,outpath,deletedpath,subpath,'Clean',4);

  