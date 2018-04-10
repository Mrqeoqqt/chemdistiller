import sys;
if sys.byteorder!='little':
    print('Only little endian machines currently supported! bye bye ....');
    quit();

sys.path.append("E:/Imperial/Metaspace_WP4_Source");

from spectra.spectrum import MSSpectrum;
from spectra.mspeak import MSPeak;
from spectralmanager.spectralmanager import SpectralManager;

spectral_manager=SpectralManager();

import os;
#%%

inpath='e:/RawData/1';
onlyfiles = [f for f in os.listdir(inpath) if os.path.isfile(os.path.join(inpath, f))];

dbs=[];

for fname in onlyfiles:
    _,fileext=os.path.splitext(fname);
    if fileext.lower()=='.msp':
        print(fname);
        spectrum=MSSpectrum();
        spectrum.parameters['dbsource']='NIST14';
        finp=open(os.path.join(inpath,fname),'r');
        while True:
            s=finp.readline();
            if s=='':
                break;
            s=s.rstrip('\n');
            if s.startswith('Name:'):
                chemname=s.split('Name: ')[1];
                spectrum.parameters['name']=chemname;
            elif s.startswith('Spectrum_type:'):
                if 'MS5' in s:
                    spectrum.parameters['level']=5;
                elif 'MS4' in s:
                    spectrum.parameters['level']=4;
                elif 'MS3' in s:
                    spectrum.parameters['level']=3;
                elif 'MS2' in s:
                    spectrum.parameters['level']=2;
            elif s.startswith('PrecursorMZ:'):
                if ',' in s:
                    ss=s.split(',');
                    s=ss.pop();
                else:
                    s=s.split('PrecursorMZ: ')[1];
                spectrum.parameters['precursor_mz']=float(s);
            elif s.startswith('Precursor_type:'):
                spectrum.parameters['precursor_ion']=s.split('Precursor_type: ')[1];
            elif s.startswith('Collision_gas: '):
                spectrum.parameters['collision_record']=s.split('Collision_gas: ')[1];
            elif s.startswith('Ion_mode: '):
                if 'P' in s.split('Ion_mode: ')[1]:
                    spectrum.parameters['mode']=1;
                else:
                    spectrum.parameters['mode']=-1;                    
            elif s.startswith('Collision_energy: '):
                try:
                    spectrum.parameters['collision_energy']=float(s.split('Collision_energy: ')[1]);
                except:
                    spectrum.parameters['collision_energy']=-1;
            elif s.startswith('InChIKey: '):
                spectrum.parameters['inchikey']=s.split('InChIKey: ')[1];
            elif s.startswith('Formula: '):
                spectrum.parameters['formula']=s.split('Formula: ')[1];
            elif s.startswith('ExactMass: '):
                spectrum.parameters['exactmass']=float(s.split('ExactMass: ')[1]);
            elif s.startswith('Num Peaks: '):
                npeaks=int(s.split('Num Peaks: ')[1]);
                for i in range(npeaks):
                    s=finp.readline();
                    s=s.rstrip('\n');
                    s=s.split(' ');
                    mz=float(s[0]);
                    intensity=float(s[1]);
                    peak=MSPeak();
                    peak.mz=mz;
                    peak.intensity=intensity;
                    peak.number=i+1;
                    spectrum.add_peak(peak);
                dbs.append(spectrum);
                spectrum=MSSpectrum();
                spectrum.parameters['dbsource']='NIST14';
        finp.close();
        #break;

#%%
for s in dbs:
    if not ('inchikey' in s.parameters.keys()):
        print(s.parameters['name']);

#%%

for l in range(4):
        for i in reversed(range(len(dbs))):
            print(i,l);
            if dbs[i].parameters['level']==5-l:
                dbindex=-1;
                for j in range(len(dbs)):
                    if dbs[j].parameters['level']==4-l:
                        if dbs[i].parameters['mode']==dbs[j].parameters['mode']:
                            for peak in dbs[j].peaks:
                                if abs(peak.mz-dbs[i].parameters['precursor_mz'])<0.1:
                                    peak.add_spectrum(dbs[i]);
                                    dbindex=i;
                                    break;
                    if dbindex<>-1:
                        break;
                if dbindex<>-1:
                    _=dbs.pop(dbindex);
#%%                    
for i in reversed(range(len(dbs))):
            if dbs[i].parameters['level']==4:
                _=dbs.pop(i);
            elif dbs[i].parameters['level']==3:
                _=dbs.pop(i);
            elif dbs[i].parameters['level']==2:
                spectrum=MSSpectrum();
                spectrum.parameters['level']=1;
                spectrum.parameters['mode']=dbs[i].parameters['mode'];
                

                spectrum.parameters['dbsource']=dbs[i].parameters['dbsource'];
                spectrum.parameters['formula']=dbs[i].parameters['formula'];
                spectrum.parameters['exactmass']=dbs[i].parameters['exactmass'];
                #spectrum.parameters['smiles']=dbs[i].parameters['smiles'];
                #spectrum.parameters['inchi']=dbs[i].parameters['inchi'];
                #spectrum.parameters['shortinchi']=dbs[i].parameters['shortinchi'];
                if 'inchikey' in dbs[i].parameters.keys():
                    spectrum.parameters['inchikey']=dbs[i].parameters['inchikey'];

                peak=MSPeak();
                peak.mz=dbs[i].parameters['precursor_mz'];
                peak.intensity=100.0;
                peak.number=1;
                peak.parameters['ion_type']=dbs[i].parameters['precursor_ion'];
                peak.add_spectrum(dbs[i]);
                spectrum.add_peak(peak);
                _=dbs.pop(i);
                dbs.append(spectrum);
            elif dbs[i].parameters['level']==1:
                accepted=False;
                for peak in dbs[i].peaks:
                    if peak.ms_spectra:
                        accepted=True;
                        break;
                if not accepted:
                    _=dbs.pop(i);                    
#%%                
for spectrum in dbs:
        spectral_manager.ms_spectra.append(spectrum);
        print('Current count: %s'%len(spectral_manager.ms_spectra));

        
spectral_manager.export_textfile_spectra_to_folder('e:/Imperial/TestDB/NIST14/');



#%%
inchikeys=[];
formulas=[];
exactmasses=[];
for spectrum in dbs:
    if 'inchikey' in spectrum.parameters.keys():
        inchikey=spectrum.parameters['inchikey'];
        if not (inchikey in inchikeys):
            inchikeys.append(inchikey);
            formulas.append(spectrum.parameters['formula']);
            exactmasses.append(spectrum.parameters['exactmass']);
            print(inchikey);
        
#%%
fout=open('e:/imperial/temp.txt','w')        ;
for i in range(len(inchikeys)):
    fout.write('%s\t%s\t%s\n'%(inchikeys[i],formulas[i],exactmasses[i]));
        
fout.close();
#%%
inchikeys=[];
formulas=[];
exactmasses=[];
inchi=[];
smiles=[];
shortinchi=[];
charges=[];

finp=open('e:/imperial/temp2.txt','r');
for s in finp:
    s=s.rstrip('\n').split('\t');
    inchikeys.append(s[0]);
    formulas.append(s[1]);
    exactmasses.append(float(s[2]));
    inchi.append(s[3]);
    shortinchi.append(s[4]);
    smiles.append(s[5]);
    charges.append(s[6]);
finp.close();
keys={};
for i in range(len(inchikeys)):
    keys[inchikeys[i]]=i;
#%%
for spectrum in dbs:
    if 'inchikey' in spectrum.parameters.keys():
        inchikey=spectrum.parameters['inchikey'];
        if inchikey in keys.keys():
            index=keys[inchikey];
            spectrum.parameters['inchi']=inchi[index];
            spectrum.parameters['shortinchi']=shortinchi[index];
            spectrum.parameters['smiles']=smiles[index];
            spectrum.parameters['charge']=charges[index];
#%%
spectral_manager.export_textfile_spectra_to_folder('e:/Imperial/TestDB/NIST14a/');
#%%
spectral_manager.import_textfile_spectra_from_folder('e:/Imperial/TestDB/HMDB_MassBank');
#%%
spectral_manager.import_textfile_spectra_from_folder('e:/Imperial/TestDB/NIST14a/');

#%%
#Condense and filter spectra here
print(len(spectral_manager.ms_spectra));
for i in reversed(range(len(spectral_manager.ms_spectra))):
    print(i);
    spectrum=spectral_manager.ms_spectra[i];
    if not ('shortinchi' in spectrum.parameters):
        _=spectral_manager.ms_spectra.pop(i);
    elif spectrum.parameters['shortinchi']=='':
        _=spectral_manager.ms_spectra.pop(i);
    if not ('charge' in spectrum.parameters):
        _=spectral_manager.ms_spectra.pop(i);
    elif spectrum.parameters['charge']=='':
        _=spectral_manager.ms_spectra.pop(i);
    if not ('inchi' in spectrum.parameters):
        _=spectral_manager.ms_spectra.pop(i);
    elif spectrum.parameters['inchi']=='':
        _=spectral_manager.ms_spectra.pop(i);
print(len(spectral_manager.ms_spectra));        


for i in reversed(range(len(spectral_manager.ms_spectra))):
    print(i);
    spectrum=spectral_manager.ms_spectra[i];
    for j in range(i):
        #print(i,j);
        prevspectrum=spectral_manager.ms_spectra[j];
        if spectrum.parameters['inchi']==prevspectrum.parameters['inchi']:
            #merge spectra;
            #print('Merging');
            for prevpeak in prevspectrum.peaks:
                #print(prevpeak.mz);
                for peakind in reversed(range(len(spectrum.peaks))):
                    #print(peakind);
                    peak=spectrum.peaks[peakind];
                    if abs(peak.mz-prevpeak.mz)<0.11:
                        for subspectrum in peak.ms_spectra:
                            prevpeak.ms_spectra.append(subspectrum);
                            subspectrum.parent_peak=prevpeak;
                        _=spectrum.peaks.pop(peakind);
                        
            for peak in spectrum.peaks:
                 prevspectrum.peaks.append(peak);
                 peak.parent_spectrum=prevspectrum;
            _=spectral_manager.ms_spectra.pop(i);
            
            break;
            
print(len(spectral_manager.ms_spectra));

sinchis={};
cc=0;
for spectrum in spectral_manager.ms_spectra:
    shortinchi=spectrum.parameters['shortinchi'];
    if shortinchi in sinchis.keys():
        spectrum.parameters['global_index']=sinchis[shortinchi];
    else:
        cc+=1;
        print(cc);
        sinchis[shortinchi]=cc;
        spectrum.parameters['global_index']=cc;


spectral_manager.export_textfile_spectra_to_folder('e:/Imperial/TestDB/Nist14_HMDB_MassBank/');        

#%%
spectral_manager.close();
