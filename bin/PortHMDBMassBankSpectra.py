import sys;
import pybel;
if sys.byteorder!='little':
    print('Only little endian machines currently supported! bye bye ....');
    quit();

sys.path.append("E:/Imperial/Metaspace_WP4_Source");

from spectra.spectrum import MSSpectrum;
from spectra.mspeak import MSPeak;
from spectralmanager.spectralmanager import SpectralManager;

spectral_manager=SpectralManager();

#%%
import mysql.connector;
cnx = mysql.connector.connect(user='ivan', password='metaspace', host='155.198.146.41');
cursor = cnx.cursor();
cnx2 = mysql.connector.connect(user='ivan', password='metaspace', host='155.198.146.41');
cursor2 = cnx2.cursor();

workDB="MSX_SPECTRA";

try:
    cnx.database = workDB;  #Setting the Database in the connection.
    cnx2.database = workDB;  #Setting
except:
    print("Cannot set DB: %s" %workDB);


command="Select distinct(FP_ID) \
from Spectra where \
Ac_MassSpecType='MS2' and \
PK_Precision >=1;";
    

res=cursor.execute(command);  
lst=[];
for (FP_ID) in cursor: 
    lst.append(FP_ID[0]);

print(len(lst));


for FP_ID in lst:
    command="Select EntryID, DB_Source, Ch_Formula, Ch_ExactMass, Ch_SMILES, \
    Ch_InChi, Ac_MassSpecType, Ac_MassSpecIonMode, Ac_CollisionEnergy, Ac_CollisionEnergyRecord, MS_FocusedIon, Ac_IonType \
    from Spectra where FP_ID=%s;"%FP_ID;
    dbs=[];    
    res=cursor.execute(command);  
    for EntryID, DB_Source, Ch_Formula, Ch_ExactMass, Ch_SMILES, \
    Ch_InChi, Ac_MassSpecType, Ac_MassSpecIonMode, Ac_CollisionEnergy, \
    Ac_CollisionEnergyRecord, MS_FocusedIon, Ac_IonType in cursor:
        print(EntryID);
        spectrum=MSSpectrum();
        if Ac_MassSpecIonMode=='P':
            spectrum.parameters['mode']=1;
        else:
            spectrum.parameters['mode']=-1;
        if Ac_MassSpecType=='MS' or Ac_MassSpecType=='MS1':
            spectrum.parameters['level']=1;
        elif Ac_MassSpecType=='MS2':
            spectrum.parameters['level']=2;
        elif Ac_MassSpecType=='MS3':
            spectrum.parameters['level']=3;
        elif Ac_MassSpecType=='MS4':
            spectrum.parameters['level']=4;
        mol=pybel.readstring('smi',str(Ch_SMILES));
        mol.addh();
        charge=mol.charge;
        Ch_ExactMass=mol.exactmass;
        spectrum.parameters['dbsource']=DB_Source;
        spectrum.parameters['formula']=Ch_Formula;
        spectrum.parameters['exactmass']=Ch_ExactMass;
        spectrum.parameters['charge']=charge;
        spectrum.parameters['smiles']=Ch_SMILES;
        spectrum.parameters['inchi']=Ch_InChi[3:];
        
        sinchi=Ch_InChi[3:].split('/');
        shortinchi=sinchi[0];
        for j in range(1,len(sinchi)):
                    if sinchi[j][0]=='c' or sinchi[j][0]=='h':
                        shortinchi+='/'+sinchi[j];
        
        spectrum.parameters['shortinchi']=shortinchi;
         
        command="Select PK_No, PK_MZ, PK_Intensity from peaks where EntryID=%s;"%EntryID;
        res=cursor2.execute(command);
        for PK_No, PK_MZ, PK_Intensity in cursor2:
            #print(PK_MZ);
            peak=MSPeak();
            peak.mz=PK_MZ;
            peak.intensity=PK_Intensity;
            peak.number=PK_No;
            spectrum.add_peak(peak);
        if spectrum.parameters['level']>1:
            #if spectrum.parameters['level']==4:
            #    fout.write('')
            spectrum.parameters['collision_energy']=float(Ac_CollisionEnergy);
            spectrum.parameters['collision_record']=Ac_CollisionEnergyRecord;
            precursor_mz=0.0;
            precursor_ion='';
            if ('PRECURSOR_M/Z' in MS_FocusedIon) and ('PRECURSOR_TYPE' in MS_FocusedIon):
                MS_FocusedIon=MS_FocusedIon.replace('\t',' ').rstrip('\n').replace('  ',' ').split(' ');
                #print(MS_FocusedIon);
                for i in range(len(MS_FocusedIon)):
                    if MS_FocusedIon[i]=='PRECURSOR_M/Z':
                        #print(MS_FocusedIon[i+1]);
                        #print(bytes(MS_FocusedIon[i+1]));
                        precursor_mz=MS_FocusedIon[i+1];
                    if MS_FocusedIon[i]=='PRECURSOR_TYPE':
                        precursor_ion=MS_FocusedIon[i+1];
                if precursor_ion<>'' and precursor_mz>0.0:
                    if not (('*' in precursor_ion) or ('/' in precursor_ion) or ('*' in precursor_mz)):
                        if '/' in precursor_mz:
                            precursor_mz=precursor_mz.split('/');
                            precursor_mz=precursor_mz.pop();
                        spectrum.parameters['precursor_mz']=float(precursor_mz);
                        spectrum.parameters['precursor_ion']=precursor_ion;
                        dbs.append(spectrum);
        else:
            dbs.append(spectrum);

    #processdbs;
    for l in range(3):
        for i in reversed(range(len(dbs))):
            if dbs[i].parameters['level']==4-l:
                dbindex=-1;
                for j in range(len(dbs)):
                    if dbs[j].parameters['level']==3-l:
                        if dbs[i].parameters['mode']==dbs[j].parameters['mode']:
                            for peak in dbs[j].peaks:
                                if abs(peak.mz-dbs[i].parameters['precursor_mz'])<0.11:
                                    peak.add_spectrum(dbs[i]);
                                    dbindex=i;
                                    break;
                    if dbindex<>-1:
                        break;
                if dbindex<>-1:
                    _=dbs.pop(dbindex);
                    
    for i in reversed(range(len(dbs))):
            if dbs[i].parameters['level']==4:
                _=dbs.pop(i);
                #print('level4');
            elif dbs[i].parameters['level']==3:
                _=dbs.pop(i);
            elif dbs[i].parameters['level']==2:
                spectrum=MSSpectrum();
                spectrum.parameters['level']=1;
                spectrum.parameters['mode']=dbs[i].parameters['mode'];
                

                spectrum.parameters['dbsource']=dbs[i].parameters['dbsource'];
                spectrum.parameters['formula']=dbs[i].parameters['formula'];
                spectrum.parameters['exactmass']=dbs[i].parameters['exactmass'];
                spectrum.parameters['charge']=dbs[i].parameters['charge'];
                spectrum.parameters['smiles']=dbs[i].parameters['smiles'];
                spectrum.parameters['inchi']=dbs[i].parameters['inchi'];
                spectrum.parameters['shortinchi']=dbs[i].parameters['shortinchi'];

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
                
    for spectrum in dbs:
        spectral_manager.ms_spectra.append(spectrum);
    print('Current count: %s'%len(spectral_manager.ms_spectra));
    #if len(spectral_manager.ms_spectra)>200:
    #    break;
        
#%%

spectral_manager.export_textfile_spectra_to_folder('e:/Imperial/TestDB/HMDB_MassBank/');
spectral_manager.close();

cursor2.close();
cursor.close();  #Tidying up at the end.
cnx.close();    
cnx2.close();    






