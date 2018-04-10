# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 14:33:47 2016

@author: Dr. Ivan Laponogov, ICL, London
"""
import math;

BondDistanceDictionary=[
["H",	1,	"H",	432,	0.74,	0,],
["H",	1,	"B",	389,	1.19,	0,],
["H",	1,	"C",	411,	1.09,	0,],
["H",	1,	"Si",	318,	1.48,	0,],
["H",	1,	"Ge",	288,	1.53,	0,],
["H",	1,	"Sn",	251,	1.70,	0,],
["H",	1,	"N",	386,	1.01,	0,],
["H",	1,	"P",	322,	1.44,	0,],
["H",	1,	"As",	247,	1.52,	0,],
["H",	1,	"O",	459,	0.96,	0,],
["H",	1,	"S",	363,	1.34,	0,],
["H",	1,	"Se",	276,	1.46,	0,],
["H",	1,	"Te",	238,	1.70,	0,],
["H",	1,	"F",	565,	0.92,	0,],
["H",	1,	"Cl",	428,	1.27,	0,],
["H",	1,	"Br",	362,	1.41,	0,],
["H",	1,	"I",	295,	1.61,	0,],
["B",	1,	"B",	293,	1.75,	1,],
["B",	1,	"O",	536,	1.75,	1,],
["B",	1,	"F",	613,	1.75,	1,],
["B",	1,	"Cl",	456,	1.75,	0,],
["B",	1,	"Br",	377,	1.75,	1,],
["C",	1,	"C",	346,	1.54,	0,],
["C",	2,	"C",	602,	1.34,	0,],
["C",	3,	"C",	835,	1.20,	0,],
["C",	4,	"C",	518,	1.40,	0,],
["C",	1,	"Si",	318,	1.85,	0,],
["C",	1,	"Ge",	238,	1.95,	0,],
["C",	1,	"Sn",	192,	2.16,	0,],
["C",	1,	"Pb",	130,	2.30,	0,],
["C",	1,	"N",	305,	1.47,	0,],
["C",	2,	"N",	615,	1.29,	0,],
["C",	3,	"N",	887,	1.16,	0,],
["C",	4,	"N",	518,	1.37,	1,],
["C",	1,	"P",	264,	1.84,	0,],
["C",	4,	"P",	500,	1.73,	1,],
["C",	1,	"O",	358,	1.43,	0,],
["C",	2,	"O",	799,	1.20,	0,],
["C",	3,	"O",	1072,	1.13,	0,],
["C",	4,	"O",	500,	1.30,	1,],
["C",	1,	"B",	356,	1.70,	1,],
["C",	1,	"S",	272,	1.82,	0,],
["C",	2,	"S",	573,	1.60,	0,],
["C",	1,	"F",	485,	1.35,	0,],
["C",	1,	"Cl",	327,	1.77,	0,],
["C",	1,	"Br",	285,	1.94,	0,],
["C",	1,	"I",	213,	2.14,	0,],
["Si",	1,	"Si",	222,	2.33,	0,],
["Si",	1,	"N",	355,	1.70,	1,],
["Si",	1,	"O",	452,	1.63,	0,],
["Si",	1,	"S",	293,	2.00,	0,],
["Si",	1,	"F",	565,	1.60,	0,],
["Si",	1,	"Cl",	381,	2.02,	0,],
["Si",	1,	"Br",	310,	2.15,	0,],
["Si",	1,	"I",	234,	2.43,	0,],
["Ge",	1,	"Ge",	188,	2.41,	0,],
["Ge",	1,	"N",	257,	2.30,	1,],
["Ge",	1,	"F",	470,	1.68,	0,],
["Ge",	1,	"Cl",	349,	2.10,	0,],
["Ge",	1,	"Br",	276,	2.30,	0,],
["Ge",	1,	"I",	212,	2.30,	1,],
["Sn",	1,	"F",	414,	2.33,	1,],
["Sn",	1,	"Cl",	323,	2.33,	0,],
["Sn",	1,	"Br",	273,	2.50,	0,],
["Sn",	1,	"I",	205,	2.70,	0,],
["Pb",	1,	"F",	331,	2.70,	1,],
["Pb",	1,	"Cl",	243,	2.42,	0,],
["Pb",	1,	"Br",	201,	2.42,	1,],
["Pb",	1,	"I",	142,	2.79,	0,],
["N",	1,	"N",	167,	1.45,	0,],
["N",	2,	"N",	418,	1.25,	0,],
["N",	3,	"N",	942,	1.10,	0,],
["N",	1,	"O",	201,	1.40,	0,],
["N",	1,	"S",	201,	1.72,	1,],
["N",	2,	"O",	607,	1.21,	0,],
["N",	1,	"F",	283,	1.36,	0,],
["N",	1,	"Cl",	313,	1.75,	0,],
["P",	1,	"P",	201,	2.21,	0,],
["P",	1,	"O",	335,	1.63,	0,],
["P",	2,	"O",	544,	1.50,	0,],
["P",	2,	"S",	335,	1.86,	0,],
["P",	1,	"F",	490,	1.54,	0,],
["P",	1,	"Cl",	326,	2.03,	0,],
["P",	1,	"Br",	264,	2.03,	1,],
["P",	1,	"I",	184,	2.03,	1,],
["As",	1,	"As",	146,	2.43,	0,],
["As",	1,	"O",	301,	1.78,	0,],
["As",	1,	"F",	484,	1.71,	0,],
["As",	1,	"Cl",	322,	2.16,	0,],
["As",	1,	"Br",	458,	2.33,	0,],
["As",	1,	"I",	200,	2.54,	0,],
["Sb",	1,	"Sb",	121,	2.32,	1,],
["Sb",	1,	"F",	440,	2.32,	1,],
["Sb",	1,	"Cl",	248,	2.32,	1,],
["Sb",	1,	"Cl",	315,	2.32,	0,],
["O",	1,	"O",	142,	1.48,	0,],
["O",	2,	"O",	494,	1.21,	0,],
["O",	1,	"F",	190,	1.42,	0,],
["S",	2,	"O",	522,	1.43,	0,],
["S",	1,	"S",	226,	2.05,	0,],
["S",	2,	"S",	425,	1.49,	0,],
["S",	1,	"F",	284,	1.56,	0,],
["S",	1,	"Cl",	255,	2.07,	0,],
["Se",	1,	"Se",	172,	2.32,	1,],
["Se",	2,	"Se",	272,	2.15,	0,],
["F",	1,	"F",	155,	1.42,	0,],
["Cl",	1,	"Cl",	240,	1.99,	0,],
["Br",	1,	"Br",	190,	2.28,	0,],
["I",	1,	"I",	148,	2.67,	0,],
["I",	1,	"O",	201,	2.00,	1,],
["I",	1,	"F",	273,	1.91,	0,],
["I",	1,	"Cl",	208,	2.32,	0,],
["I",	1,	"Br",	175,	2.30,	1,],
["Kr",	1,	"F",	50,	1.90,	0,],
["Xe",	1,	"O",	84,	1.75,	0,],
["Xe",	1,	"F",	130,	1.95,	0,],
];

def PredictBondsFromXYZ(atoms,atomsXYZ):
    global BondDistanceDictionary;
    bonds=[];
    if len(atoms)==0 or len(atoms)!=len(atomsXYZ):
        raise NameError('Wrong input!');
    for i in range(0,len(atoms)-1):
        c1=atoms[i][1];
        x1,y1,z1=atomsXYZ[i];
        for j in range(i+1,len(atoms)):
            c2=atoms[j][1];
            x2,y2,z2=atomsXYZ[j];
            dx=x1-x2;dy=y1-y2;dz=z1-z2;
            dd=math.sqrt(dx*dx+dy*dy+dz*dz);
            if dd<3.2:
                candid=[];
                for k in range(0,len(BondDistanceDictionary)):
                    a1=BondDistanceDictionary[k][0];
                    a2=BondDistanceDictionary[k][2];
                    if (a1==c1 and a2==c2) or (a1==c2 and a2==c1):
                        candid.append(BondDistanceDictionary[k]);
                mindd=5.0;
                enn=0;
                idd=0;
                    
                for k in range(0,len(candid)):
                    ddd=abs(dd-candid[k][4])/candid[k][4];
                    if ddd<0.15:
                        if mindd>ddd:
                            idd=candid[k][1];
                            mindd=ddd;
                            enn=candid[k][3];
                if idd>0:
                    bonds.append([i,j,idd,1.0/enn]);
                    
    
    return bonds;
    


def GetUniqueBonds(bonds):
        uniquebonds=[];
        for i in bonds:
            if len(i)==2:
                i.append(1);
                i.append(0.0);
            elif len(i)==3:
                i.append(0.0);
            if not (([i[0],i[1],i[2],i[3]] in uniquebonds) or ([i[1],i[0],i[2],i[3]] in uniquebonds)):
                uniquebonds.append(i);
        return uniquebonds;
    

class FragmenterBasic:
    
    def __PropagateAtomGroup(self,index,groupindex):
        for i in range(0,self.AtomCount):
            if self.BondsArray[i][index]>0:
                if self.AtomGroups[i]==0:
                    self.AtomGroups[i]=groupindex;
                    self.__PropagateAtomGroup(i,groupindex);
    

    def GenerateSDFFromAtomNumbers(self,AtomNumbers):
        SDF='AtomSet %s\nBasicFingerprinter\n\n'%AtomNumbers;
        AtomBonds=[];
        for i in range(0,len(AtomNumbers)-1):
            for j in range(i+1,len(AtomNumbers)):
                if self.BondsArray[AtomNumbers[i]][AtomNumbers[j]]>0:
                    AtomBonds.append([i,j,self.BondsArray[AtomNumbers[i]][AtomNumbers[j]]]);
        SDF+='%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d V2000\n'%(len(AtomNumbers),len(AtomBonds),0,0,0,0,0,0,0,0,999);
        atoms=[];
        for i in AtomNumbers:
            atoms.append([i,self.Atoms[i],0.0,0.0,0.0]);
        
        if hasattr(self,'AtomXYZ'):
            for atom in atoms:
                atom[2]=self.AtomXYZ[atom[0]][0];
                atom[3]=self.AtomXYZ[atom[0]][1];
                atom[4]=self.AtomXYZ[atom[0]][2];
            
        for atom in atoms:
            charge=self.Atoms[atom[0]][3];
            if charge==3:
                ccc=1;
            elif charge==2:
                ccc=2;
            elif charge==1:
                ccc=3;
            elif charge==-1:
                ccc=5;
            elif charge==-2:
                ccc=6;
            elif charge==-3:
                ccc=7;
            else:
                ccc=0;
            v=0;
            for j in AtomNumbers:
                if self.BondsArray[atom[0]][j]>0:
                    v+=1;
            SDF+='%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d\n'%(atom[2],atom[3],atom[4],\
            atom[1][1],atom[1][2],\
            ccc,\
            0,1,0,v,0,0,0,0,0,0);
        for bond in AtomBonds:
            SDF+='%3d%3d%3d%3d%3d%3d%3d\n'%(bond[0]+1,bond[1]+1,bond[2],0,0,0,0);
        charge=0.0;
        for i in AtomNumbers:
            charge+=self.Atoms[i][3];
        SDF+='>  <TOTAL_CHARGE>\n%s\n'%charge;
        mass=0.0;
        for i in AtomNumbers:
            mass+=self.Atoms[i][0];
        SDF+='>  <EXACT_MASS>\n%s\n'%mass;
        SDF+='$$$$';
        return SDF;        

    def __GenerateGroups(self,current_level):
        for i in range(0,self.AtomCount):
            self.AtomGroups[i]=0;
        #self.AtomGroups=[0]*self.AtomCount;
        groupcount=0;
        for i in range(0,self.AtomCount):
            if self.AtomGroups[i]==0:
                groupcount+=1;
                self.AtomGroups[i]=groupcount;
                self.__PropagateAtomGroup(i,groupcount);
                
        for i in range(1, groupcount+1):
            mass=0.0;
            for j in range(0, self.AtomCount):
                if self.AtomGroups[j]==i:
                    mass+=self.Atoms[j][0];
            if not (mass in self.Fragments):
                self.Fragments.append(mass);
                if self.ReturnAtomNumbers:
                    self.GroupAtomNumbers.append([]);
                    FragmentIndex=len(self.Fragments)-1;
            else:
                if self.ReturnAtomNumbers:
                    FragmentIndex=self.Fragments.index(mass);
            if self.ReturnAtomNumbers:
                    numbers=[];
                    for j in range(0, self.AtomCount):
                        if self.AtomGroups[j]==i:
                            numbers.append(j);
                    if self.ReturnAtomNumbers and (not (numbers in self.GroupAtomNumbers[FragmentIndex])):
                        self.GroupAtomNumbers[FragmentIndex].append(numbers);
                    
        

    def __ProcessFragment(self,current_level):
        for i in range(0,len(self.BreakableBonds)):
            if self.BreakableBonds[i][3]==True:
                bond=self.BreakableBonds[i][1];
                self.BondsArray[bond[0]][bond[1]]=0;
                self.BondsArray[bond[1]][bond[0]]=0;
        self.__GenerateGroups(current_level);
        for i in range(0,len(self.BreakableBonds)):
            if self.BreakableBonds[i][3]==True:
                bond=self.BreakableBonds[i][1];
                if len(bond)>2:
                    self.BondsArray[bond[0]][bond[1]]=bond[2];
                    self.BondsArray[bond[1]][bond[0]]=bond[2];
                else:
                    self.BondsArray[bond[0]][bond[1]]=1;
                    self.BondsArray[bond[1]][bond[0]]=1;
            
        
    def __ProcessFragmentationLevel(self,fragmentation_depth,current_level):
        for i in range(0,len(self.BreakableBonds)):
            if self.BreakableBonds[i][3]==False:
                    self.BreakableBonds[i][3]=True;
                    self.BrokenBonds.append(self.BreakableBonds[i]);
                    self.__ProcessFragment(current_level);
                    
                    if current_level<fragmentation_depth:
                        self.__ProcessFragmentationLevel(fragmentation_depth,current_level+1);

                    self.BreakableBonds[i][3]=False;
                    #print(self.BrokenBonds);
                    self.BrokenBonds.pop();
    
    
    
    def SortResults(self,Ascending=True):
            if Ascending:
                self.Fragments=sorted(self.Fragments,key=lambda sortit: sortit[0]);
            else:
                self.Fragments=sorted(self.Fragments,key=lambda sortit: sortit[0], reverse=True);
            
        
        
        
    def __init__(self,atoms,bonds,atomsXYZ=[],breakage_threshold=0.0,breakage_threshold_absolute=False,fragmentation_depth=2,Group_Atom_Numbers=False,forces=[],masspower=0, use_top_energy=False, top_energy_count=1):
        '''
        Preparation of the fragmentation inputs
        '''
        if len(atomsXYZ)>0 and len(atomsXYZ)==len(atoms):
            self.AtomXYZ=atomsXYZ;
        elif len(atomsXYZ)>0:
            raise NameError('Atoms and AtomsXYZ have different lengths!');
        if len(atoms)==0: 
            raise NameError('Empty Atoms!');
        self.ReturnAtomNumbers=Group_Atom_Numbers;
        self.AtomCount=len(atoms);
        self.Atoms=atoms;
        self.AtomGroups=[0]*self.AtomCount;
        self.BondsArray= [[0] * self.AtomCount for i in range(self.AtomCount)];
        self.Bonds=bonds;
        self.BreakableBonds=[];
        self.BrokenBonds=[];
        self.BondForces=[];
        for bond in self.Bonds:
            if len(bond)>2:
                self.BondsArray[bond[0]][bond[1]]=bond[2];
                self.BondsArray[bond[1]][bond[0]]=bond[2];
            else:
                self.BondsArray[bond[0]][bond[1]]=1;
                self.BondsArray[bond[1]][bond[0]]=1;
            
            if len(bond)>3:
                self.BondForces.append(bond[3]);
            else:
                self.BondForces.append(0.0);
        self.Forces=forces;
        if len(forces)==len(atoms):
            print('Making weighted bond forces');
            for i  in range(0,len(self.Bonds)):
                bond=self.Bonds[i];
                f1=forces[bond[0]];
                f2=forces[bond[1]];
                m1=self.Atoms[bond[0]][0];
                m2=self.Atoms[bond[1]][0];
                m1=math.pow(m1,masspower);
                m2=math.pow(m2,masspower);
                self.BondForces[i]=math.sqrt(f1*f1*m1+f2*f2*m2);
                
        if breakage_threshold_absolute==False:
            minf=1.0e8;
            maxf=0.0;
            for bondforce in self.BondForces:
                if bondforce>maxf:
                    maxf=bondforce;
                if bondforce<minf:
                    minf=bondforce;
            ForceRange=maxf-minf;
            if abs(ForceRange)<1.0e-9:
                print('Forces Identical! Cannot normalise! Setting all to 0.0');
                for i in range(0,len(self.BondForces)):
                    self.BondForces[i]=0.0;
            else:
                for i in range(0,len(self.BondForces)):
                    self.BondForces[i]=(self.BondForces[i]-minf)/ForceRange;
        
        if use_top_energy:
            bf=sorted(self.BondForces,reverse=True);
            if top_energy_count>len(bf):
                breakage_threshold=bf[len(bf)-1];
            else:
                breakage_threshold=bf[top_energy_count-1];
                
        self.BreakageThreshold=breakage_threshold;
        for i in range(0,len(self.BondForces)):
            if self.BondForces[i]>=breakage_threshold:
                self.BreakableBonds.append([i,bonds[i],self.BondForces[i],False]);
        
        
        fragmentation_depth=min(fragmentation_depth,len(self.BreakableBonds));
        #print('FragDepth %s'%fragmentation_depth);
        self.Fragments=[];
        if self.ReturnAtomNumbers:
            self.GroupAtomNumbers=[];
        self.__ProcessFragment(0);
        if fragmentation_depth>0:
            self.__ProcessFragmentationLevel(fragmentation_depth,1);
        for i in range(0,len(self.Fragments)):
            self.Fragments[i]=[self.Fragments[i]];
        if self.ReturnAtomNumbers:
            for i in range(0,len(self.Fragments)):
                self.Fragments[i].append(self.GroupAtomNumbers[i]);
            
                



if __name__=='__main__':

    '''
    import timeit;
     
    print(timeit.timeit("atoms=[[12.0,'C',0,0,0], [12.0,'C',0,0,0], [12.0,'C',0,0,0], [12.0,'C',0,0,0], [12.0,'C',0,0,0], [12.0,'C',0,0,0], [14.003074,'N',0,0,0], [12.0,'C',0,0,0], [14.003074,'N',0,0,0], [14.003074,'N',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0]];\
    bonds=[[0, 1, 4], [0, 5, 4], [0, 11, 1], [1, 2, 4], [1, 12, 1], [2, 3, 4], [2, 13, 1], [3, 4, 4], [3, 9, 1], [4, 5, 4], [4, 6, 1], [5, 14, 1], [6, 7, 1], [6, 17, 1], [6, 18, 1], [7, 8, 1], [7, 9, 1], [8, 15, 1], [8, 16, 1], [9, 10, 1]];\
    frag=FragmenterBasic(atoms,bonds,breakage_threshold=0,\
    breakage_threshold_absolute=True,\
    fragmentation_depth=2,\
    Group_Atom_Numbers=True);", \
    setup="from __main__ import FragmenterBasic",number=100));

    '''

    atomsXYZ=[\
    [0.04945500,   -0.11684300,    0.15115200],\
    [1.41715100,   -0.09888600,    0.15645200],\
    [2.16194500,    0.05035500,    1.32240700],\
    [1.43010800,    0.11899500,    2.50812800],\
    [0.00591200,    0.12587300,    2.49966400],\
    [-0.71147800,    0.03869400,    1.30628400],\
    [-0.41136200,    0.14424700,    3.80919700],\
    [0.71067400,   0.01790900,    4.63422000],\
    [0.78900500,   -0.22026900,    5.93421500],\
    [1.81935900,    0.16876000,    3.82026600],\
    [2.77095700,   -0.00156300,    4.17986200],\
    [-0.51240000,   -0.20911300,   -0.78243100],\
    [1.97578000,   -0.13222500,   -0.78170600],\
    [3.23306900,    0.10908800,    1.35820900],\
    [-1.77494500,    0.14621100,    1.32983300],\
    [0.00939887,   -0.77833871,    6.21841999],\
    [1.63866877,   -0.70868697,    6.13300979],\
    [-1.36890700,    0.16187100,    4.15897600]];

    atoms=[[12.0,'C',0,0,0], [12.0,'C',0,0,0], [12.0,'C',0,0,0], [12.0,'C',0,0,0], [12.0,'C',0,0,0], [12.0,'C',0,0,0], [14.003074,'N',0,0,0], [12.0,'C',0,0,0], [14.003074,'N',0,0,0], [14.003074,'N',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0], [1.007825,'H',0,0,0]];
    bonds=[[0, 1, 4], [0, 5, 4], [0, 11, 1], [1, 2, 4], [1, 12, 1], [2, 3, 4], [2, 13, 1], [3, 4, 4], [3, 9, 1], [4, 5, 4], [4, 6, 1], [5, 14, 1], [6, 7, 1], [6, 17, 1], [7, 8, 1], [7, 9, 1], [8, 15, 1], [8, 16, 1], [9, 10, 1]];
    
    frag=FragmenterBasic(atoms,bonds,atomsXYZ,breakage_threshold=0.0,breakage_threshold_absolute=True,fragmentation_depth=2,Group_Atom_Numbers=True);
    frag.SortResults(Ascending=False);
    print(frag.Fragments);
    if frag.ReturnAtomNumbers:
        print(frag.GenerateSDFFromAtomNumbers(frag.Fragments[0][1][0]));





