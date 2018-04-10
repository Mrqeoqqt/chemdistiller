# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 17:56:58 2016

@author: ilaponog
"""


class Fragmenter:
    def GetCharges(self,AtomGroups):
        result=[];
        for group in AtomGroups:
            charge=0;
            for i in group:
                charge+=self.Atoms[i][1];
            if not (charge in result):
                result.append(charge);
        return result;
        
    
    def __PropagateAtomGroup(self,index,groupindex):
        recd=self.Atoms[index][2];
        for i in range(0,len(recd)):
            if recd[i][1]==True:
                j=recd[i][0];
                if self.AtomGroups[j]==0:
                    self.AtomGroups[j]=groupindex;
                    self.__PropagateAtomGroup(j,groupindex);
    

    def __GenerateGroups(self,current_level):
        for i in range(0,self.AtomCount):
            self.AtomGroups[i]=0;
        groupcount=0;
        for i in range(0,self.AtomCount):
            if self.AtomGroups[i]==0:
                groupcount+=1;
                self.AtomGroups[i]=groupcount;
                self.__PropagateAtomGroup(i,groupcount);
                
        mass=[0.0]*groupcount;
        
        for j in range(0, self.AtomCount):
            mass[self.AtomGroups[j]-1]+=self.Atoms[j][0];
        
        for i in range(0, groupcount):
            if not (mass[i] in self.Fragments):
                self.Fragments.append(mass[i]);
                if self.ReturnAtomNumbers:
                    self.GroupAtomNumbers.append([]);
                    FragmentIndex=len(self.Fragments)-1;
            else:
                if self.ReturnAtomNumbers:
                    FragmentIndex=self.Fragments.index(mass[i]);
            
            if self.ReturnAtomNumbers:
                    numbers=[];
                    for j in range(0, self.AtomCount):
                        if self.AtomGroups[j]==i:
                            numbers.append(j);
                    if self.ReturnAtomNumbers and (not (numbers in self.GroupAtomNumbers[FragmentIndex])):
                        self.GroupAtomNumbers[FragmentIndex].append(numbers);
                    
        
    def __ProcessFragmentationLevel(self,fragmentation_depth,current_level,startbond):
        for i in range(startbond,len(self.BreakableBonds)):
                    breakbond=self.BreakableBonds[i];
                    self.Atoms[breakbond[0]][2][breakbond[1]][1]=False;
                    self.Atoms[breakbond[2]][2][breakbond[3]][1]=False;
                    
                    #self.BondsArray[bond[0]][bond[1]]=0;
                    #self.BondsArray[bond[1]][bond[0]]=0;
        
                    self.__GenerateGroups(current_level);
                    
                    if current_level<fragmentation_depth:
                        self.__ProcessFragmentationLevel(fragmentation_depth,current_level+1,i+1);
                    self.Atoms[breakbond[0]][2][breakbond[1]][1]=True;
                    self.Atoms[breakbond[2]][2][breakbond[3]][1]=True;
                    
                    #self.BondsArray[bond[0]][bond[1]]=1;
                    #self.BondsArray[bond[1]][bond[0]]=1;

    
    
    
    def SortResults(self,Ascending=True):
            if Ascending:
                self.Fragments=sorted(self.Fragments,key=lambda sortit: sortit[0]);
            else:
                self.Fragments=sorted(self.Fragments,key=lambda sortit: sortit[0], reverse=True);
            
        
        
        
    def __init__(self,atoms,bonds,fragmentation_depth=2,Group_Atom_Numbers=False):
        '''
        Preparation of the fragmentation inputs
        '''
        if len(atoms)==0: 
            raise NameError('Empty Atoms!');
        self.ReturnAtomNumbers=Group_Atom_Numbers;
        self.AtomCount=len(atoms);
        self.Atoms=atoms;
        self.AtomGroups=[0]*self.AtomCount;
        #self.BondsArray= [[0] * self.AtomCount for i in range(self.AtomCount)];
        self.Bonds=bonds;
        self.BreakableBonds=[];
        for i in self.Atoms:
            i.append([]);
        #self.BrokenBonds=[];
        #self.BondForces=[];
        #for bond in self.Bonds:
        #    if len(bond)>2:
        #        self.BondsArray[bond[0]][bond[1]]=bond[2];
        #        self.BondsArray[bond[1]][bond[0]]=bond[2];
        #    else:
        #        self.BondsArray[bond[0]][bond[1]]=1;
        #        self.BondsArray[bond[1]][bond[0]]=1;
        #    if len(bond)>3:
        #        self.BondForces.append(bond[3]);
        #    else:
        #        self.BondForces.append(0.0);
        #self.BreakageThreshold=breakage_threshold;
        for i in range(0,len(self.Bonds)):
            bond=self.Bonds[i];
            self.Atoms[bond[0]][2].append([bond[1],True]);
            self.Atoms[bond[1]][2].append([bond[0],True]);
            if (self.Atoms[bond[0]][0]>1.1) and (self.Atoms[bond[1]][0]>1.1):
                self.BreakableBonds.append([bond[0],len(self.Atoms[bond[0]][2])-1,bond[1],len(self.Atoms[bond[1]][2])-1]);
        #print(len(self.BondForces),len(self.BreakableBonds));
        fragmentation_depth=min(fragmentation_depth,len(self.BreakableBonds));
        self.Fragments=[];
        if self.ReturnAtomNumbers:
            self.GroupAtomNumbers=[];
        #self.__ProcessFragment(0);
        self.__GenerateGroups(0);
        if fragmentation_depth>0:
            self.__ProcessFragmentationLevel(fragmentation_depth,1,0);
        for i in range(0,len(self.Fragments)):
            self.Fragments[i]=[self.Fragments[i]];
        if self.ReturnAtomNumbers:
            for i in range(0,len(self.Fragments)):
                self.Fragments[i].append(self.GroupAtomNumbers[i]);
            

def ProcessFragments(idx,atomlist,bondlist,fout):
    atoms=[];
    fout.write('%s\t'%idx);
    for s in atomlist:
        s=s.split(',');
        #s[0]=float(s[0]);
        #s[1]=int(s[3]);
        #s[3]=int(s[3]);
        #s[4]=int(s[4]);
        atoms.append([float(s[0]),int(s[3])]);
    bonds=[];
    for s in bondlist:
        s=s.split(',');
        #s[2]=float(s[3]);
        #s[0]=;
        #s[1]=int(s[1]);
        #s.pop();s.pop();
        #s[2]=int(s[2]);
        bonds.append([int(s[0]),int(s[1])]);
    if 'Neutral' in idx:
        frg=Fragmenter(atoms,bonds,2,False);
    else:
        frg=Fragmenter(atoms,bonds,2,True);
    frg.SortResults(False);
    fout.write('%s'%frg.Fragments[0][0]);
    for i in range(1,len(frg.Fragments)):
        fout.write(',%s'%frg.Fragments[i][0]);
    if not ('Neutral' in idx):
        fout.write('\t');
        fout.write('%s'%frg.GetCharges(frg.Fragments[0][1]));
        for i in range(1,len(frg.Fragments)):
            fout.write(',%s'%frg.GetCharges(frg.Fragments[i][1]));
    fout.write('\n');
        
        
        


def DoDefaultProcessing():

    inputfile='fraginput.txt';
    outputfile='fragments.frg';
    
    finp=open(inputfile,'r');
    inlist=[];
    for s in finp:
        s.replace('\n','');
        inlist.append(s);
    finp.close();
    
    
    i=0;
    fout=open(outputfile,'w');
    while i<len(inlist):
        idx=inlist[i];
        print(idx);
        i+=1;
        atomlist=inlist[i].split('\t');
        i+=1;
        if not (('Negative' in inlist[i]) or  ('Positive' in inlist[i]) or ('Neutral' in inlist[i])):
            bondlist=inlist[i].split('\t');
            ProcessFragments(idx,atomlist,bondlist,fout);
            i+=1;
        else:
            fout.write('%s\t'%idx);
            fout.write('%s'%atomlist[0].split(',')[0]);
            for j in range(1,len(atomlist)):
                fout.write(',%s'%atomlist[j].split(',')[0]);
            if  ('Negative' in inlist[i]) or  ('Positive' in inlist[i]):
                fout.write('\t');
                fout.write('%s'%atomlist[0].split(',')[3]);
                for j in range(1,len(atomlist)):
                    fout.write(',%s'%atomlist[j].split(',')[3]);
            fout.write('\n');
                
    fout.close();

if __name__=='__main__':
    #DoDefaultProcessing();
    #8_8_0_5_15688
    #atoms=[[12.0,'C',0,0,0],[12.0,'C',0,0,0],[12.0,'C',0,0,0],[31.972071,'S',0,0,0],[12.0,'C',0,0,0],[1.007825032,'H',0,0,0],[1.007825032,'H',0,0,0],[1.007825032,'H',0,0,0],[1.007825032,'H',0,0,0],[1.007825032,'H',0,0,0],[1.007825032,'H',0,0,0],[1.007825032,'H',0,0,0],[1.007825032,'H',0,0,0]];
    #bonds=[[0,1,1,1.0],[1,2,2,1.0],[2,3,1,1.0],[3,4,1,1.0],[2,5,1,0.0],[1,6,1,0.0],[0,7,1,0.0],[0,8,1,0.0],[0,9,1,0.0],[4,10,1,0.0],[4,11,1,0.0],[4,12,1,0.0]];
    
    atoms=[[12.0,0],[12.0,0],[15.99491462,0],[12.0,0],[12.0,0],[15.99491462,0],[1.007825032,0],[1.007825032,0],[1.007825032,0],[1.007825032,0],[1.007825032,0],[1.007825032,0],[1.007825032,0],[1.007825032,0]];
    bonds=[[0,1],[1,2],[1,3],[3,4],[3,5],[0,6],[0,7],[0,8],[1,9],[2,10],[4,11],[4,12],[4,13]];

    ff=Fragmenter(atoms,bonds,2,False);
    ff.SortResults(False);
    print(ff.Fragments);


    #'88.052429496,87.044604464,73.0289544,72.057514876,72.021129368,71.049689844,58.005479304,57.03403978,56.026214748,55.054775224,45.03403978,44.026214748,43.018389716,30.010564684,28.031300128,27.99491462,27.023475096,17.002739652,15.99491462,15.023475096,1.007825032'
    