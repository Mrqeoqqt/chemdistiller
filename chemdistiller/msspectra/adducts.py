# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 15:06:32 2016

@author: Dr. Ivan Laponogov
"""

class Adduct:
    #mode, text definition, extra mass, multiplier, charge before adduct formation
    def __init__(self, adduct_definition=[1, '[M+H]+', 1.007276, 1.0, 0]):
        self.definition=adduct_definition[1];
        self.mode=adduct_definition[0];
        self.__substr=adduct_definition[2];
        self.__mult=adduct_definition[3];
        self.charge=adduct_definition[4];
        
        
    def get_mass(self, mz):
        return (mz-self.__substr)*self.__mult;

    def get_mz(self, mass):
        return mass/self.__mult+self.__substr;
        
        

# http://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/
#
#Huang N.; Siegel M.M.1; Kruppa G.H.; Laukien F.H.; J Am Soc Mass Spectrom 1999, 10, 1166–1173; 
#Automation of a Fourier transform ion cyclotron resonance mass spectrometer for acquisition, 
#analysis, and e-mailing of high-resolution exact-mass electrospray ionization mass spectral data.
#
#Interferences and contaminants encountered in modern mass spectrometry (Bernd O. Keller, Jie Sui, 
#Alex B. Young and Randy M. Whittal, ANALYTICA CHIMICA ACTA, 627 (1): 71-81)

    
positive_mode_adducts=[
Adduct([1,'[M]+', 0.0, 1.0, 1]),\
Adduct([1,'[M+3H]3+', 1.007276466879, 3.0, 0]),\
Adduct([1,'[M+2H+Na]3+',8.334590, 3.0, 0]),\
Adduct([1,'[M+H+2Na]3+',15.7661904, 3.0, 0]),\
Adduct([1,'[M+3Na]3+',22.989218, 3.0, 0]),\
Adduct([1,'[M+2H]2+',1.007276466879, 2.0, 0]),\
Adduct([1,'[M+H+NH4]2+',9.520550, 2.0, 0]),\
Adduct([1,'[M+H+Na]2+',11.998247, 2.0, 0]),\
Adduct([1,'[M+H+K]2+',19.985217, 2.0, 0]),\
Adduct([1,'[M+CH3CN+2H]2+',21.520550, 2.0, 0]),\
Adduct([1,'[M+2Na]2+',22.989218, 2.0, 0]),\
Adduct([1,'[M+2CH3CN+2H]2+',42.033823, 2.0, 0]),\
Adduct([1,'[M+3CH3CN+2H]2+',62.547097, 2.0, 0]),\
Adduct([1,'[M+H]+',1.007276466879, 1.0, 0]),\
Adduct([1,'[M+NH4]+',18.033823, 1.0, 0]),\
Adduct([1,'[M+Na]+',22.989218, 1.0, 0]),\
Adduct([1,'[M+CH3OH+H]+',33.033489, 1.0, 0]),\
Adduct([1,'[M+K]+',38.963158, 1.0, 0]),\
Adduct([1,'[M+CH3CN+H]+',42.033823, 1.0, 0]),\
Adduct([1,'[M+2Na-H]+',44.971160, 1.0, 0]),\
Adduct([1,'[M+C3H7OH+H]+',61.06534, 1.0, 0]),\
Adduct([1,'[M+CH3CN+Na]+',64.015765, 1.0, 0]),\
Adduct([1,'[M+2K-H]+',76.919040, 1.0, 0]),\
Adduct([1,'[M+C2H6OS+H]+',79.02122, 1.0, 0]),\
Adduct([1,'[M+2CH3CN+H]+',83.060370, 1.0, 0]),\
Adduct([1,'[M+C3H7OH+Na+H]+',84.05511, 1.0, 0]),\
Adduct([1,'[2M+H]+',1.007276466879, 0.5, 0]),\
Adduct([1,'[2M+NH4]+',18.033823, 0.5, 0]),\
Adduct([1,'[2M+Na]+',22.989218, 0.5, 0]),\
Adduct([1,'[2M+K]+',38.963158, 0.5, 0]),\
Adduct([1,'[2M+CH3CN+H]+',42.033823, 0.5, 0]),\
Adduct([1,'[2M+CH3CN+Na]+',64.015765, 0.5, 0]),\
];

negative_mode_adducts=[\
Adduct([-1,'[M]-', 0.0, 1.0, -1]),\
Adduct([-1,'[M-H]-', -1.007276466879, 1.0, 0]),\
Adduct([-1,'[M-3H]3-', -1.007276466879, 3.0, 0]),\
Adduct([-1,'[M-2H]2-', -1.007276466879, 2.0, 0]),\
Adduct([-1,'[M-H2O-H]-', -19.01839, 1.0, 0]),\
Adduct([-1,'[M+Na-2H]-', 20.974666, 1.0, 0]),\
Adduct([-1,'[M+Cl]-', 34.969402, 1.0, 0]),\
Adduct([-1,'[M+K-2H]-', 36.948606, 1.0, 0]),\
Adduct([-1,'[M+HCOOH-H]-', 44.998201, 1.0, 0]),\
Adduct([-1,'[M+CH3COOH-H]-', 59.013851, 1.0, 0]),\
Adduct([-1,'[M+Br]-', 78.918885, 1.0, 0]),\
Adduct([-1,'[M+C2HF3O2-H]-', 112.985586, 1.0, 0]),\
Adduct([-1,'[2M-H]-', -1.007276466879, 0.5, 0]),\
Adduct([-1,'[2M+HCOOH-H]-', 44.998201, 0.5, 0]),\
Adduct([-1,'[2M+CH3COOH-H]-', 59.013851, 0.5, 0]),\
Adduct([-1,'[3M-H]-', -1.007276466879, 1.0/3.0, 0]),\
];


def get_adduct_by_name(adduct_name, mode=0):
    result=None;
    lst=None;
    if mode==-1:
        lst=negative_mode_adducts;
    elif mode==1:
        lst=positive_mode_adducts;
    else:
        lst=negative_mode_adducts+positive_mode_adducts;
    if not (lst is None):
        for i in lst:
            if i.definition.upper()==adduct_name.upper():
                result=i;
                break;
        
    return result;
    

def get_positive_mode_adducts(adducts=[]):
    result=[];
    for adduct in positive_mode_adducts:
        if adduct.definition in adducts:
            result.append(adduct);
    return result;        


def get_negative_mode_adducts(adducts=[]):
    result=[];
    for adduct in negative_mode_adducts:
        if adduct.definition in adducts:
            result.append(adduct);
    return result;        


global_supported_adducts=[];

if __name__=='__main__':
    
    for i in range(len(positive_mode_adducts)):
        print(positive_mode_adducts[i].definition);
        
    for i in range(len(negative_mode_adducts)):
        print(negative_mode_adducts[i].definition);        
        
    print(get_adduct_by_name('[M+H]+',1).definition);
    