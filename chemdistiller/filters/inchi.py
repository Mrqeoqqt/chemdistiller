# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 15:03:03 2016

@author: Dr. Ivan Laponogov
"""

from chemdistiller.settings import test_as_single_process, registered_scorers, registered_filters;
#match_type - 0 - soft, 1 - subclass of ref, 2 - subclass of test, 3 - strict masked, 4 - strict 

class InChIFilter:
    name='InChIFilter';
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
                if s[0].lower().startswith('inchi'):
                    self.ref_inchi=s[1].split('/');
                elif s[0].lower().startswith('match_type'):
                    self.match_type=int(s[1]);
                
            elif s.lower().startswith('use_short_inchi'):
                self.use_short_inchi=True;
                self.required_fields=['ShortInChI'];

            elif s.lower().startswith('end'):
                return;

    def _pipe_to_textfile(self, fout, indent=''):                
        fout.write('%s%s=%s\n'%(indent,'inchi',"/".join(self.ref_inchi)));
        fout.write('%s%s=%s\n'%(indent,'match_type',self.match_type));
        if self.use_short_inchi:
            fout.write('%suse_short_inchi\n'%indent);
        else:
            fout.write('%suse_full_inchi\n'%indent);

    def __init__(self, ref_inchi='', use_short_inchi=False, match_type=0):
        
        self.ref_inchi=ref_inchi.split('/');
        if '1S'==self.ref_inchi[0]:
            del self.ref_inchi[0];
        
        self.use_short_inchi=use_short_inchi;
        self.match_type=match_type;
        if use_short_inchi:
            self.required_fields=['ShortInChI'];
        else:
            self.required_fields=['InChI'];            
    
    def rejected(self, testrecord):
        if self.use_short_inchi:
            testinchi=testrecord['ShortInChI'];
        else:
            testinchi=testrecord['InChI'];
        result=False;
        testinchi=testinchi.split('/');
        if self.match_type==0:
            ml=min(len(testinchi),len(self.ref_inchi));
            for i in range(ml):
                if testinchi[i]!=self.ref_inchi[i]:
                    s1=testinchi[i];
                    s2=self.ref_inchi[i];
                    if len(s1)!=len(s2):
                        return True;
                    else:
                        for j in range(len(s1)):
                            if s1[j]!=s2[j]:
                                if s1[j]!='?' and s2[j]!='?':
                                    return True;
        elif self.match_type==1:
            if len(testinchi)<len(self.ref_inchi):
                return True;
            ml=min(len(testinchi),len(self.ref_inchi));
            for i in range(ml):
                if testinchi[i]!=self.ref_inchi[i]:
                    s1=testinchi[i];
                    s2=self.ref_inchi[i];
                    if len(s1)!=len(s2):
                        return True;
                    else:
                        for j in range(len(s1)):
                            if s1[j]!=s2[j]:
                                if s2[j]!='?':
                                    return True;            
        elif self.match_type==2:
            if len(testinchi)>len(self.ref_inchi):
                return True;
            ml=min(len(testinchi),len(self.ref_inchi));
            for i in range(ml):
                if testinchi[i]!=self.ref_inchi[i]:
                    s1=testinchi[i];
                    s2=self.ref_inchi[i];
                    if len(s1)!=len(s2):
                        return True;
                    else:
                        for j in range(len(s1)):
                            if s1[j]!=s2[j]:
                                if s1[j]!='?':
                                    return True;            
        elif self.match_type==3:
            if len(testinchi)!=len(self.ref_inchi):
                return True;
            else:
                for i in range(len(testinchi)):
                    if testinchi[i]!=self.ref_inchi[i]:
                        s1=testinchi[i];
                        s2=self.ref_inchi[i];
                        if len(s1)!=len(s2):
                            return True;
                        else:
                            for j in range(len(s1)):
                                if s1[j]!=s2[j]:
                                    if s1[j]!='?' and s2[j]!='?':
                                        return True;
                    
        else:
            if len(testinchi)!=len(self.ref_inchi):
                return True;
            else:
                for i in range(len(testinchi)):
                    if testinchi[i]!=self.ref_inchi[i]:
                        return True;
                
        
        
        return result;
        
        
    def __repr__(self):
        s='/'.join(self.ref_inchi);
        if self.match_type==0:
            m='soft';
        elif self.match_type==1:            
            m='subclass of ref';
        elif self.match_type==2:
            m='subclass of test';
        elif self.match_type==3:
            m='strict masked';
        elif self.match_type==4:
            m='strict';
        return 'InChi_Filter(ref="%s", %s)'%(s,m);

registered_filters[InChIFilter.__name__]=InChIFilter;

if __name__=='__main__':
    print('Soft, InChI');
    test=InChIFilter('C2H5OH/c1-2');  
    print(test.rejected({'InChI':'C2H5O'}));
    print(test.rejected({'InChI':'C2H5OH'}));
    print(test.rejected({'InChI':'C2H5O?'}));
    print(test.rejected({'InChI':'C2H5O?K'}));
    print(test.rejected({'InChI':'C2H5OH/c1-?'}));
    print(test.rejected({'InChI':'C2H5OH/c1-3'}));
    
    print('Strict, InChI');
    test=InChIFilter('C2H5OH/c1-2', match_type=4);  
    print(test.rejected({'InChI':'C2H5O'}));
    print(test.rejected({'InChI':'C2H5OH'}));
    print(test.rejected({'InChI':'C2H5O?'}));
    print(test.rejected({'InChI':'C2H5O?K'}));
    print(test.rejected({'InChI':'C2H5OH/c1-?'}));
    print(test.rejected({'InChI':'C2H5OH/c1-3'}));
    print(test.rejected({'InChI':'C2H5OH/c1-2'}));    
    
    print('Strict Masked, InChI');
    test=InChIFilter('C2H5OH/c1-2', match_type=3);  
    print(test.rejected({'InChI':'C2H5O'}));
    print(test.rejected({'InChI':'C2H5OH'}));
    print(test.rejected({'InChI':'C2H5O?'}));
    print(test.rejected({'InChI':'C2H5O?K'}));
    print(test.rejected({'InChI':'C2H5OH/c1-?'}));
    print(test.rejected({'InChI':'C2H5OH/c1-3'}));
    print(test.rejected({'InChI':'C2H5OH/c1-2'}));    
    
    print('SubClass of ref, ShortInChI');
    test=InChIFilter('C2H5OH/c1-2', use_short_inchi=True, match_type=1);  
    print(test.rejected({'ShortInChI':'C2H5O'}));
    print(test.rejected({'ShortInChI':'C2H5OH'}));
    print(test.rejected({'ShortInChI':'C2H5O?'}));
    print(test.rejected({'ShortInChI':'C2H5O?K'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-?'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-3'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-2'}));        
    test=InChIFilter('C2H5OH/c1-?', use_short_inchi=True, match_type=1);  
    print(test.rejected({'ShortInChI':'C2H5O'}));
    print(test.rejected({'ShortInChI':'C2H5OH'}));
    print(test.rejected({'ShortInChI':'C2H5O?'}));
    print(test.rejected({'ShortInChI':'C2H5O?K'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-?'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-3'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-2'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-2/h1'}));
    
    print('SubClass of test, ShortInChI');
    test=InChIFilter('C2H5OH/c1-2', use_short_inchi=True, match_type=2);  
    print(test.rejected({'ShortInChI':'C2H5O'}));
    print(test.rejected({'ShortInChI':'C2H5OH'}));
    print(test.rejected({'ShortInChI':'C2H5O?'}));
    print(test.rejected({'ShortInChI':'C2H5O?K'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-?'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-3'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-2'}));        
    test=InChIFilter('C2H5OH/c1-?', use_short_inchi=True, match_type=2);  
    print(test.rejected({'ShortInChI':'C2H5O'}));
    print(test.rejected({'ShortInChI':'C2H5OH'}));
    print(test.rejected({'ShortInChI':'C2H5O?'}));
    print(test.rejected({'ShortInChI':'C2H5O?K'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-?'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-3'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-2'}));
    print(test.rejected({'ShortInChI':'C2H5OH/c1-2/h1'}));
    
    print(test)