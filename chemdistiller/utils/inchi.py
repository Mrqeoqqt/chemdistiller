# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:01:56 2016

@author: Dr. Ivan Laponogov
"""

import hashlib;
import numpy as np;

#initialize the look up tables as described in ikey_base26

__base_string='ABCDEFGHIJKLMNOPQRSTUVWXYZ';

__triplets=[];

for i in range(26):
    if __base_string[i]!='E':
        for j in range(26):
            for k in range(26):
                c='%s%s%s'%(__base_string[i],__base_string[j],__base_string[k]);
                if c<'TAA' or c>'TTV':
                    __triplets.append(c);
                
__duplets=[];

for i in range(26):
    for j in range(26):
        __duplets.append('%s%s'%(__base_string[i],__base_string[j]));
        
__checksum_weights=[1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25];

#print(__base_string)
#print(__triplets)

def get_inchi_sha256_dig(inchi):
    s = bytearray(hashlib.sha256(inchi.encode()).digest());
    val=np.array(s,dtype=np.uint8);
    return val;

def get_short_inchi_from_full_inchi(inchi):
    return parse_inchi(inchi)[0];

def parse_inchi(inchi):
    #print(inchi)
    if inchi.upper().startswith('INCHI='):
        inchi=inchi[6:];
    #print(inchi)
    inchi=inchi.split('/');
    shortinchi='';
    restinchi='';
    prot='';
    mainblock=True;
    if inchi[0]=='1S':
        i=1;
    else:
        i=0;
    #print(inchi)
    if len(inchi)>i:
        shortinchi+=inchi[i];
                
        for j in range(i+1,len(inchi)):
            ij=inchi[j];
            if len(ij)>0:
                if (ij[0]=='c' or ij[0]=='h' or ij[0]=='q') and mainblock:
                    shortinchi+='/'+ij;
                elif ij[0]=='p':
                    prot=ij;
                    mainblock=False;
                else:
                    restinchi+='/'+ij;
                    mainblock=False;
            else:
                shortinchi+='/';
    else:
        shortinchi+='/';
        
    return [shortinchi,restinchi,prot];

    
def inchikey_from_inchi(inchi):
    global __base_string;
    shortinchi,restinchi,prot=parse_inchi(inchi);

    if prot!='':
        charge=int(prot.lstrip('p'));
        if abs(charge)>12:
            pstring='A';
        else:
            pstring=__base_string[charge+13];
    else:
        pstring='N';
    
    value=get_inchi_sha256_dig(shortinchi);
    
    t1=((value[0] & 0b11111111))    | ((value[1] & 0b00111111)<<8);
    t2=((value[1] & 0b11000000)>>6) | ((value[2] & 0b11111111)<<2) | ((value[3] & 0b00001111)<<10);
    t3=((value[3] & 0b11110000)>>4) | ((value[4] & 0b11111111)<<4) | ((value[5] & 0b00000011)<<12);
    t4=((value[5] & 0b11111100)>>2) | ((value[6] & 0b11111111)<<6);
    d1=((value[7] & 0b11111111))    | ((value[8] & 0b00000001)<<8);
    
    if (len(restinchi)>0) and (len(restinchi)<255):  
        value=get_inchi_sha256_dig('%s%s'%(restinchi, restinchi));
    else:
        value=get_inchi_sha256_dig(restinchi);
    
    r1=((value[0] & 0b11111111))    | ((value[1] & 0b00111111)<<8);
    r2=((value[1] & 0b11000000)>>6) | ((value[2] & 0b11111111)<<2) | ((value[3] & 0b00001111)<<10);
    rd1=((value[3] & 0b11110000)>>4)    | ((value[4] & 0b00011111)<<4);
    
    inchikey='%s-%sSA-%s'%(\
        ('%s%s%s%s%s'%(__triplets[t1],__triplets[t2],__triplets[t3],__triplets[t4],__duplets[d1])),\
        ('%s%s%s'%(__triplets[r1],__triplets[r2],__duplets[rd1])),\
        pstring);
    return inchikey;


def inchikeyvalues_from_inchi(inchi):
    global __base_string;
    shortinchi,restinchi,prot=parse_inchi(inchi);

    if prot!='':
        charge=int(prot.lstrip('p'));
    else:
        charge=0;
    
    value=get_inchi_sha256_dig(shortinchi);

    result=[];
    for i in range(8):
        result.append(value[i]);
    result.append(value[8] & 0b00000001);
    
    value=get_inchi_sha256_dig('%s%s'%(restinchi, restinchi));
    
    for i in range(4):
        result.append(value[i]);
    result.append(value[4] & 0b00011111);
    result.append(charge);
    
    return result;
    
def inchikey_from_inchikeyvalues(values):    
    if len(values)!=15:
        raise ValueError('Input InChIKeyValues array should be 15 elements long!');
    
    charge=values[-1];
    
    if abs(charge)>12:
        pstring='A';
    else:
        pstring=__base_string[charge+13];


    t1=((values[0] & 0b11111111))    | ((values[1] & 0b00111111)<<8);
    t2=((values[1] & 0b11000000)>>6) | ((values[2] & 0b11111111)<<2) | ((values[3] & 0b00001111)<<10);
    t3=((values[3] & 0b11110000)>>4) | ((values[4] & 0b11111111)<<4) | ((values[5] & 0b00000011)<<12);
    t4=((values[5] & 0b11111100)>>2) | ((values[6] & 0b11111111)<<6);
    d1=((values[7] & 0b11111111))    | ((values[8] & 0b00000001)<<8);

    r1=((values[9] & 0b11111111))    | ((values[10] & 0b00111111)<<8);
    r2=((values[10] & 0b11000000)>>6) | ((values[11] & 0b11111111)<<2) | ((values[12] & 0b00001111)<<10);
    rd1=((values[12] & 0b11110000)>>4)    | ((values[13] & 0b00011111)<<4);


    inchikey='%s-%sSA-%s'%(\
        ('%s%s%s%s%s'%(__triplets[t1],__triplets[t2],__triplets[t3],__triplets[t4],__duplets[d1])),\
        ('%s%s%s'%(__triplets[r1],__triplets[r2],__duplets[rd1])),\
        pstring);
    return inchikey;

def inchikeyvalues_from_inchikey(inchikey):
    
    inchikey=inchikey.rstrip('InChIKey=');
    if len(inchikey)!=27:
        raise ValueError('Input InChIKey should be 26 elements long!');
    
    major, minor, charge=inchikey.split('-');
    
    if minor[-2:]!='SA':
        raise ValueError('Unsupported InChiKey format');
        
    if charge=='A':
        charge=13;
    else:
        charge=__base_string.index(charge)-13;

    t1=__triplets.index(major[0:3]);
    t2=__triplets.index(major[3:6]);
    t3=__triplets.index(major[6:9]);
    t4=__triplets.index(major[9:12]);
    d1=__duplets.index(major[12:14]);
    
    r1=__triplets.index(minor[0:3]);
    r2=__triplets.index(minor[3:6]);
    rd1=__duplets.index(minor[6:8]);
    
    result=[t1 & 0b11111111,\
            ((t1>>8) & 0b00111111) | ((t2<<6) & 0b11000000),\
            ((t2>>2) & 0b11111111),\
            ((t2>>10) & 0b00001111) | ((t3<<4) & 0b11110000),\
            ((t3>>4) & 0b11111111),\
            ((t3>>12) & 0b00000011) | ((t4<<2) & 0b11111100),\
            ((t4>>6) & 0b11111111),\
            d1 & 0b11111111,\
            ((d1>>8) & 0b00000001),\
            r1 & 0b11111111,\
            ((r1>>8) & 0b00111111) | ((r2<<6) & 0b11000000),\
            ((r2>>2) & 0b11111111),\
            ((r2>>10) & 0b00001111) | ((rd1<<4) & 0b11110000),\
            ((rd1>>4) & 0b00011111),\
            charge];
    
    return result;
    
    
if __name__=='__main__':
    
    inchi='InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+/p+3';
    inchikey='IAQRGUVFOMOMEM-ONEGZZNKSA-N';
    print(inchi);
    print(inchikey_from_inchi(inchi));
    print(inchikeyvalues_from_inchi(inchi));
    print(inchikey_from_inchikeyvalues(inchikeyvalues_from_inchi(inchi)));
    print(inchikey);
    print(inchikeyvalues_from_inchikey(inchikey));
    print(parse_inchi(inchi))
    print(len(inchikeyvalues_from_inchikey(inchikey)))
    '''

    inchi='ClH/h1H';
    #inchikey='IAQRGUVFOMOMEM-ONEGZZNKSA-N';
    print(inchi);
    #print(inchikey_from_inchi(inchi));
    #print(inchikeyvalues_from_inchi(inchi));
    #print(inchikey_from_inchikeyvalues(inchikeyvalues_from_inchi(inchi)));
    #print(inchikey);
    #print(inchikeyvalues_from_inchikey(inchikey));
    print(parse_inchi(inchi))
    #print(len(inchikeyvalues_from_inchikey(inchikey)))
    '''
  
        
        
    


