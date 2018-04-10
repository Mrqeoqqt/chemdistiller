# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 19:56:30 2016

@author: Dr. Ivan Laponogov
"""

import sys;
import os;

if __name__=='__main__':
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    chemdistiller_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."));
    sys.path.append(chemdistiller_path);


from chemdistiller.utils.sysutils import print_in_the_same_line;

class TextProgressBar:
        def __init__(self, title='Progress: ', minvalue=0, maxvalue=100, position=0):
            self.title=title;
            self.minvalue=minvalue;
            self.maxvalue=maxvalue;
            self.set_position(position);
            
        
        def set_position(self, position):
            self.position=position;
            if self.position>self.maxvalue:
                self.position=self.maxvalue;
            if self.position<self.minvalue:
                self.position=self.minvalue;
            percent=(self.position-self.minvalue)*100/(self.maxvalue-self.minvalue);
            bar='\r%s[%s%s] %7.2f%%'%(self.title, '='*int(percent/2), ' '*int((100-percent)/2),percent);
            print_in_the_same_line(bar);
            #print(bar),

#%%   
if __name__=='__main__':
    
    import time;     
    progressbar=TextProgressBar();
    time.sleep(1);
    progressbar.set_position(10);
    time.sleep(1);
    progressbar.set_position(20);
    time.sleep(1);
    progressbar.set_position(30);    
    time.sleep(1);
    progressbar.set_position(40);    
    time.sleep(1);
    progressbar.set_position(50);    
    time.sleep(1);
    progressbar.set_position(60);    
    time.sleep(1);
    progressbar.set_position(70);    
    time.sleep(1);
    progressbar.set_position(80);
    time.sleep(1);
    progressbar.set_position(90);
    time.sleep(1);
    progressbar.set_position(100);