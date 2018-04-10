# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 16:11:40 2016

@author: ilaponog
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

class ProgressDisplay:
    def __init__(self, message, nproc, percent=True):
        self.message=message;
        self.nproc=nproc;
        self.percent=percent;
        self.counts=[0]*nproc;
        self.positions=[0]*nproc;
        print(message);
        if self.percent==False:
            for i in range(nproc):
                print_in_the_same_line('%5i          '%i);
            print('');
            print('----------------'*nproc);
        else:
            for i in range(nproc):
                print_in_the_same_line('%5i  '%i);
            print('');
            print('--------'*nproc);
            
        self.display_progress_from_queue((0,0,0));
        
        
    def display_progress_from_queue(self, progress_update):
        self.positions[progress_update[0]]=progress_update[1];
        self.counts[progress_update[0]]=progress_update[2];
        print_in_the_same_line('\r|');
        if self.percent:
            for i in range(self.nproc):
                if self.counts[i]>0:
                    p=float(self.positions[i])*100/self.counts[i];
                else:
                    p=0.0;
                print_in_the_same_line('%6.2f|'%(p));
        else:
            for i in range(self.nproc):
                print_in_the_same_line('%5i of %5i|'%(self.positions[i],self.counts[i]));


if __name__=='__main__':          
    import time;
    progress=ProgressDisplay('Test',5);
    for i in range(11):
        for j in range(5):
            progress.display_progress_from_queue((j,i,10));
        time.sleep(1);
        