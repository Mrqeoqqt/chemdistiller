# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 18:28:46 2016

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


from chemdistiller.msspectra.peak import MSPeak, MSBasicPeak;
from chemdistiller.annotations.peakannotation import MSPeakAnnotation;
from chemdistiller.chemdb.results import MolecularRecord;


class MSSpectrum:
    def __init__(self):
        self.peaks=[];
        self.parameters={};
        self.parent_peak=None;

    def close(self):
        for peak in self.peaks:
            if isinstance(peak, MSPeak):
                peak.close();
        
    def __enter__(self):
        return self;

    def __exit__(self, exc_type, exc_value, traceback):
        self.close();
    
    def add_peak(self, peak):
        self.peaks.append(peak);
        if isinstance(peak, MSPeak):
            peak.parent_spectrum=self;
    
    def load_from_textfile(self, fname):
        if os.path.isfile(fname):
            if self.peaks:
                self.close();
            with open(fname,'r') as finp:
                self._pipe_from_textfile(finp);
        else:
            raise IOError('File %s not found!'%fname);
            
    def _pipe_from_textfile(self, finp):
        current_peak=None;
        current_peak_number=0;
        while True:
            s=finp.readline();
            if s=='':
                 return;
            s=s.rstrip('\n').lstrip();
            if '##' in s:
                s=s[:s.index('##')];
                            
            if '=' in s:
                s=s.split('=',1);
                if s[0].lower().startswith('mode'):
                    self.parameters['mode']=int(s[1]);
                elif s[0].lower().startswith('collision_energy'):
                    self.parameters['collision_energy']=float(s[1]);
                elif s[0].lower().startswith('level'):
                    self.parameters['level']=int(s[1]);
                else:
                    self.parameters[s[0].lower()]=s[1];
            elif s.lower().startswith('peaks'):
                finished=False;
                while not finished:
                    s=finp.readline();
                    if s=='':
                         return;
                    s=s.rstrip('\n').lstrip();
                    if '##' in s:
                        s=s[:s.index('##')];
                    if '=' in s:
                        s=s.split('=', 1);
                        if not hasattr(current_peak,'parameters'):
                            if isinstance(current_peak, MSBasicPeak):
                                current_peak=MSPeak(current_peak);
                                current_peak.number=current_peak_number;
                                current_peak.parent_spectrum=self;
                                self.peaks[-1]=current_peak;
                            current_peak.parameters={};
                        current_peak.parameters[s[0].lower()]=s[1];
                    elif s.lower().startswith('spectrum'):
                        sub_spectrum=MSSpectrum();
                        sub_spectrum._pipe_from_textfile(finp);
                        if hasattr(current_peak,'ms_spectra'):
                            current_peak.ms_spectra.append(sub_spectrum);
                        else:
                            if isinstance(current_peak, MSBasicPeak):
                                current_peak=MSPeak(current_peak);
                                current_peak.number=current_peak_number;
                                current_peak.parent_spectrum=self;
                                self.peaks[-1]=current_peak;

                            current_peak.ms_spectra=[sub_spectrum];
                        sub_spectrum.parent_peak=current_peak;
                    elif s.lower().startswith('annotation'):
                        annotation=MSPeakAnnotation(current_peak);
                        annotation._pipe_from_textfile(finp);
                        if hasattr(current_peak,'annotations'):
                            current_peak.annotations.append(annotation);
                        else:
                            if isinstance(current_peak, MSBasicPeak):
                                current_peak=MSPeak(current_peak);
                                current_peak.number=current_peak_number;
                                current_peak.parent_spectrum=self;
                                self.peaks[-1]=current_peak;

                            current_peak.annotations=[annotation];
                    elif s.lower().startswith('merged_annotation'):
                        merged_annotation=MolecularRecord();
                        merged_annotation._pipe_from_textfile(finp);
                        
                        if hasattr(current_peak,'merged_annotations'):
                            current_peak.merged_annotations.append(merged_annotation);
                        else:
                            if isinstance(current_peak, MSBasicPeak):
                                current_peak=MSPeak(current_peak);
                                current_peak.number=current_peak_number;
                                current_peak.parent_spectrum=self;
                                self.peaks[-1]=current_peak;
                            current_peak.merged_annotations=[merged_annotation];
                        
                    elif (s.lower().startswith('end')):
                        finished=True;
                    else:
                        s=s.split(',');
                        if self.parameters['level']>1:
                            current_peak=MSBasicPeak();
                        else:
                            current_peak=MSPeak();
                            current_peak.number=int(s[0]);
                            current_peak.parent_spectrum=self;

                        current_peak.mz=float(s[1]);
                        current_peak.intensity=float(s[2]);
                        current_peak_number=int(s[0]);
                        if len(s)>3:
                            if isinstance(current_peak, MSBasicPeak):
                                current_peak=MSPeak(current_peak);
                                current_peak.number=current_peak_number;
                                current_peak.parent_spectrum=self;
                            current_peak.ppm=float(s[3]);

                        self.peaks.append(current_peak);
            elif s.lower().startswith('end'):
                return;
                
    def _pipe_to_textfile(self, fout, indent=''):                
        subindent='%s\t'%indent;
        for key in sorted(self.parameters.keys()):
            fout.write('%s%s=%s\n'%(indent,key,self.parameters[key]));
        if self.peaks:
            fout.write('%speaks\n'%indent);
            for index in range(len(self.peaks)):
                peak=self.peaks[index];
                if hasattr(peak,'ppm'):
                    if not (peak.ppm is None):
                        fout.write('%s%s,%s,%s,%s\n'%(subindent,peak.number,peak.mz,peak.intensity,peak.ppm));
                    else:
                        fout.write('%s%s,%s,%s\n'%(subindent,peak.number,peak.mz,peak.intensity));
                else:
                    if isinstance(peak, MSBasicPeak):
                        fout.write('%s%s,%s,%s\n'%(subindent,index,peak.mz,peak.intensity));
                    else:
                        fout.write('%s%s,%s,%s\n'%(subindent,peak.number,peak.mz,peak.intensity));
                
                if hasattr(peak,'parameters'):
                    for key in peak.parameters:
                        fout.write('%s\t%s=%s\n'%(subindent,key,peak.parameters[key]));

                if hasattr(peak,'annotations'):
                    if peak.annotations:
                        for annotation in peak.annotations:
                            fout.write('%s\tannotation\n'%subindent);
                            annotation._pipe_to_textfile(fout, '%s\t\t'%subindent);
                            fout.write('%s\tend\n'%subindent);
                            
                if hasattr(peak,'merged_annotations'):
                    for merged_annotation in peak.merged_annotations:
                        fout.write('%s\tmerged_annotation\n'%subindent);
                        merged_annotation._pipe_to_textfile(fout, '%s\t\t'%subindent);
                        fout.write('%s\tend\n'%subindent);

                if hasattr(peak,'ms_spectra'):                
                    if peak.ms_spectra:
                        for spectrum in peak.ms_spectra:
                            fout.write('%s\tspectrum\n'%subindent);
                            spectrum._pipe_to_textfile(fout, '%s\t\t'%subindent);
                            fout.write('%s\tend\n'%subindent);
            fout.write('%send\n'%indent);
            
    def save_to_textfile(self, fname):
        dirname,filename=os.path.split(os.path.abspath(fname));
        if not os.path.exists(dirname):
            os.makedirs(dirname);
        with open(fname,'w') as fout:
            self._pipe_to_textfile(fout);
        
    def normalize_to_one(self):
        if 'normalization' in self.parameters:
            if self.parameters['normalization']=='sum_one':
                return;
        total=0.0;
        for peak in self.peaks:
            total+=peak.intensity;
        if total>1e-15:    
            for peak in self.peaks:
                peak.intensity=peak.intensity/total;
        self.parameters['normalization']='sum_one';
    
    
