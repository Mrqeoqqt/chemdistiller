# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 17:19:42 2017

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
from chemdistiller.utils.inchi import inchikey_from_inchi;
from chemdistiller.utils.periodictable import parse_formula, \
                                              formula_to_element_vector, \
                                              encode_formula_to_array, \
                                              decode_formula_from_array;


    
try:
    from rdkit.Chem import AllChem;
    from rdkit.Chem import Draw;
    rdkit_available=True;
except:
    rdkit_available=False;
    

def generate_candidate_image(fname, smiles):
    global rdkit_available;
    if rdkit_available:
        try:
            mol = AllChem.MolFromSmiles(smiles);
            AllChem.Compute2DCoords(mol);
            Draw.MolToFile(mol, fname, size=(400,400));
            return True;
        except:
            return False;
    else:
        return False;

def generate_spectrum_view(spectrum, min_x, max_x, max_y):
    width=100;#percent
    height=200;
    
    left=4.0;
    top=30;
    bottom=height-30;
    right=width-2.0;
    
    result='<svg width="%s%%" height="%s">\n'%(width, height);            
    
    result+='<line x1="%s%%" y1="%s" x2="%s%%" y2="%s" style="stroke:rgb(0,0,20);stroke-width:2" />\n'%(left, top, left, bottom);
    result+='<line x1="%s%%" y1="%s" x2="%s%%" y2="%s" style="stroke:rgb(0,0,20);stroke-width:2" />\n'%(left, bottom, right, bottom);
    result+='<text x="%s%%" y="%s" fill="black">%.3f</text>\n'%(left-3.5,top,max_y);
    result+='<text x="%s%%" y="%s" fill="black">%.3f</text>\n'%(left-3.5,bottom,0.0);
    #result+='<text x="%s%%" y="%s" fill="black">%.1f</text>\n'%(left,bottom+25,min_x);
    #result+='<text x="%s%%" y="%s" fill="black">%.1f</text>\n'%(right,bottom+25,max_x);
    result+='<line x1="%s%%" y1="%s" x2="%s%%" y2="%s" style="stroke:rgb(0,0,20);stroke-width:1" />\n'%(left-0.4,top,left,top);
    result+='<line x1="%s%%" y1="%s" x2="%s%%" y2="%s" style="stroke:rgb(0,0,20);stroke-width:1" />\n'%(left-0.4,bottom,left,bottom);
    #result+='<line x1="%s" y1="%s" x2="%s" y2="%s" style="stroke:rgb(0,0,20);stroke-width:1" />\n'%(left,bottom,left,bottom+5);
    #result+='<line x1="%s" y1="%s" x2="%s" y2="%s" style="stroke:rgb(0,0,20);stroke-width:1" />\n'%(right,bottom,right,bottom+5);
    
    #result+='<text x="%s%%" y="%s" fill="black" transform="rotate(90 %s%%,%s)">%s</text>\n'%(left-3.5,int((top+bottom)/2), left-3.5,int((top+bottom)/2), 'Intensity');    
    
    result+='<text x="%s%%" y="%s" fill="black">%s</text>\n'%(right+0.2,bottom,'Da');        
    
    cnt=int((max_x-min_x)/10);
    
    for i in range(1, cnt+1):
        x_pos=(float(i*10+int(min_x/10)*10)-min_x)/(max_x-min_x)*(right-left)+left;
        result+='<line x1="%s%%" y1="%s" x2="%s%%" y2="%s" style="stroke:rgb(0,0,20);stroke-width:1" />\n'%(x_pos,bottom,x_pos,bottom+5);
        result+='<text x="%s%%" y="%s" fill="black">%s</text>\n'%(x_pos-0.2,bottom+25,i*10+int(min_x/10)*10);
    
    for peak in spectrum.peaks:
        x_pos=(peak.mz-min_x)/(max_x-min_x)*(right-left)+left;
        y_pos=peak.intensity/max_y*(bottom-top);
        result+='<line x1="%s%%" y1="%s" x2="%s%%" y2="%s" style="stroke:rgb(0,0,200);stroke-width:3" />\n'%(x_pos,bottom,x_pos,bottom-y_pos);
    
    
    return result+'Sorry, your browser does not support inline SVG.</svg>';

def generate_spectrum_peak_subspectra(peakfolder, index, peak):
    with open(peakfolder+'/peak_%s.html'%index,'w') as mainHTML:    
        mainHTML.write('<!DOCTYPE html>\n');
        mainHTML.write('<html>\n');
        mainHTML.write('<head>\n');
        mainHTML.write('<title>MS2 Spectra</title>\n');
        mainHTML.write('</head>\n');
        mainHTML.write('\n');
        mainHTML.write('<body>\n');
        mainHTML.write('<p><b>MS2 Spectra for peak N %s (m/z %.5f)</b></p>\n'%(index+1, peak.mz));
        mainHTML.write('<table border=1 width="100%%">\n');

        min_x=1e8;
        max_x=-1e8;
        
        max_y=0.0;

        for i in range(len(peak.ms_spectra)):
            spectrum=peak.ms_spectra[i];
            for subpeak in spectrum.peaks:
                x=subpeak.mz;
                y=subpeak.intensity;
                if x>max_x:
                    max_x=x;
                if x<min_x:
                    min_x=x;
                if y>max_y:
                    max_y=y;
                
        min_x=min_x-20;
        max_x=max_x+20;
        if max_y==0.0:
            max_y=1.0;
        
        for i in range(len(peak.ms_spectra)):
            spectrum=peak.ms_spectra[i];
            energy='Unknown';
            
            if 'collision_energy' in spectrum.parameters:
                if spectrum.parameters['collision_energy']>-1.0:
                    energy='%s eV'%spectrum.parameters['collision_energy'];
                else:
                    if 'collision_record' in spectrum.parameters:
                        energy=spectrum.parameters['collision_record'];


            svg=generate_spectrum_view(spectrum, min_x, max_x, max_y);
            
            mainHTML.write('<tr><th width=50>%s</th><td align="center" valign="middle">%s</td></tr>\n'%(energy, svg));
            
        mainHTML.write('</table>\n');
        mainHTML.write('</body>\n');


def generate_candidate_html(mainHTML, candidate, outpath, index, subpath, merged=False):
    if 'Annotation' in candidate:
        if candidate['Annotation']=='Correct':
            mainHTML.write('<tr><td bgcolor="#ccffcc">\n');
        else:
            mainHTML.write('<tr><td bgcolor="#e0e0e0">\n');

    else:        
        mainHTML.write('<tr><td>\n');

    if 'SMILES' in candidate:
        smiles=candidate['SMILES'];
    
    has_image=False;    
    if smiles!='':
        fname=os.path.join(outpath,'candidate_%s.png'%index)
        if generate_candidate_image(fname, smiles):
            has_image=True;
    
    if has_image:
        mainHTML.write('<table><tr><td><img src="%s/candidate_%s.png"></td><td>\n'%(subpath, index));
        
    mainHTML.write('<b>Total Score:</b> %s\n'%candidate['TotalScore']);
    mainHTML.write('<br><b>Individual Scores:</b>\n');
    for score in candidate['Scores']:
        mainHTML.write('<br><b>%s:</b> %s\n'%(score, candidate['Scores'][score]));
        
    mainHTML.write('<br><b>Mass:</b> %s, <b>Charge:</b> %s'%(candidate['Mass'],candidate['Charge']));
    
    smiles='';
    inchi='';
    formula='';
    if 'Formula' in candidate:
        formula=candidate['Formula'].formula_to_string();
    elif 'FormulaVector' in candidate:
        formula=decode_formula_from_array(candidate['FormulaVector']).formula_to_string();
        
    if formula!='':
        mainHTML.write(', <b>Formula:</b> %s'%formula);
        
    if merged:
        adduct='';isotope='';
        if 'Adduct' in candidate:
            adduct=candidate['Adduct'];
        if 'Isotope' in candidate:
            isotope=candidate['Isotope'];
        mainHTML.write('<br><b>Adduct:</b> %s, <b>Isotope:</b> %s\n'%(adduct, isotope));    
            
        
    if 'SMILES' in candidate:
        mainHTML.write('<br><b>SMILES:</b> %s\n'%candidate['SMILES']);

    if 'InChI' in candidate:
        mainHTML.write('<br><b>InChI</b>=1S/%s\n'%candidate['InChI']);
        inchi='InChI=1S/%s'%candidate['InChI'];
        inchikey=inchikey_from_inchi(inchi);
        mainHTML.write('<br><b>InChIKey:</b> %s\n'%inchikey);
        mainHTML.write('<br><button onclick="CallPubChem(this.id)" id="%s">PubChem</button><br>\n'%(inchikey));
            
    
    if has_image:
        mainHTML.write('</td></tr></table>\n');
    
    
    mainHTML.write('</td></tr>\n');
    

                
def generate_spectrum_peak_annotations(peakfolder, index, peak):
    with open(peakfolder+'/peak_a_%s.html'%index,'w') as mainHTML:    
        mainHTML.write('<!DOCTYPE html>\n');
        mainHTML.write('<html>\n');
        mainHTML.write('<head>\n');
        mainHTML.write('<title>Peak Annotations</title>\n');
        mainHTML.write('</head>\n');
        mainHTML.write('\n');
        
        
        mainHTML.write('<script>\n');
        mainHTML.write('function CallPubChem(inchikey){\n');
            
        mainHTML.write('document.getElementById(inchikey).innerHTML="Loading...";\n');
        mainHTML.write('var s="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"+inchikey+"/cids/TXT";\n');
            
        mainHTML.write('var xmlHttp = new XMLHttpRequest();\n');
        mainHTML.write('xmlHttp.onreadystatechange = function() { \n');
        mainHTML.write('if (xmlHttp.readyState == 4 && xmlHttp.status == 200) {\n');
        mainHTML.write('var res = xmlHttp.responseText;\n');
            
        mainHTML.write('ss="https://pubchem.ncbi.nlm.nih.gov/compound/"+res;\n');
        mainHTML.write('document.getElementById(inchikey).innerHTML=ss;  \n');
        mainHTML.write('window.open(ss,"_blank");\n');

        mainHTML.write('};};\n');
        mainHTML.write('xmlHttp.open("GET", s, true); // true for asynchronous \n');
        mainHTML.write('xmlHttp.send(null);\n');
    
        mainHTML.write('}\n');
        mainHTML.write('</script>');
        
        mainHTML.write('<body>\n');
        mainHTML.write('<p><b>Separated per adduct annotations for peak N %s (m/z %.5f)</b></p>\n'%(index+1, peak.mz));
        mainHTML.write('<table border=1 width="100%">\n');
        ii=-1;

        for annotation in peak.annotations:
            ii+=1;
            if not os.path.exists(os.path.join(peakfolder,'peak_a_%s_%s'%(index, ii))):
                os.makedirs(os.path.join(peakfolder,'peak_a_%s_%s'%(index, ii)));
            outpath=os.path.join(peakfolder,'peak_a_%s_%s'%(index, ii));    
            mainHTML.write('<tr><td>\n');
            if not (annotation.adduct is None):
                mainHTML.write('<b>Adduct:</b> %s\n'%annotation.adduct.definition);
            mainHTML.write('<br><b>Isotope:</b> %s\n'%annotation.isotope);
            
            if not (annotation.formula_scorer is None):
                mainHTML.write('<br><b>Formula Scorer:</b> %s\n'%annotation.formula_scorer);
            if not (annotation.element_scorer is None):    
                mainHTML.write('<br><b>Element Scorer:</b> %s\n'%annotation.element_scorer);
            
            if annotation.filters:
                for i in annotation.filters:
                    mainHTML.write('<br><b>Filter:</b> %s\n'%i);
                
            if not (annotation.scores is None):
                for score in annotation.scores:
                    mainHTML.write('<br><b>Score %s:</b> %s\n'%(score, annotation.scores[score]));
        
            for parameter in annotation.parameters:
                mainHTML.write('<br><b>Parameter %s:</b> %s\n'%(parameter, annotation.parameters[parameter]));
                
            if not (annotation.mol_candidates is None):
                mainHTML.write('<br><b>Candidate molecules:</b><br>\n');
                mainHTML.write('<table border=1 width="100%">\n');
                jj=-1;
                for candidate in annotation.mol_candidates.mol_list:
                    jj+=1;
                    generate_candidate_html(mainHTML, candidate, outpath, jj, 'peak_a_%s_%s'%(index, ii), merged=False);
                
                mainHTML.write('</table>\n');

                
            mainHTML.write('</td></tr>\n');
        
        #peak annotations here

        mainHTML.write('</table>\n');
    
        mainHTML.write('</body>\n');

    
    
         
def generate_spectrum_peak_merged_annotations(peakfolder, index, peak):
    with open(peakfolder+'/peak_ma_%s.html'%index,'w') as mainHTML:    
        mainHTML.write('<!DOCTYPE html>\n');
        mainHTML.write('<html>\n');
        mainHTML.write('<head>\n');
        mainHTML.write('<title>MS2 Spectra</title>\n');
        mainHTML.write('</head>\n');
        mainHTML.write('\n');
        mainHTML.write('<body>\n');
        mainHTML.write('<p><b>Merged annotations for peak N %s (m/z %.5f)</b></p>\n'%(index+1, peak.mz));
        
        if not (peak.merged_annotations is None):
            mainHTML.write('<table border=1 width="100%">\n');

            if not os.path.exists(os.path.join(peakfolder,'peak_ma_%s'%index)):
                os.makedirs(os.path.join(peakfolder,'peak_ma_%s'%index));
            outpath=os.path.join(peakfolder,'peak_ma_%s'%index);    
 
            jj=-1;
            for candidate in peak.merged_annotations:
                jj+=1;
                generate_candidate_html(mainHTML, candidate, outpath, jj, 'peak_ma_%s'%(index), merged=True);
                
            mainHTML.write('</table>\n');
    
        mainHTML.write('</body>\n');

    
    


def generate_index(outpath):
    with open(outpath+'/index.html','w') as frameHTML:
        frameHTML.write('<!DOCTYPE html>\n');
        frameHTML.write('<html>\n');
        frameHTML.write('<frameset cols="10%,90%">\n');
        frameHTML.write('<frame src="main.html">\n');
        frameHTML.write('<frame src="Spectra/spec_0.html" name="MS1Spectrum">\n');
        frameHTML.write('</frameset>\n');
        frameHTML.write('</html>\n');
        

def generate_main(outpath, specmanager):
    with open(outpath+'/main.html','w') as mainHTML:
        mainHTML.write('<!DOCTYPE html>\n');
        mainHTML.write('<html>\n');
        mainHTML.write('<head>\n');
        mainHTML.write('<title>MS1 Spectra</title>\n');
        mainHTML.write('</head>\n');
        mainHTML.write('\n');
        mainHTML.write('<body>\n');
        mainHTML.write('<table border=1>\n');
        mainHTML.write('<tr><th>MS1 Spectra</th></tr>\n');

        for index in range(len(specmanager.ms_spectra)):
            spectrum=specmanager.ms_spectra[index];
            if 'title' in spectrum.parameters:
                title=spectrum.parameters['title'];
            elif 'db_molecule_name' in spectrum.parameters:
                title=spectrum.parameters['db_molecule_name'];
            elif 'smiles' in spectrum.parameters:
                title=spectrum.parameters['smiles'];
            elif 'filename' in spectrum.parameters:
                title=spectrum.parameters['filename'];
            else:
                title='%s'%index;
            spectrum._title=title;
            mainHTML.write('<tr><td>%s. <a href="Spectra/spec_%s.html" target="MS1Spectrum">%s</a></td></tr>\n'%(index+1, index, title));
        
        mainHTML.write('</table>\n');
        mainHTML.write('</body>\n');


    
def generate_spectrum(outpath, index, spectrum):
    with open(outpath+'/spec_%s.html'%index,'w') as frameHTML:
        frameHTML.write('<!DOCTYPE html>\n');
        frameHTML.write('<html>\n');
        frameHTML.write('<frameset rows="40%,60%">\n');
        frameHTML.write('<frame src="specdesc_%s.html">\n'%index);
        frameHTML.write('<frameset cols="20%,80%">\n');
        frameHTML.write('<frame src="specpeaks_%s.html">\n'%index);
        frameHTML.write('<frame src="SpecPeaks_%s/peak_0.html" name="SpectrumPeak">\n'%index);
        frameHTML.write('</frameset>\n');
        frameHTML.write('</frameset>\n');
        frameHTML.write('</html>\n');
        
    with open(outpath+'/specdesc_%s.html'%index,'w') as mainHTML:
        mainHTML.write('<!DOCTYPE html>\n');
        mainHTML.write('<html>\n');
        mainHTML.write('<head>\n');
        mainHTML.write('<title>%s Spectrum</title>\n'%spectrum._title);
        mainHTML.write('</head>\n');
        mainHTML.write('\n');
        mainHTML.write('<body>\n');
        mainHTML.write('<p><b>%s spectrum</b></p>\n'%spectrum._title);
        mainHTML.write('<table border=1 width="100%">\n');
        
        
        min_x=1e8;
        max_x=-1e8;
        
        max_y=0.0;
        
        for peak in spectrum.peaks:
            x=peak.mz;
            y=peak.intensity;
            if x>max_x:
                max_x=x;
            if x<min_x:
                min_x=x;
            if y>max_y:
                max_y=y;
                
                
        min_x=min_x-20;
        max_x=max_x+20;
        if max_y==0.0:
            max_y=1.0;
        
        svg=generate_spectrum_view(spectrum, min_x, max_x, max_y);
        mainHTML.write('<tr><td colspan="2" align="center" valign="middle">%s</td></tr>\n'%(svg));
        l=[];
        
        if 'db_molecule_name' in spectrum.parameters:
            l.append(('Molecule name', spectrum.parameters['db_molecule_name']));

        if 'exactmass' in spectrum.parameters:
            l.append(('Exact Mass, Da', spectrum.parameters['exactmass']));

        if 'formula' in spectrum.parameters:
            l.append(('Formula', spectrum.parameters['formula']));

        if 'inchi' in spectrum.parameters:
            l.append(('InChI', spectrum.parameters['inchi']));
        
        if 'mode' in spectrum.parameters:
            if spectrum.parameters['mode']==1:
                l.append(('Mode','positive'));
            else:
                l.append(('Mode','negative'));

        if 'db_source' in spectrum.parameters:
            l.append(('Source database', spectrum.parameters['db_source']));

        if 'db_entry_id' in spectrum.parameters:
            l.append(('Entry IDs', spectrum.parameters['db_entry_id']));
        
        if 'db_authors' in spectrum.parameters:
            l.append(('Database authors', spectrum.parameters['db_authors']));
            
        if 'db_copyright' in spectrum.parameters:
            l.append(('Database copyright', spectrum.parameters['db_copyright']));
            
        if 'db_license' in spectrum.parameters:
            l.append(('Database license', spectrum.parameters['db_license']));

        for i in l:
            mainHTML.write('<tr><th>%s</th><td>%s</td></tr>\n'%i);
        
        mainHTML.write('</table>\n');
        mainHTML.write('</body>\n');

    with open(outpath+'/specpeaks_%s.html'%index,'w') as mainHTML:
        mainHTML.write('<!DOCTYPE html>\n');
        mainHTML.write('<html>\n');
        mainHTML.write('<head>\n');
        mainHTML.write('<title>%s Spectrum</title>\n'%spectrum._title);
        mainHTML.write('</head>\n');
        mainHTML.write('\n');
        mainHTML.write('<body>\n');
        mainHTML.write('<p><b>MS1 peaks</b></p>\n');
        mainHTML.write('<table border=1>\n');
        mainHTML.write('<tr><th>N</th><th>m/z</th><th>Int.</th><th>Ion</th><th>MS2</th><th>Annot.</th></tr>\n');
        for i in range(len(spectrum.peaks)):
            
            
            peak=spectrum.peaks[i];
            if 'ion_type' in peak.parameters:
                ion=peak.parameters['ion_type'];
            else:
                ion='';
            if hasattr(peak,'ms_spectra'):
                if peak.ms_spectra:
                    spec=len(peak.ms_spectra);
            else:
                spec='';
                
            if hasattr(peak,'annotations'):
                if peak.annotations:
                    ann1=len(peak.annotations);
            else:
                ann1='';

            if hasattr(peak,'merged_annotations'):
                if peak.merged_annotations:
                    ann2=len(peak.merged_annotations);
            else:
                ann2='';

            if ann1=='' and ann2=='':
                annot='';
            elif ann1=='' and ann2!='':
                annot='<a href="SpecPeaks_%s/peak_ma_%s.html" target="SpectrumPeak">Merged: %s</a>'%(index, i, ann2);
            elif ann1!='' and ann2=='':
                annot='<a href="SpecPeaks_%s/peak_a_%s.html" target="SpectrumPeak">Indiv.:%s</a>'%(index, i, ann1);
            else:
                annot='<a href="SpecPeaks_%s/peak_a_%s.html" target="SpectrumPeak">Indiv.:%s</a><br><a href="SpecPeaks_%s/peak_ma_%s.html" target="SpectrumPeak">Merged:%s</a>'%(index, i, ann1, index, i, ann2);
                
            mainHTML.write('<tr><td>%s</td><td>%.5f</td><td>%.5f</td><td>%s</td><td><a href="SpecPeaks_%s/peak_%s.html" target="SpectrumPeak">%s</a></td><td>%s</td></tr>\n'%(i+1, peak.mz, peak.intensity, ion, index, i, spec, annot));
            if spec!='' or ann1!='' or ann2!='':
                peakfolder=os.path.join(outpath, 'SpecPeaks_%s'%index);
                if not os.path.exists(peakfolder):
                    os.makedirs(peakfolder);
                if spec!='':
                    generate_spectrum_peak_subspectra(peakfolder, i, spectrum.peaks[i]);
                if ann1!='':
                    generate_spectrum_peak_annotations(peakfolder, i, spectrum.peaks[i]);
                if ann2!='':
                    generate_spectrum_peak_merged_annotations(peakfolder, i, spectrum.peaks[i]);
                    
        mainHTML.write('</table>\n');
        mainHTML.write('</body>\n');



def generate_HTML_report(outpath, specmanager):
    if not os.path.exists(outpath):
        os.makedirs(outpath);
    
    specfolder=os.path.join(outpath,'Spectra');
    if not os.path.exists(specfolder):
        os.makedirs(specfolder);
    
    generate_index(outpath);
    
    generate_main(outpath, specmanager);
    print('Generating reports:\n');
    
    for i in range(len(specmanager.ms_spectra)):
        print_in_the_same_line('\rSpectrum %s of %s'%(i+1,len(specmanager.ms_spectra)));
        spectrum=specmanager.ms_spectra[i];
        generate_spectrum(specfolder, i, spectrum);
        if hasattr(spectrum,'_title'):
            del spectrum._title;
    print('\nFinished');
    
if __name__=='__main__':
    
    spectral_input_path='e:/Imperial/TestDB/totest';
    outpath='e:/Imperial/TestDB/totestHTMLReport';
    from chemdistiller.msspectra.manager import SpectralManager;
    specmanager=SpectralManager();
    specmanager.import_textfile_spectra_from_folder(spectral_input_path);
    generate_HTML_report(outpath, specmanager);
    specmanager.close();
    
    
    
    
    
    