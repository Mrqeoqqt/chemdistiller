"""
Mendeleev's table related utilities


@author: Dr. Ivan Laponogov


"""

import numpy as np;
import copy;
from operator import itemgetter;
from math import exp, log;

#Isotopic data derived from: https://www.ncsu.edu/chemistry/msf/pdf/IsotopicMass_NaturalAbundance.pdf
#Where no isotopic info available, the natural abandance is set to 100%
#Refs. G. Audi, A. H. Wapstra Nucl. Phys A 1993,565, 1-65 and G. Audi, A. H. Wapstra Nucl. Phys A 1995,595, 409-480
#Refs. The 1997 report of the IUPAC Subcommittee for Isotopic Abundance Measurements by K.J.R. Rosman, P.D.P. Taylor
#      Pure Appl. Chem. 1999, 71, 1593-1607

#elements_list_full:
#Element No, Element Name, Element Full Name, Natural isotopes, Natural isotope masses (AU), Natural isotope abundance (%)

proton_mass=1.007276466879;
electron_mass=5.48579909070e-04;

hydrogen_mass=proton_mass+electron_mass;

elements_list_full=[\
[1, 'H', 'Hydrogen', [1, 2], [1.007825, 2.014102], [99.9885, 0.0115]],\
[2, 'He', 'Helium', [4, 3], [4.002603, 3.016029], [99.999863, 0.000137]],\
[3, 'Li', 'Lithium', [7, 6], [7.016004, 6.015122], [92.41, 7.59]],\
[4, 'Be', 'Beryllium', [9], [9.012182], [100.0]],\
[5, 'B', 'Boron', [11, 10], [11.009305, 10.012937], [80.1, 19.9]],\
[6, 'C', 'Carbon', [12, 13], [12.0, 13.003355], [98.93, 1.07]],\
[7, 'N', 'Nitrogen', [14, 15], [14.003074, 15.000109], [99.632, 0.368]],\
[8, 'O', 'Oxygen', [16, 18, 17], [15.994915, 17.99916, 16.999132], [99.757, 0.205, 0.038]],\
[9, 'F', 'Fluorine', [19], [18.998403], [100.0]],\
[10, 'Ne', 'Neon', [20, 22, 21], [19.99244, 21.991386, 20.993847], [90.48, 9.25, 0.27]],\
[11, 'Na', 'Sodium', [23], [22.98977], [100.0]],\
[12, 'Mg', 'Magnesium', [24, 26, 25], [23.985042, 25.982593, 24.985837], [78.99, 11.01, 10]],\
[13, 'Al', 'Aluminum', [27], [26.981538], [100.0]],\
[14, 'Si', 'Silicon', [28, 29, 30], [27.976927, 28.976495, 29.97377], [92.2297, 4.6832, 3.0872]],\
[15, 'P', 'Phosphorus', [31], [30.973762], [100.0]],\
[16, 'S', 'Sulfur', [32, 34, 33, 36], [31.972071, 33.967867, 32.971458, 35.967081], [94.93, 4.29, 0.76, 0.02]],\
[17, 'Cl', 'Chlorine', [35, 37], [34.968853, 36.965903], [75.78, 24.22]],\
[18, 'Ar', 'Argon', [40, 36, 38], [39.962383, 35.967546, 37.962732], [99.6003, 0.3365, 0.0632]],\
[19, 'K', 'Potassium', [39, 41, 40], [38.963707, 40.961826, 39.963999], [93.2581, 6.7302, 0.0117]],\
[20, 'Ca', 'Calcium', [40, 44, 42, 48, 43, 46], [39.962591, 43.955481, 41.958618, 47.952534, 42.958767, 45.953693], [96.941, 2.086, 0.647, 0.187, 0.135, 0.004]],\
[21, 'Sc', 'Scandium', [45], [44.95591], [100.0]],\
[22, 'Ti', 'Titanium', [48, 46, 47, 49, 50], [47.947947, 45.952629, 46.951764, 48.947871, 49.944792], [73.72, 8.25, 7.44, 5.41, 5.18]],\
[23, 'V', 'Vanadium', [51, 50], [50.943964, 49.947163], [99.75, 0.25]],\
[24, 'Cr', 'Chromium', [52, 53, 50, 54], [51.940512, 52.940654, 49.94605, 53.938885], [83.789, 9.501, 4.345, 2.365]],\
[25, 'Mn', 'Manganese', [55], [54.93805], [100.0]],\
[26, 'Fe', 'Iron', [56, 54, 57, 58], [55.934942, 53.939615, 56.935399, 57.93328], [91.754, 5.845, 2.119, 0.282]],\
[27, 'Co', 'Cobalt', [59], [58.9332], [100.0]],\
[28, 'Ni', 'Nickel', [58, 60, 62, 61, 64], [57.935348, 59.930791, 61.928349, 60.93106, 63.92797], [68.0769, 26.2231, 3.6345, 1.1399, 0.9256]],\
[29, 'Cu', 'Copper', [63, 65], [62.929601, 64.927794], [69.17, 30.83]],\
[30, 'Zn', 'Zinc', [64, 66, 68, 67, 70], [63.929147, 65.926037, 67.924848, 66.927131, 69.925325], [48.63, 27.9, 18.75, 4.1, 0.62]],\
[31, 'Ga', 'Gallium', [69, 71], [68.925581, 70.924705], [60.108, 39.892]],\
[32, 'Ge', 'Germanium', [74, 72, 70, 73, 76], [73.921178, 71.922076, 69.92425, 72.923459, 75.921403], [36.28, 27.54, 20.84, 7.73, 7.61]],\
[33, 'As', 'Arsenic', [75], [74.921596], [100.0]],\
[34, 'Se', 'Selenium', [80, 78, 76, 82, 77, 74], [79.916522, 77.91731, 75.919214, 81.9167, 76.919915, 73.922477], [49.61, 23.77, 9.37, 8.73, 7.63, 0.89]],\
[35, 'Br', 'Bromine', [79, 81], [78.918338, 80.916291], [50.69, 49.31]],\
[36, 'Kr', 'Krypton', [84, 86, 82, 83, 80, 78], [83.911507, 85.91061, 81.913485, 82.914136, 79.916378, 77.920386], [57, 17.3, 11.58, 11.49, 2.28, 0.35]],\
[37, 'Rb', 'Rubidium', [85, 87], [84.911789, 86.909183], [72.17, 27.83]],\
[38, 'Sr', 'Strontium', [88, 86, 87, 84], [87.905614, 85.909262, 86.908879, 83.913425], [82.58, 9.86, 7, 0.56]],\
[39, 'Y', 'Yttrium', [89], [88.905848], [100.0]],\
[40, 'Zr', 'Zirconium', [90, 94, 92, 91, 96], [89.904704, 93.906316, 91.90504, 90.905645, 95.908276], [51.45, 17.38, 17.15, 11.22, 2.8]],\
[41, 'Nb', 'Niobium', [93], [92.906378], [100.0]],\
[42, 'Mo', 'Molybdenum', [98, 96, 95, 92, 100, 97, 94], [97.905408, 95.904679, 94.905841, 91.90681, 99.907477, 96.906021, 93.905088], [24.13, 16.68, 15.92, 14.84, 9.63, 9.55, 9.25]],\
[43, 'Tc', 'Technetium', [98], [97.907216], [100.0]],\
[44, 'Ru', 'Ruthenium', [102, 104, 101, 99, 100, 96, 98], [101.90435, 103.90543, 100.905582, 98.905939, 99.90422, 95.907598, 97.905287], [31.55, 18.62, 17.06, 12.76, 12.6, 5.54, 1.87]],\
[45, 'Rh', 'Rhodium', [103], [102.905504], [100.0]],\
[46, 'Pd', 'Palladium', [106, 108, 105, 110, 104, 102], [105.903483, 107.903894, 104.905084, 109.905152, 103.904035, 101.905608], [27.33, 26.46, 22.33, 11.72, 11.14, 1.02]],\
[47, 'Ag', 'Silver', [107, 109], [106.905093, 108.904756], [51.839, 48.161]],\
[48, 'Cd', 'Cadmium', [114, 112, 111, 110, 113, 116, 106, 108], [113.903358, 111.902757, 110.904182, 109.903006, 112.904401, 115.904755, 105.906458, 107.904183], [28.73, 24.13, 12.8, 12.49, 12.22, 7.49, 1.25, 0.89]],\
[49, 'In', 'Indium', [115, 113], [114.903878, 112.904061], [95.71, 4.29]],\
[50, 'Sn', 'Tin', [120, 118, 116, 119, 117, 124, 122, 112, 114, 115], [119.902197, 117.901606, 115.901744, 118.903309, 116.902954, 123.905275, 121.90344, 111.904821, 113.902782, 114.903346], [32.58, 24.22, 14.54, 8.59, 7.68, 5.79, 4.63, 0.97, 0.66, 0.34]],\
[51, 'Sb', 'Antimony', [121, 123], [120.903818, 122.904216], [57.21, 42.79]],\
[52, 'Te', 'Tellurium', [130, 128, 126, 125, 124, 122, 123, 120], [129.906223, 127.904461, 125.903306, 124.904425, 123.902819, 121.903047, 122.904273, 119.90402], [34.08, 31.74, 18.84, 7.07, 4.74, 2.55, 0.89, 0.09]],\
[53, 'I', 'Iodine', [127], [126.904468], [100.0]],\
[54, 'Xe', 'Xenon', [132, 129, 131, 134, 136, 130, 128, 124, 126], [131.904154, 128.904779, 130.905082, 133.905395, 135.90722, 129.903508, 127.90353, 123.905896, 125.904269], [26.89, 26.44, 21.18, 10.44, 8.87, 4.08, 1.92, 0.09, 0.09]],\
[55, 'Cs', 'Cesium', [133], [132.905447], [100.0]],\
[56, 'Ba', 'Barium', [138, 137, 136, 135, 134, 130, 132], [137.905241, 136.905821, 135.90457, 134.905683, 133.904503, 129.90631, 131.905056], [71.698, 11.232, 7.854, 6.592, 2.417, 0.106, 0.101]],\
[57, 'La', 'Lanthanum', [139, 138], [138.906348, 137.907107], [99.91, 0.09]],\
[58, 'Ce', 'Cerium', [140, 142, 138, 136], [139.905434, 141.90924, 137.905986, 135.907144], [88.45, 11.114, 0.251, 0.185]],\
[59, 'Pr', 'Praseodymium', [141], [140.907648], [100.0]],\
[60, 'Nd', 'Neodymium', [142, 144, 146, 143, 145, 148, 150], [141.907719, 143.910083, 145.913112, 142.90981, 144.912569, 147.916889, 149.920887], [27.2, 23.8, 17.2, 12.2, 8.3, 5.7, 5.6]],\
[61, 'Pm', 'Promethium', [145], [144.912744], [100.0]],\
[62, 'Sm', 'Samarium', [152, 154, 147, 149, 148, 150, 144], [151.919728, 153.922205, 146.914893, 148.91718, 147.914818, 149.917271, 143.911995], [26.75, 22.75, 14.99, 13.82, 11.24, 7.38, 3.07]],\
[63, 'Eu', 'Europium', [153, 151], [152.921226, 150.919846], [52.19, 47.81]],\
[64, 'Gd', 'Gadolinium', [158, 160, 156, 157, 155, 154, 152], [157.924101, 159.927051, 155.92212, 156.923957, 154.922619, 153.920862, 151.919788], [24.84, 21.86, 20.47, 15.65, 14.8, 2.18, 0.2]],\
[65, 'Tb', 'Terbium', [159], [158.925343], [100.0]],\
[66, 'Dy', 'Dysprosium', [164, 162, 163, 161, 160, 158, 156], [163.929171, 161.926795, 162.928728, 160.92693, 159.925194, 157.924405, 155.924278], [28.18, 25.51, 24.9, 18.91, 2.34, 0.1, 0.06]],\
[67, 'Ho', 'Holmium', [165], [164.930319], [100.0]],\
[68, 'Er', 'Erbium', [166, 168, 167, 170, 164, 162], [165.93029, 167.932368, 166.932045, 169.93546, 163.929197, 161.928775], [33.61, 26.78, 22.93, 14.93, 1.61, 0.14]],\
[69, 'Tm', 'Thulium', [169], [168.934211], [100.0]],\
[70, 'Yb', 'Ytterbium', [174, 172, 173, 171, 176, 170, 168], [173.938858, 171.936378, 172.938207, 170.936322, 175.942568, 169.934759, 167.933894], [31.83, 21.83, 16.13, 14.28, 12.76, 3.04, 0.13]],\
[71, 'Lu', 'Lutetium', [175, 176], [174.940768, 175.942682], [97.41, 2.59]],\
[72, 'Hf', 'Hafnium', [180, 178, 177, 179, 176, 174], [179.946549, 177.943698, 176.94322, 178.945815, 175.941402, 173.94004], [35.08, 27.28, 18.6, 13.62, 5.26, 0.16]],\
[73, 'Ta', 'Tantalum', [181, 180], [180.947996, 179.947466], [99.988, 0.012]],\
[74, 'W', 'Tungsten', [184, 186, 182, 183, 180], [183.950933, 185.954362, 181.948206, 182.950224, 179.946706], [30.64, 28.43, 26.5, 14.31, 0.12]],\
[75, 'Re', 'Rhenium', [187, 185], [186.955751, 184.952956], [62.6, 37.4]],\
[76, 'Os', 'Osmium', [192, 190, 189, 188, 187, 186, 184], [191.961479, 189.958445, 188.958145, 187.955836, 186.955748, 185.953838, 183.952491], [40.78, 26.26, 16.15, 13.24, 1.96, 1.59, 0.02]],\
[77, 'Ir', 'Iridium', [193, 191], [192.962924, 190.960591], [62.7, 37.3]],\
[78, 'Pt', 'Platinum', [195, 194, 196, 198, 192, 190], [194.964774, 193.962664, 195.964935, 197.967876, 191.961035, 189.95993], [33.832, 32.967, 25.242, 7.163, 0.782, 0.014]],\
[79, 'Au', 'Gold', [197], [196.966552], [100.0]],\
[80, 'Hg', 'Mercury', [202, 200, 199, 201, 198, 204, 196], [201.970626, 199.968309, 198.968262, 200.970285, 197.966752, 203.973476, 195.965815], [29.86, 23.1, 16.87, 13.18, 9.97, 6.87, 0.15]],\
[81, 'Tl', 'Thallium', [205, 203], [204.974412, 202.972329], [70.476, 29.524]],\
[82, 'Pb', 'Lead', [208, 206, 207, 204], [207.976636, 205.974449, 206.975881, 203.973029], [52.4, 24.1, 22.1, 1.4]],\
[83, 'Bi', 'Bismuth', [209], [208.980383], [100.0]],\
[84, 'Po', 'Polonium', [209], [208.982416], [100.0]],\
[85, 'At', 'Astatine', [210], [209.987131], [100.0]],\
[86, 'Rn', 'Radon', [222], [222.01757], [100.0]],\
[87, 'Fr', 'Francium', [223], [223.019731], [100.0]],\
[88, 'Ra', 'Radium', [226], [226.025403], [100.0]],\
[89, 'Ac', 'Actinium', [227], [227.027747], [100.0]],\
[90, 'Th', 'Thorium', [232], [232.03805], [100.0]],\
[91, 'Pa', 'Protactinium', [231], [231.035879], [100.0]],\
[92, 'U', 'Uranium', [238, 235, 234], [238.050783, 235.043923, 234.040946], [99.2745, 0.72, 0.0055]],\
[93, 'Np', 'Neptunium', [237], [237.048167], [100.0]],\
[94, 'Pu', 'Plutonium', [244], [244.064198], [100.0]],\
[95, 'Am', 'Americium', [243], [243.061373], [100.0]],\
[96, 'Cm', 'Curium', [247], [247.070347], [100.0]]];



element_max_valences=[
[	1	,	-1	,	1	],
[	2	,	0	,	0	],
[	3	,	0	,	1	],
[	4	,	0	,	2	],
[	5	,	-5	,	3	],
[	6	,	-4	,	4	],
[	7	,	-3	,	5	],
[	8	,	-2	,	2	],
[	9	,	-1	,	0	],
[	10	,	0	,	0	],
[	11	,	-1	,	1	],
[	12	,	0	,	2	],
[	13	,	-2	,	3	],
[	14	,	-4	,	4	],
[	15	,	-3	,	5	],
[	16	,	-2	,	6	],
[	17	,	-1	,	7	],
[	18	,	0	,	0	],
[	19	,	-1	,	1	],
[	20	,	0	,	2	],
[	21	,	0	,	3	],
[	22	,	-2	,	4	],
[	23	,	-3	,	5	],
[	24	,	-4	,	6	],
[	25	,	-3	,	7	],
[	26	,	-4	,	7	],
[	27	,	-3	,	5	],
[	28	,	-2	,	4	],
[	29	,	-2	,	4	],
[	30	,	-2	,	2	],
[	31	,	-5	,	3	],
[	32	,	-4	,	4	],
[	33	,	-3	,	5	],
[	34	,	-2	,	6	],
[	35	,	-1	,	7	],
[	36	,	0	,	2	],
[	37	,	-1	,	1	],
[	38	,	0	,	2	],
[	39	,	0	,	3	],
[	40	,	-2	,	4	],
[	41	,	-3	,	5	],
[	42	,	-4	,	6	],
[	43	,	-3	,	7	],
[	44	,	-4	,	8	],
[	45	,	-3	,	6	],
[	46	,	0	,	6	],
[	47	,	-2	,	4	],
[	48	,	-2	,	2	],
[	49	,	-5	,	3	],
[	50	,	-4	,	4	],
[	51	,	-3	,	5	],
[	52	,	-2	,	6	],
[	53	,	-1	,	7	],
[	54	,	0	,	8	],
[	55	,	-1	,	1	],
[	56	,	0	,	2	],
[	57	,	0	,	3	],
[	58	,	0	,	4	],
[	59	,	0	,	5	],
[	60	,	0	,	4	],
[	61	,	0	,	3	],
[	62	,	0	,	3	],
[	63	,	0	,	3	],
[	64	,	0	,	3	],
[	65	,	0	,	4	],
[	66	,	0	,	4	],
[	67	,	0	,	3	],
[	68	,	0	,	3	],
[	69	,	0	,	3	],
[	70	,	0	,	3	],
[	71	,	0	,	3	],
[	72	,	-2	,	4	],
[	73	,	-3	,	5	],
[	74	,	-4	,	6	],
[	75	,	-3	,	7	],
[	76	,	-4	,	8	],
[	77	,	-3	,	9	],
[	78	,	-3	,	6	],
[	79	,	-3	,	5	],
[	80	,	-2	,	2	],
[	81	,	-5	,	3	],
[	82	,	-4	,	4	],
[	83	,	-3	,	5	],
[	84	,	-2	,	6	],
[	85	,	-1	,	7	],
[	86	,	0	,	6	],
[	87	,	0	,	1	],
[	88	,	0	,	2	],
[	89	,	0	,	3	],
[	90	,	0	,	4	],
[	91	,	0	,	5	],
[	92	,	0	,	6	],
[	93	,	0	,	7	],
[	94	,	0	,	7	],
[	95	,	0	,	7	],
[	96	,	0	,	6	],
];

for s in elements_list_full:
    s.append([]);
    for i in range(len(s[5])):
        s[5][i]=s[5][i]/100.0;
        s[6].append((s[3][i]-s[3][0]));

def get_elements_list_copy():
    return copy.deepcopy(elements_list_full);

elements_list=[];
for s in elements_list_full:
    elements_list.append(s[1]);

elements_dict={};
for s in elements_list_full:
    elements_dict[s[1]]=s[0]-1;

def ppm_to_Da(value):
    return exp(value*1.0e-6)*hydrogen_mass;

def Da_to_ppm(value):
    return log(value/hydrogen_mass)*1.0e6;
    

class Formula(dict):
    def formula_to_string(self):
        s='';
        for key in sorted(self.keys()):
            val=self[key];
            if val>1:
                s='%s%s%s'%(s,key,val);
            else:
                s='%s%s'%(s,key);
                
        return s;

    def __repr__(self):
        return 'Formula("%s")'%self.formula_to_string();
        
    
def parse_formula(formula):
    lst=Formula();
    formula=formula.split('.');
    for fla in formula:
        state=0;
        n=[]; #count element
        s=[]; #name element
        mul=1; #total multiple
        ss=''; #temp substring
        nn=''; #temp subnumber
        for i in fla:
            if i>='A' and i<='Z':
                if state==0:
                    if ss!='':
                        mul=int(ss);
                    else:
                        mul=1;
                elif state==1:
                        s.append(ss);
                        n.append(1);
                elif state==2:
                        s.append(ss);
                        n.append(int(nn));
                state=1;
                ss=i;
                nn='';
            elif i>='a' and i<='z':
                ss+=i;
            elif i>='0' and i<='9':                
                if state==0:
                    ss+=i;
                elif state==1:
                    nn+=i;
                    state=2;
                elif state==2:
                    nn+=i;
                    #state=2;
        if state==0:
            if ss!='':
                mul=int(ss);
            else:
                mul=1;
        elif state==1:
            s.append(ss);
            n.append(1);
        elif state==2:
            s.append(ss);
            n.append(int(nn));            
        for i in range(len(s)):
            if s[i] in lst:
                lst[s[i]]+=n[i]*mul;
            else:
                lst[s[i]]=n[i]*mul;
    return lst;    

def encode_formula_to_array(formula):
    result=np.zeros((96,),dtype=np.uint16);

    for s in formula.keys():
        try:
            result[elements_dict[s]]=formula[s];
        except:
            print('Unknown Element: %s'%s);
    
    return result;

def decode_formula_from_array(formula_array):
    result=Formula();

    for i in range(96):
        if formula_array[i]>0:
            result[elements_list[i]]=formula_array[i];
        
    return result;
            
def formula_to_element_vector(formula):
    elements=np.zeros((96,),dtype=np.uint8);
    for s in formula.keys():
        try:
            elements[elements_dict[s]]=1;
        except:
            print('Unknown Element: %s'%s);
    return np.packbits(elements);

def set_to_element_vector(elementset):
    elements=np.zeros((96,),dtype=np.uint8);
    for s in elementset:
        try:
            elements[elements_dict[s]]=1;
        except:
            print('Unknown Element: %s'%s);
    return np.packbits(elements);

def list_to_element_vector(elementlist):
    elements=np.zeros((96,),dtype=np.uint8);
    for s in elementlist:
        try:
            elements[elements_dict[s]]=1;
        except:
            print('Unknown Element: %s'%s);
    return np.packbits(elements);

def element_vector_to_list(elements):
    elements=np.unpackbits(elements);
    result=[];
    for i in range(96):
        if elements[i]==1:
            result.append(elements_list[i]);
    return result;

def is_within_elements_binary(testelements, refelements):
    m=np.bitwise_or(refelements,testelements);
    result=True;
    for i in range(len(refelements)):
        if m[i]!=refelements[i]:
            result=False;
            break;
    return result;

def not_within_elements_binary(testelements, refelements):
    m=np.bitwise_or(refelements,testelements);
    result=False;
    for i in range(len(refelements)):
        if m[i]!=refelements[i]:
            result=True;
            break;
    return result;
    
    
def get_isotopic_pattern_integer(formula, work_element_list=elements_list_full, rejection_threshold=0.0001):
    result={0:1.0};
    rejection_threshold2=0.0000001;
    for element in formula.keys():
        ecount=formula[element];
        index=elements_dict[element];
        isotopes=work_element_list[index][6];
        abundances=work_element_list[index][5];
        
        for i in range(ecount):
            nextresult={};
            for res in result.keys():
                prevab=result[res];
                for isotope in range(len(isotopes)):
                    ab=prevab*abundances[isotope];
                    if ab>=rejection_threshold2:
                        nextmass=res+isotopes[isotope];
                        if nextmass in nextresult:
                            nextresult[nextmass]=ab+nextresult[nextmass];
                        else:
                            nextresult[nextmass]=ab;                            
            result=nextresult;
    nextresult=[];
    
    for mass in result.keys():
        nextresult.append([mass,result[mass]]);
    
    result=sorted(nextresult,key=itemgetter(0));                        
    
   
    for i in reversed(range(len(result))):
        if result[i][1]<rejection_threshold:    
            del result[i];
    
    return result;




def get_isotopic_pattern(formula, work_element_list=elements_list_full, rejection_threshold=0.0001, condense=False, ppm=20):
    result={0:1.0};
    rejection_threshold2=0.0000001;
    for element in formula.keys():
        ecount=formula[element];
        index=elements_dict[element];
        isotopes=work_element_list[index][4];
        abundances=work_element_list[index][5];
        
        for i in range(ecount):
            nextresult={};
            for res in result.keys():
                prevab=result[res];
                for isotope in range(len(isotopes)):
                    ab=prevab*abundances[isotope];
                    if ab>=rejection_threshold2:
                        nextmass=res+isotopes[isotope];
                        if nextmass in nextresult:
                            nextresult[nextmass]=ab+nextresult[nextmass];
                        else:
                            nextresult[nextmass]=ab;                            
            result=nextresult;
    nextresult=[];
    
    for mass in result.keys():
          nextresult.append([mass,result[mass]]);
    
    result=sorted(nextresult,key=itemgetter(0));                        
    
   
    for i in reversed(range(len(result))):
        if result[i][1]<rejection_threshold:    
            del result[i];
            
    if condense:
        result=sorted(result,key=itemgetter(1),reverse=True);                        
        nextresult=[];
        while result:
            m=0.0;
            intens=0.0;
            for i in reversed(range(len(result))):
                if abs((result[i][0]-result[0][0]))<=result[0][0]/1000000.0*ppm:
                    m+=result[i][0]*result[i][1];
                    intens+=result[i][1];
                    del result[i];
            nextresult.append([m/intens,intens]);
        result=sorted(nextresult,key=itemgetter(0));                        
    
    return result;

  
def get_all_atom_combination_masses(formula, ppm=1.0):
    if 'H' in formula:
        maxH=formula['H'];
    else:
        maxH=0;
    
    workatoms=[];
    workatom_masses=[];
    workatom_max_valence=[];
    
    for key in formula.keys():
        if key!='H':
            workatoms.append(formula[key]);
            index=elements_dict[key];
            workatom_masses.append(elements_list_full[index][4][0]);
            workatom_max_valence.append(abs(element_max_valences[index][1]));

    if maxH>0:
        workatoms.insert(0, maxH);
        workatom_masses.insert(0, hydrogen_mass);
        workatom_max_valence.insert(0, 1);
        
    mass=0.0;
    counts=[];
    for i in range(len(workatoms)):
        mass+=workatom_masses[i]*workatoms[i];
        counts.append(workatoms[i]);
    
    maxppm=Da_to_ppm(mass);
    regions_count=int(maxppm/ppm)+1;
    masses=[0]*regions_count;

    masses[int(Da_to_ppm(mass)/ppm)]+=1;
    totlen=len(workatoms);

    counting=True;
    total_valence=0;
    atcount=0;
 
    if maxH>0:
        for i in range(1, totlen):
            total_valence+=workatom_max_valence[i]*counts[i];
            atcount+=counts[i];
            
    #print('Total valence: %s'%total_valence);
            
    while counting:
        cc=True;
        nextupd=False;
        for i in range(totlen):
            if counts[i]>0:
                counts[i]-=1;
                if maxH>0 and i>0:
                    total_valence-=workatom_max_valence[i];
                    atcount-=1;
                cc=False;
                if nextupd:
                    if maxH>0:                    
                        for j in range(1, i):
                            counts[j]=workatoms[j];
                            atcount+=counts[j];
                            total_valence+=workatom_max_valence[j]*counts[j];
                            
                        
                        #print(total_valence)
                        #print(atcount)
                        if atcount>0:
                            counts[0]=total_valence-(atcount-1);
                            
                        else:
                            counts[0]=total_valence;
                            
                        
                        if counts[0]>0:
                            counts[0]=max(min(counts[0], maxH), 1);
                        else:
                            counts[0]=2;
                        #print(counts[0])
                    else:
                        for j in range(0, i):
                            counts[j]=workatoms[j];
                        
                    
                    nextupd=False
                break;
            else:
                nextupd=True;
        
        if cc:
            counting=False;
        else:
            mass=0.0;
            for i in range(totlen):
                mass+=workatom_masses[i]*counts[i];
            if mass>0.0:
                masses[int(Da_to_ppm(mass)/ppm)]+=1;
            #print(counts)
            
    return masses        
        
    
    
def print_ppm(masses, ppm=1.0):
    print('Masses, rounded down')
    prevppm=0;
    for i in range(len(masses)):
        if masses[i]>0:
            mz=ppm_to_Da(i*ppm);
            print('%8.5f\t%s\t%s'%(mz, masses[i], (i-prevppm)*ppm));
            prevppm=i
    

def save_ppm(fname, masses, ppm=1.0):
    with open(fname, 'w') as fout:
        fout.write('Masses, rounded down\n');
        prevppm=0;
        for i in range(len(masses)):
            if masses[i]>0:
                mz=ppm_to_Da(i*ppm);
                fout.write('%8.5f\t%s\t%s\n'%(mz, masses[i], (i-prevppm)*ppm));
                prevppm=i


if __name__=='__main__':
    '''
    print(parse_formula('2Mn.4O5H3C10'));
    C2H5OH=formula_to_element_vector(parse_formula('C2H5OH'));
    print(C2H5OH);
    C2H5OH=set_to_element_vector(set(['C','H','O','H']));
    print(C2H5OH);
    C2H5OH=list_to_element_vector(['C','H','O','H']);
    print(C2H5OH);
    
    arr=encode_formula_to_array(parse_formula('C2H5OH'));

    print(arr);
    
    print(decode_formula_from_array(arr))
    print(type(decode_formula_from_array(arr)['C']))
    print((parse_formula('C2H5OH').formula_to_string()));
    '''
    
    print(get_isotopic_pattern(parse_formula('C2H5OHK'), work_element_list=elements_list_full, rejection_threshold=0.0001, condense=True, ppm=20));
    
    
    
    #masses=get_all_atom_combination_masses(parse_formula('C2H5OH'), 20.0);
    
    #print_ppm(masses, 20.0);
    #save_ppm('e:/Imperial/testmasses1.list', masses, 20.0);
    
    
    #masses=get_all_atom_combination_masses(parse_formula('C21H44NO7P'), 20.0);
    #print_ppm(masses, 20.0);
    #save_ppm('e:/Imperial/testmasses2.list', masses, 20.0);

    #masses=get_all_atom_combination_masses(parse_formula('C55H75N17O13'), 20.0);
    #print_ppm(masses, 20.0);
    #save_ppm('e:/Imperial/testmasses3.list', masses, 20.0);







