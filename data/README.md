# Test Set of MS2 spectra for ChemDistiller

* This is a set of 867 compounds with their respective MS2 spectra from MassBank (http://www.massbank.jp/). 
* The spectra were ported into the ChemDistiller (https://bitbucket.org/iAnalytica/chemdistillerpython) tree-like text format for spectral data and can be used directly to test its functionality. 
* MS2 spectra had their corresponding MS1 spectra generated with ideal precursor MS1 peak m/z values. Otherwise, the MS2 spectra were not altered.

* These spectra and compounds were pre-selected based on the conditions of its license being CC-BY (https://creativecommons.org/licenses/by/3.0/), thus making them free to re-use, even for commercial purposes, provided the proper attribution is made. 
* These spectra were additionally prefiltered to have estimated precision of at least 0.005 Da for the peak with the highest m/z value.

To use these spectra you will need to unpack them into an _empty_ folder and load into the Spectral Manager in ChemDistiller. 

Example:
```python
#Test spectra are unpacked into the test_spectral_database_path folder
from chemdistiller.msspectra.manager import SpectralManager;
sp_manager = SpectralManager();
sp_manager.import_textfile_spectra_from_folder(test_spectral_database_path);
```

One can sequentially load spectra from several locations into the Spectral Manager in ChemDistiller:
```python
from chemdistiller.msspectra.manager import SpectralManager;
sp_manager = SpectralManager();
sp_manager.import_textfile_spectra_from_folder(path1);
sp_manager.import_textfile_spectra_from_folder(path2);
sp_manager.import_textfile_spectra_from_folder(path3);
```
__Note!__ Current algorithm searches for any __*.txt__ files in the folder and tries to import them. To avoid conflicts the folder with spectra should not contain any other files with the extension __*.txt__ ! For the future version compatability, ideally, there should _only_ be spectral files in the spectral folder. Please do not mix these files with others in the same folder.


For more information on the use of the engine check out ChemDistiller documentation (https://bitbucket.org/iAnalytica/chemdistillerpython)

---
#### The following [MassBank][1] entries were used to generate this set of data files:
[1]: http://www.massbank.jp/

Authors: K. Wilkinson, S. Miranda   
CopyRight: UPAO   

	Molecule: FLAVOXATE
    Entries: UPA00001

Authors: Otto J, Stravs M, Schymanski E, Singer H, Department of Environmental Chemistry, Eawag   
CopyRight: Copyright (C) 2015 Eawag, Duebendorf, Switzerland   

    Molecule: Niclosamide 5-chloro-N-(2-chloro-4-nitrophenyl)-2-hydroxybenzamide
    Entries: EQ332001, EQ332002, EQ332003, EQ332004, EQ332005, EQ332006, EQ332051, EQ332052, EQ332053, EQ332054, EQ332055, EQ332056
    Molecule: Zidovudine 1-[(2R,4S,5S)-4-azido-5-(hydroxymethyl)oxolan-2-yl]-5-methylpyrimidine-2,4-dione
    Entries: EQ328701, EQ328702, EQ328703, EQ328704, EQ328705, EQ328706, EQ328751, EQ328752, EQ328753, EQ328754, EQ328755, EQ328756
    Molecule: Doxazosin [4-(4-amino-6,7-dimethoxyquinazolin-2-yl)piperazin-1-yl]-(2,3-dihydro-1,4-benzodioxin-3-yl)methanone
    Entries: EQ329301, EQ329302, EQ329303, EQ329304, EQ329305, EQ329306, EQ329351, EQ329352, EQ329353, EQ329354, EQ329355, EQ329356
    Molecule: Linezolid N-[[(5S)-3-(3-fluoro-4-morpholin-4-ylphenyl)-2-oxo-1,3-oxazolidin-5-yl]methyl]acetamide
    Entries: EQ329101, EQ329102, EQ329103, EQ329104, EQ329105, EQ329106
    Molecule: N-ethyl-4-methoxybenzamide
    Entries: EQ333701, EQ333702, EQ333703, EQ333704, EQ333705, EQ333706
    Molecule: Bufexamac 2-(4-butoxyphenyl)-N-hydroxyacetamide
    Entries: EQ332201, EQ332202, EQ332203, EQ332204, EQ332205, EQ332206, EQ332251, EQ332252, EQ332253, EQ332254, EQ332255, EQ332256
    Molecule: Tilidine ethyl (2R)-2-(dimethylamino)-1-phenylcyclohex-3-ene-1-carboxylate
    Entries: EQ332901, EQ332902, EQ332903, EQ332904, EQ332905, EQ332906
    Molecule: Olopatadine 2-[(11Z)-11-[3-(dimethylazaniumyl)propylidene]-6H-benzo[c][1]benzoxepin-2-yl]acetate
    Entries: EQ332301, EQ332302, EQ332303, EQ332304, EQ332305, EQ332306, EQ332351, EQ332352, EQ332353, EQ332354, EQ332355, EQ332356
    Molecule: Dienogest 2-[(8S,13S,14S,17R)-17-hydroxy-13-methyl-3-oxo-1,2,6,7,8,11,12,14,15,16-decahydrocyclopenta[a]phenanthren-17-yl]acetonitrile
    Entries: EQ328401, EQ328402, EQ328403, EQ328404, EQ328405, EQ328406, EQ328451, EQ328452, EQ328453, EQ328454, EQ328455, EQ328456
    Molecule: Hydrocodone (4R,4aR,7aR,12bS)-9-methoxy-3-methyl-1,2,4,4a,5,6,7a,13-octahydro-4,12-methanobenzofuro[3,2-e]isoquinoline-7-one
    Entries: EQ333201, EQ333202, EQ333203, EQ333204, EQ333205, EQ333206
    Molecule: Maprotiline
    Entries: EQ331701, EQ331702, EQ331703, EQ331704, EQ331705, EQ331706
    Molecule: Carisoprodol [2-(carbamoyloxymethyl)-2-methylpentyl] N-propan-2-ylcarbamate
    Entries: EQ332701, EQ332702, EQ332703, EQ332704, EQ332705, EQ332706
    Molecule: Benproperine 1-[1-(2-benzylphenoxy)propan-2-yl]piperidine
    Entries: EQ329601, EQ329602, EQ329603, EQ329604, EQ329605, EQ329606
    Molecule: Bexarotene 4-[1-(3,5,5,8,8-pentamethyl-6,7-dihydronaphthalen-2-yl)ethenyl]benzoic acid
    Entries: EQ329501, EQ329502, EQ329503, EQ329504, EQ329505, EQ329506, EQ329551, EQ329552, EQ329553, EQ329554
    Molecule: Dropropizine 3-(4-phenylpiperazin-1-yl)propane-1,2-diol
    Entries: EQ330701, EQ330702, EQ330703, EQ330704, EQ330705, EQ330706
    Molecule: Atropine (1R,5S)-8-Methyl-8-azabicyclo[3.2.1]oct-3-yl tropate
    Entries: EQ333401, EQ333402, EQ333403, EQ333404, EQ333405, EQ333406
    Molecule: Meptazinol 3-(3-ethyl-1-methylazepan-3-yl)phenol
    Entries: EQ329001, EQ329002, EQ329003, EQ329004, EQ329005, EQ329006
    Molecule: Pioglitazone 5-[[4-[2-(5-ethylpyridin-2-yl)ethoxy]phenyl]methyl]-1,3-thiazolidine-2,4-dione
    Entries: EQ328601, EQ328602, EQ328603, EQ328604, EQ328605, EQ328606, EQ328651, EQ328652, EQ328653, EQ328654, EQ328655, EQ328656
    Molecule: Nateglinide
    Entries: EQ328901, EQ328902, EQ328903, EQ328904, EQ328905, EQ328906, EQ328951, EQ328952, EQ328953, EQ328954, EQ328955, EQ328956
    Molecule: Efavirenz (4S)-6-chloro-4-(2-cyclopropylethynyl)-4-(trifluoromethyl)-1H-3,1-benzoxazin-2-one
    Entries: EQ329201, EQ329202, EQ329203, EQ329204, EQ329205, EQ329206, EQ329251, EQ329252, EQ329253, EQ329254, EQ329255, EQ329256
    Molecule: Praziquantel 2-(cyclohexanecarbonyl)-3,6,7,11b-tetrahydro-1H-pyrazino[2,1-a]isoquinolin-4-one
    Entries: EQ327201, EQ327202, EQ327203, EQ327204, EQ327205, EQ327206
    Molecule: Meperidine ethyl 1-methyl-4-phenylpiperidine-4-carboxylate
    Entries: EQ333101, EQ333102, EQ333103, EQ333104, EQ333105, EQ333106
    Molecule: Penciclovir 2-amino-9-[4-hydroxy-3-(hydroxymethyl)butyl]-3H-purin-6-one
    Entries: EQ328801, EQ328802, EQ328803, EQ328804, EQ328805, EQ328806, EQ328851, EQ328852, EQ328853, EQ328854, EQ328855, EQ328856
    Molecule: Dibucaine 2-butoxy-N-[2-(diethylamino)ethyl]quinoline-4-carboxamide
    Entries: EQ329401, EQ329402, EQ329403, EQ329404, EQ329405, EQ329406
    Molecule: Bromochlorophen 2-bromo-6-[(3-bromo-5-chloro-2-hydroxyphenyl)methyl]-4-chlorophenol
    Entries: EQ329851, EQ329852, EQ329853, EQ329854, EQ329855, EQ329856
    Molecule: Tetrazepam 7-chloro-5-(cyclohexen-1-yl)-1-methyl-3H-1,4-benzodiazepin-2-one
    Entries: EQ333001, EQ333002, EQ333003, EQ333004, EQ333005, EQ333006
    Molecule: Bupivacaine 1-butyl-N-(2,6-dimethylphenyl)piperidine-2-carboxamide
    Entries: EQ330501, EQ330502, EQ330503, EQ330504, EQ330505, EQ330506
    Molecule: Etodolac 2-(1,8-diethyl-4,9-dihydro-3H-pyrano[3,4-b]indol-1-yl)acetic acid
    Entries: EQ330801, EQ330802, EQ330803, EQ330804, EQ330805, EQ330806, EQ330851, EQ330852, EQ330853, EQ330854, EQ330855, EQ330856
    Molecule: Dihydrocodeine (4R,4aR,7S,7aR,12bS)-9-methoxy-3-methyl-2,4,4a,5,6,7,7a,13-octahydro-1H-4,12-methanobenzofuro[3,2-e]isoquinoline-7-ol
    Entries: EQ333301, EQ333302, EQ333303, EQ333304, EQ333305, EQ333306
    Molecule: Allopurinol 1,2-dihydropyrazolo[3,4-d]pyrimidin-4-one
    Entries: EQ333601, EQ333602, EQ333603, EQ333604, EQ333605, EQ333606, EQ333651, EQ333652, EQ333653, EQ333654, EQ333655, EQ333656
    Molecule: Cathine (1S,2S)-2-amino-1-phenylpropan-1-ol
    Entries: EQ333501, EQ333502, EQ333503, EQ333504, EQ333505, EQ333506
    Molecule: Linoleic acid (9Z,12Z)-octadeca-9,12-dienoic acid
    Entries: EQ331601, EQ331602, EQ331603, EQ331604, EQ331605, EQ331606, EQ331651, EQ331652, EQ331653, EQ331654, EQ331655

Authors: Nikiforos Alygizakis, Anna Bletsou, Nikolaos Thomaidis, University of Athens   
CopyRight: Copyright (C) 2015 Department of Chemistry, University of Athens   

    Molecule: Sulfamethizole 4-amino-N-(5-methyl-1,3,4-thiadiazol-2-yl)benzenesulfonamide
    Entries: AU101702
    Molecule: Sulfachlorpyridazine 4-amino-N-(6-chloropyridazin-3-yl)benzenesulfonamide
    Entries: AU100701
    Molecule: Atorvastatin (3R,5R)-7-[2-(4-fluorophenyl)-3-phenyl-4-(phenylcarbamoyl)-5-propan-2-ylpyrrol-1-yl]-3,5-dihydroxyheptanoic acid
    Entries: AU112901
    Molecule: Danofloxacin 1-cyclopropyl-6-fluoro-7-[(1S,4S)-5-methyl-2,5-diazabicyclo[2.2.1]heptan-2-yl]-4-oxoquinoline-3-carboxylic acid
    Entries: AU102701
    Molecule: Sarafloxacin 6-fluoro-1-(4-fluorophenyl)-4-oxo-7-piperazin-4-ium-1-ylquinoline-3-carboxylate
    Entries: AU103501
    Molecule: Mabuterol 1-[4-amino-3-chloro-5-(trifluoromethyl)phenyl]-2-(tert-butylamino)ethanol
    Entries: AU110201
    Molecule: Flubendazole methyl N-[6-(4-fluorobenzoyl)-1H-benzimidazol-2-yl]carbamate
    Entries: AU116304
    Molecule: Sulfadimidine 4-amino-N-(4,6-dimethylpyrimidin-2-yl)benzenesulfonamide
    Entries: AU100801
    Molecule: Sulfadoxine 4-amino-N-(5,6-dimethoxypyrimidin-4-yl)benzenesulfonamide
    Entries: AU101001
    Molecule: Robenidine 1,2-bis[(E)-(4-chlorophenyl)methylideneamino]guanidine
    Entries: AU119301
    Molecule: Mebendazole methyl N-(6-benzoyl-1H-benzimidazol-2-yl)carbamate
    Entries: AU116401
    Molecule: Ofloxacin 9-Fluoro-3-methyl-10-(4-methyl-1-piperazinyl)-7-oxo-2,3-dihydro-7H-[1,4]oxazino[2,3,4-ij]quinoline-6-carboxylic acid
    Entries: AU103301
    Molecule: Oxfendazole methyl N-[6-(benzenesulfinyl)-1H-benzimidazol-2-yl]carbamate
    Entries: AU116502
    Molecule: Sulfamerazine 4-amino-N-(4-methylpyrimidin-2-yl)benzenesulfonamide
    Entries: AU101601
    Molecule: Ceftiofur (6R,7R)-7-[[(2Z)-2-(2-amino-1,3-thiazol-4-yl)-2-methoxyiminoacetyl]amino]-3-(furan-2-carbonylsulfanylmethyl)-8-oxo-5-thia-1-azabicyclo[4.2.0]oct-2-ene-2-carboxylic acid
    Entries: AU105101
    Molecule: Sulfaguanidine 2-(4-aminophenyl)sulfonylguanidine
    Entries: AU101201
    Molecule: Niflumic acid 2-[3-(trifluoromethyl)anilino]pyridine-3-carboxylic acid
    Entries: AU115401
    Molecule: Triclabendazole 6-chloro-5-(2,3-dichlorophenoxy)-2-methylsulfanyl-1H-benzimidazole
    Entries: AU106702
    Molecule: Lincomycin (2S,4R)-N-[(1R,2R)-2-hydroxy-1-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-methylsulfanyloxan-2-yl]propyl]-1-methyl-4-propylpyrrolidine-2-carboxamide
    Entries: AU118001
    Molecule: Ciprofloxacin 1-cyclopropyl-6-fluoro-4-oxo-7-piperazin-1-ylquinoline-3-carboxylic acid
    Entries: AU102603
    Molecule: Sulfisoxazole 4-amino-N-(3,4-dimethyl-1,2-oxazol-5-yl)benzenesulfonamide
    Entries: AU101302
    Molecule: Thiabendazole 4-(1H-benzimidazol-2-yl)-1,3-thiazole
    Entries: AU106601
    Molecule: Oxolinic acid 5-ethyl-8-oxo-[1,3]dioxolo[4,5-g]quinoline-7-carboxylic acid
    Entries: AU103401
    Molecule: Atenolol 2-[4-[2-hydroxy-3-(propan-2-ylamino)propoxy]phenyl]acetamide
    Entries: AU110601
    Molecule: Diaveridine 5-[(3,4-dimethoxyphenyl)methyl]pyrimidine-2,4-diamine
    Entries: AU108501
    Molecule: Albendazole methyl N-(6-propylsulfanyl-1H-benzimidazol-2-yl)carbamate
    Entries: AU105801
    Molecule: Ronidazole (1-methyl-5-nitroimidazol-2-yl)methyl carbamate
    Entries: AU107101
    Molecule: Valsartan (2S)-3-methyl-2-[pentanoyl-[[4-[2-(2H-tetrazol-5-yl)phenyl]phenyl]methyl]amino]butanoic acid
    Entries: AU111202
    Molecule: Difloxacin 6-fluoro-1-(4-fluorophenyl)-7-(4-methylpiperazin-4-ium-1-yl)-4-oxoquinoline-3-carboxylate
    Entries: AU102801
    Molecule: Paracetamol N-(4-hydroxyphenyl)acetamide
    Entries: AU112601
    Molecule: Sulfamoxole 4-amino-N-(4,5-dimethyl-1,3-oxazol-2-yl)benzenesulfonamide
    Entries: AU101902
    Molecule: Sulfamethoxazole 4-amino-N-(5-methyl-1,2-oxazol-3-yl)benzenesulfonamide
    Entries: AU101801
    Molecule: Tiamulin
    Entries: AU105501
    Molecule: Clopidol 3,5-dichloro-2,6-dimethyl-1H-pyridin-4-one
    Entries: AU118301
    Molecule: Clenbuterol 1-(4-amino-3,5-dichlorophenyl)-2-(tert-butylamino)ethanol
    Entries: AU109901
    Molecule: Ractopamine 4-[3-[[2-hydroxy-2-(4-hydroxyphenyl)ethyl]amino]butyl]phenol
    Entries: AU110301
    Molecule: Morantel 1-methyl-2-[(E)-2-(3-methylthiophen-2-yl)ethenyl]-5,6-dihydro-4H-pyrimidine
    Entries: AU113801
    Molecule: Ternidazole 3-(2-methyl-5-nitroimidazol-1-yl)propan-1-ol
    Entries: AU107001
    Molecule: Progesterone (8S,9S,10R,13S,14S,17S)-17-acetyl-10,13-dimethyl-1,2,6,7,8,9,11,12,14,15,16,17-dodecahydrocyclopenta[a]phenanthren-3-one
    Entries: AU107701
    Molecule: Enrofloxacin 1-cyclopropyl-7-(4-ethylpiperazin-4-ium-1-yl)-6-fluoro-4-oxoquinoline-3-carboxylate
    Entries: AU102901
    Molecule: Ethopabate methyl 4-acetamido-2-ethoxybenzoate
    Entries: AU118402
    Molecule: Propranolol 1-naphthalen-1-yloxy-3-(propan-2-ylamino)propan-2-ol
    Entries: AU110801
    Molecule: Imidocarb 1,3-bis[3-(4,5-dihydro-1H-imidazol-2-yl)phenyl]urea
    Entries: AU117003
    Molecule: Norfloxacin 1-ethyl-6-fluoro-4-oxo-7-piperazin-4-ium-1-ylquinoline-3-carboxylate
    Entries: AU103202
    Molecule: Metoprolol 1-[4-(2-methoxyethyl)phenoxy]-3-(propan-2-ylamino)propan-2-ol
    Entries: AU110701
    Molecule: Sulfaclozine 4-amino-N-(6-chloropyrazin-2-yl)benzenesulfonamide
    Entries: AU100601
    Molecule: Ranitidine 1-N'-[2-[[5-[(dimethylamino)methyl]furan-2-yl]methylsulfanyl]ethyl]-1-N-methyl-2-nitroethene-1,1-diamine
    Entries: AU111501
    Molecule: Sulfameter 4-amino-N-(5-methoxypyrimidin-2-yl)benzenesulfonamide
    Entries: AU102501
    Molecule: Flumequine 9-Fluoro-5-methyl-1-oxo-6,7-dihydro-1H,5H-pyrido[3,2,1-ij]quinoline-2-carboxylic acid
    Entries: AU103001
    Molecule: Marbofloxacin 9-Fluoro-3-methyl-10-(4-methyl-1-piperazinyl)-7-oxo-2,3-dihydro-7H-[1,3,4]oxadiazino[6,5,4-ij]quinoline-6-carboxylic acid
    Entries: AU103102
    Molecule: Albendazole Sulfone methyl N-(6-propylsulfonyl-1H-benzimidazol-2-yl)carbamate
    Entries: AU105901
    Molecule: Salbutamol 4-[2-(tert-butylamino)-1-hydroxyethyl]-2-(hydroxymethyl)phenol
    Entries: AU110002

Authors: K. A. Wilkinson, S. N. Miranda   
CopyRight: UPAO   

    Molecule: RIVAROXABAN
    Entries: UPA00012

Authors: K.A. Wilkinson & S.N. Miranda   
CopyRight: UPAO   

    Molecule: TELMISARTAN
    Entries: UPA00022
    Molecule: TRIMETHOPRIM
    Entries: UPA00020
    Molecule: CANDESARTAN CILEXETIL
    Entries: UPA00014
    Molecule: VALSARTAN
    Entries: UPA00017
    Molecule: OLMESARTAN MEDOXOMIL
    Entries: UPA00016
    Molecule: LORATADINE
    Entries: UPA00015
    Molecule: 8-CHLOROTHEOPHYLLINE
    Entries: UPA00023
    Molecule: ROSUVASTATIN
    Entries: UPA00013
    Molecule: KETOROLAC
    Entries: UPA00021

Authors: Stravs M, Schymanski E, Singer H, Department of Environmental Chemistry, Eawag and Thomaidis N, University of Athens   
CopyRight: Copyright (C) 2015 Eawag, Duebendorf, Switzerland   

    Molecule: Chlorthal-dimethyl Dimethyl 2,3,5,6-tetrachloroterephthalate Dimethyl 2,3,5,6-tetrachlorobenzene-1,4-dicarboxylate
    Entries: EQ370101, EQ370102, EQ370103, EQ370104, EQ370105, EQ370106, EQ370107, EQ370108, EQ370109
    Molecule: Triallate S-(2,3,3-trichloroprop-2-enyl) N,N-di(propan-2-yl)carbamothioate
    Entries: EQ372501, EQ372502, EQ372503, EQ372504, EQ372505, EQ372506, EQ372507, EQ372508, EQ372509
    Molecule: Tetraconazole 1-[2-(2,4-dichlorophenyl)-3-(1,1,2,2-tetrafluoroethoxy)propyl]-1,2,4-triazole
    Entries: EQ372301, EQ372302, EQ372303, EQ372304, EQ372305, EQ372306, EQ372307, EQ372308, EQ372309
    Molecule: Phenylbutazone 4-butyl-1,2-diphenylpyrazolidine-3,5-dione
    Entries: EQ361501, EQ361502, EQ361503, EQ361504, EQ361505, EQ361506, EQ361507, EQ361508, EQ361509, EQ361551, EQ361552, EQ361553, EQ361554, EQ361555, EQ361556, EQ361557, EQ361558, EQ361559
    Molecule: Bromadiolone 3-[3-[4-(4-bromophenyl)phenyl]-3-hydroxy-1-phenylpropyl]-4-hydroxychromen-2-one
    Entries: EQ365651, EQ365652, EQ365653, EQ365654, EQ365655, EQ365656, EQ365657, EQ365658, EQ365659
    Molecule: Imazaquin 2-(4-methyl-5-oxo-4-propan-2-yl-1H-imidazol-2-yl)quinoline-3-carboxylic acid
    Entries: EQ371101, EQ371102, EQ371103, EQ371104, EQ371105, EQ371106, EQ371107, EQ371108, EQ371109, EQ371151, EQ371152, EQ371153, EQ371154, EQ371155, EQ371156, EQ371157, EQ371158
    Molecule: Meclofenamic Acid 2-(2,6-dichloro-3-methylanilino)benzoic acid
    Entries: EQ369001, EQ369002, EQ369003, EQ369004, EQ369005, EQ369006, EQ369007, EQ369008, EQ369009, EQ369051, EQ369052, EQ369053, EQ369054, EQ369055, EQ369056, EQ369057
    Molecule: Profoxydim 2-[1-[2-(4-chlorophenoxy)propoxyamino]butylidene]-5-(thian-3-yl)cyclohexane-1,3-dione
    Entries: EQ360451, EQ360452, EQ360453, EQ360454, EQ360455, EQ360456, EQ360457, EQ360458, EQ360459
    Molecule: Nortriptyline 3-(5,6-dihydrodibenzo[2,1-b:2`,1`-f][7]annulen-11-ylidene)-N-methylpropan-1-amine
    Entries: EQ369201, EQ369202, EQ369203, EQ369204, EQ369205, EQ369206, EQ369207, EQ369208, EQ369209
    Molecule: Difloxacin 6-fluoro-1-(4-fluorophenyl)-7-(4-methylpiperazin-4-ium-1-yl)-4-oxoquinoline-3-carboxylate
    Entries: EQ366601, EQ366602, EQ366603, EQ366604, EQ366605, EQ366606, EQ366607, EQ366608, EQ366609
    Molecule: Topiramate Topamax (2,2,7,7-tetramethyl-5,5a,8a,8b-tetrahydrodi[1,3]dioxolo[4,5-a:5`,3`-d]pyran-3a-yl)methyl sulfamate
    Entries: EQ363501, EQ363502, EQ363503, EQ363504, EQ363505, EQ363506, EQ363507, EQ363508, EQ363509, EQ363551, EQ363552, EQ363553, EQ363554, EQ363555, EQ363556, EQ363557, EQ363558, EQ363559
    Molecule: Norclozapine 3-chloro-6-piperazin-1-yl-5H-benzo[b][1,4]benzodiazepine
    Entries: EQ360701, EQ360702, EQ360703, EQ360704, EQ360705, EQ360706, EQ360707, EQ360708, EQ360709
    Molecule: Proquinazid 6-Iodo-2-propoxy-3-propylquinazolin-4-one
    Entries: EQ349101, EQ349102, EQ349103, EQ349104, EQ349105, EQ349106, EQ349107, EQ349108, EQ349109
    Molecule: Molinate S-ethyl azepane-1-carbothioate
    Entries: EQ371401, EQ371402, EQ371403, EQ371404, EQ371405, EQ371406, EQ371407, EQ371408, EQ371409
    Molecule: Di-n-butyl phthalate Dibutyl phthalate dibutyl benzene-1,2-dicarboxylate
    Entries: EQ367001, EQ367002, EQ367003, EQ367004, EQ367005, EQ367006, EQ367007, EQ367008, EQ367009, EQ367051, EQ367052, EQ367053, EQ367054, EQ367055, EQ367056, EQ367057, EQ367058, EQ367059
    Molecule: Triamterene 6-phenylpteridine-2,4,7-triamine
    Entries: EQ363701, EQ363702, EQ363703, EQ363704, EQ363705, EQ363706, EQ363707, EQ363708, EQ363709
    Molecule: 4-Toluenesulfonamide p-Toluenesulfonamide 4-methylbenzenesulfonimidic acid
    Entries: EQ361801, EQ361802, EQ361803, EQ361804, EQ361805, EQ361806, EQ361807, EQ361808, EQ361809, EQ361851, EQ361852, EQ361853, EQ361854, EQ361855, EQ361856, EQ361857, EQ361858, EQ361859
    Molecule: Fluometuron 1,1-dimethyl-3-[3-(trifluoromethyl)phenyl]urea
    Entries: EQ370901, EQ370902, EQ370903, EQ370904, EQ370905, EQ370906, EQ370907, EQ370908, EQ370909, EQ370951, EQ370952, EQ370953, EQ370954, EQ370955, EQ370956, EQ370957, EQ370958, EQ370959
    Molecule: Thiamphenicol 2,2-dichloro-N-[1,3-dihydroxy-1-(4-methylsulfonylphenyl)propan-2-yl]acetamide
    Entries: EQ363401, EQ363402, EQ363403, EQ363404, EQ363405, EQ363406, EQ363407, EQ363408, EQ363409, EQ363451, EQ363452, EQ363453, EQ363454, EQ363455, EQ363456, EQ363457, EQ363458, EQ363459
    Molecule: Olanzapine 2-methyl-4-(4-methylpiperazin-1-yl)-5H-thieno[3,2-c][1,5]benzodiazepine
    Entries: EQ369701, EQ369702, EQ369703, EQ369704, EQ369705, EQ369706, EQ369707, EQ369708, EQ369709
    Molecule: Isoxaben 2,6-dimethoxy-N-[3-(3-methylpentan-3-yl)-1,2-oxazol-5-yl]benzamide
    Entries: EQ360201, EQ360202, EQ360203, EQ360204, EQ360205, EQ360206, EQ360207, EQ360208, EQ360209, EQ360251, EQ360252, EQ360253, EQ360254, EQ360255, EQ360256, EQ360257, EQ360258, EQ360259
    Molecule: Norfentanyl N-phenyl-N-piperidin-4-ylpropanamide
    Entries: EQ361201, EQ361202, EQ361203, EQ361204, EQ361205, EQ361206, EQ361207, EQ361208, EQ361209
    Molecule: Perfluorobutane sulfonic acid (PFBuS) Perfluorobutanesulfonic acid 1,1,2,2,3,3,4,4,4-nonafluorobutane-1-sulfonic acid
    Entries: EQ367551, EQ367552, EQ367553, EQ367554, EQ367555, EQ367556, EQ367557, EQ367558, EQ367559
    Molecule: Niflumic acid 2-[3-(trifluoromethyl)anilino]pyridine-3-carboxylic acid
    Entries: EQ369101, EQ369102, EQ369103, EQ369104, EQ369105, EQ369106, EQ369107, EQ369108, EQ369109, EQ369151, EQ369152, EQ369153, EQ369154, EQ369155, EQ369156, EQ369157, EQ369158, EQ369159
    Molecule: Flurtamone 5-(methylamino)-2-phenyl-4-[3-(trifluoromethyl)phenyl]furan-3-one
    Entries: EQ360001, EQ360002, EQ360003, EQ360004, EQ360005, EQ360006, EQ360007, EQ360008, EQ360009, EQ360051, EQ360052, EQ360053, EQ360054, EQ360055, EQ360056, EQ360057, EQ360058, EQ360059
    Molecule: Flubendazole Methyl N-[6-(4-fluorobenzoyl)-1H-benzimidazol-2-yl]carbamate
    Entries: EQ365901, EQ365902, EQ365903, EQ365904, EQ365905, EQ365906, EQ365907, EQ365908, EQ365909, EQ365951, EQ365952, EQ365953, EQ365954, EQ365955, EQ365956, EQ365957, EQ365958, EQ365959
    Molecule: Fenthion-sulfone Fenthione sulphone O,O-Dimethyl O-[3-methyl-4-(methylsulfonyl)phenyl] phosphorothioate
    Entries: EQ359901, EQ359902, EQ359903, EQ359904, EQ359905, EQ359906, EQ359907, EQ359908, EQ359909
    Molecule: Ethoprop 1-[ethoxy(propylsulfanyl)phosphoryl]sulfanylpropane
    Entries: EQ365101, EQ365102, EQ365103, EQ365104, EQ365105, EQ365106, EQ365107, EQ365108, EQ365109
    Molecule: Fenthion-sulfoxide Mesulfenos O,O-Dimethyl O-[3-methyl-4-(methylsulfinyl)phenyl] phosphorothioate
    Entries: EQ370701, EQ370702, EQ370703, EQ370704, EQ370705, EQ370706, EQ370707, EQ370708, EQ370709
    Molecule: 2,3,4,6-Tetrachlorophenol 1-Hydroxy-2,3,4, 6-tetrachlorobenzene
    Entries: EQ357551, EQ357552, EQ357553, EQ357554, EQ357555, EQ357556
    Molecule: Cimetidine 1-cyano-2-methyl-3-[2-[(5-methyl-1H-imidazol-4-yl)methylsulfanyl]ethyl]guanidine
    Entries: EQ372201, EQ372202, EQ372203, EQ372204, EQ372205, EQ372206, EQ372207, EQ372208, EQ372209, EQ372251, EQ372252, EQ372253, EQ372254, EQ372255, EQ372256, EQ372257, EQ372258, EQ372259
    Molecule: Paclobutrazole Paclobutrazol 1-(4-chlorophenyl)-4,4-dimethyl-2-(1,2,4-triazol-1-yl)pentan-3-ol
    Entries: EQ370501, EQ370502, EQ370503, EQ370504, EQ370505, EQ370506, EQ370507, EQ370508, EQ370509
    Molecule: Droperidol 3-[1-[4-(4-fluorophenyl)-4-oxobutyl]-3,6-dihydro-2H-pyridin-4-yl]-1H-benzimidazol-2-one
    Entries: EQ367901, EQ367902, EQ367903, EQ367904, EQ367905, EQ367906, EQ367907, EQ367908, EQ367909
    Molecule: Clenbuterol 1-(4-amino-3,5-dichlorophenyl)-2-(tert-butylamino)ethanol
    Entries: EQ359201, EQ359202, EQ359203, EQ359204, EQ359205, EQ359206, EQ359207, EQ359208, EQ359209
    Molecule: Haloxyfop 2-[4-[3-chloro-5-(trifluoromethyl)pyridin-2-yl]oxyphenoxy]propanoic acid
    Entries: EQ356201, EQ356202, EQ356203, EQ356204, EQ356205, EQ356206, EQ356207, EQ356208, EQ356209, EQ356251, EQ356252, EQ356253, EQ356254, EQ356255, EQ356256, EQ356257, EQ356258, EQ356259
    Molecule: Quinmerac 7-chloro-3-methyl-quinoline-8-carboxylic acid 7-chloro-3-methylquinoline-8-carboxylic acid
    Entries: EQ372001, EQ372002, EQ372003, EQ372004, EQ372005, EQ372006, EQ372007, EQ372008, EQ372009
    Molecule: Clomipramine 3-(2-chloro-5,6-dihydrobenzo[b][1]benzazepin-11-yl)-N,N-dimethylpropan-1-amine
    Entries: EQ366901, EQ366902, EQ366903, EQ366904, EQ366905, EQ366906, EQ366907, EQ366908, EQ366909
    Molecule: Enrofloxacin 1-cyclopropyl-7-(4-ethylpiperazin-4-ium-1-yl)-6-fluoro-4-oxoquinoline-3-carboxylate
    Entries: EQ366701, EQ366702, EQ366703, EQ366704, EQ366705, EQ366706, EQ366707, EQ366708, EQ366709
    Molecule: 2-Naphthoxyacetic acid beta-Naphthoxyacetic acid (2-Naphthyloxy)acetic acid
    Entries: EQ371601, EQ371602, EQ371603, EQ371604, EQ371605, EQ371606, EQ371607, EQ371608, EQ371609, EQ371651, EQ371652, EQ371653, EQ371654, EQ371655, EQ371656, EQ371657, EQ371658, EQ371659
    Molecule: 2-Toluenesulfonamide o-Toluenesulfonamide 2-methylbenzenesulfonamide
    Entries: EQ360501, EQ360502, EQ360503, EQ360504, EQ360505, EQ360506, EQ360507, EQ360508, EQ360509, EQ360551, EQ360552, EQ360553, EQ360554, EQ360555, EQ360556, EQ360557, EQ360558, EQ360559
    Molecule: Diaveridine 5-[(3,4-dimethoxyphenyl)methyl]pyrimidine-2,4-diamine
    Entries: EQ364501, EQ364502, EQ364503, EQ364504, EQ364505, EQ364506, EQ364507, EQ364508, EQ364509
    Molecule: Flutriafol 1-(2-fluorophenyl)-1-(4-fluorophenyl)-2-(1,2,4-triazol-1-yl)ethanol
    Entries: EQ371001, EQ371002, EQ371003, EQ371004, EQ371005, EQ371006, EQ371007, EQ371008, EQ371009
    Molecule: Nitrazepam 7-nitro-5-phenyl-1,3-dihydro-1,4-benzodiazepin-2-one
    Entries: EQ368301, EQ368302, EQ368303, EQ368304, EQ368305, EQ368306, EQ368307, EQ368308, EQ368309, EQ368351, EQ368352, EQ368353, EQ368354, EQ368355, EQ368356, EQ368357, EQ368358, EQ368359
    Molecule: Danofloxacin 1-cyclopropyl-6-fluoro-7-(5-methyl-2-aza-5-azoniabicyclo[2.2.1]heptan-2-yl)-4-oxoquinoline-3-carboxylate
    Entries: EQ366401, EQ366402, EQ366403, EQ366404, EQ366405, EQ366406, EQ366407, EQ366408, EQ366409
    Molecule: Ethoxyquin 6-ethoxy-2,2,4-trimethyl-1H-quinoline
    Entries: EQ365201, EQ365202, EQ365203, EQ365204, EQ365205, EQ365206, EQ365207, EQ365208, EQ365209
    Molecule: Benzenesulfonamide Benzenesulfonimidic acid
    Entries: EQ359501, EQ359502, EQ359503, EQ359504, EQ359505, EQ359506, EQ359507, EQ359508, EQ359509, EQ359551, EQ359552, EQ359553, EQ359554, EQ359555, EQ359556, EQ359557, EQ359558, EQ359559
    Molecule: Doxepine Doxepin 3-(6H-benzo[c][1]benzoxepin-11-ylidene)-N,N-dimethylpropan-1-amine
    Entries: EQ367601, EQ367602, EQ367603, EQ367604, EQ367605, EQ367606, EQ367607, EQ367608, EQ367609
    Molecule: Albendazole Methyl N-(6-propylsulfanyl-1H-benzimidazol-2-yl)carbamate
    Entries: EQ358001, EQ358002, EQ358003, EQ358004, EQ358005, EQ358006, EQ358007, EQ358008, EQ358009, EQ358051, EQ358052, EQ358053, EQ358054, EQ358055, EQ358056, EQ358057, EQ358058, EQ358059
    Molecule: delta9-Tetrahydrocannabinol delta9-THC 6,6,9-trimethyl-3-pentyl-6a,7,8,10a-tetrahydrobenzo[c]chromen-1-ol
    Entries: EQ362901, EQ362902, EQ362903, EQ362904, EQ362905, EQ362906, EQ362907, EQ362908, EQ362909
    Molecule: Florfenicol 2,2-dichloro-N-{3-fluoro-1-hydroxy-1-[4-(methylsulfonyl)phenyl]propan-2-yl}acetamide 2,2-dichloro-N-[3-fluoro-1-hydroxy-1-(4-methylsulfonylphenyl)propan-2-yl]acetamide
    Entries: EQ364101, EQ364102, EQ364103, EQ364104, EQ364105, EQ364106, EQ364107, EQ364108, EQ364109, EQ364151, EQ364152, EQ364153, EQ364154, EQ364155, EQ364156, EQ364157, EQ364158, EQ364159
    Molecule: Isoxaflutole (5-cyclopropyl-1,2-oxazol-4-yl)-[2-methylsulfonyl-4-(trifluoromethyl)phenyl]methanone
    Entries: EQ348801, EQ348802, EQ348803, EQ348804, EQ348805, EQ348806, EQ348807, EQ348808, EQ348809, EQ348851, EQ348852, EQ348853, EQ348854, EQ348855, EQ348856, EQ348857, EQ348858, EQ348859
    Molecule: 3,4-Methylenedioxyamphetamine (MDA) 1-(3,4-Methylenedioxyphenyl)-2-aminopropane 1-(1,3-benzodioxol-5-yl)propan-2-amine
    Entries: EQ371501, EQ371502, EQ371503, EQ371504, EQ371505, EQ371506, EQ371507, EQ371508, EQ371509
    Molecule: Albendazole sulfone methyl N-(6-propylsulfonyl-1H-benzimidazol-2-yl)carbamate
    Entries: EQ364751, EQ364752, EQ364753, EQ364754, EQ364755, EQ364756, EQ364757, EQ364758, EQ364759
    Molecule: PCP Pentachlorophenol 2,3,4,5,6-pentachlorophenol
    Entries: EQ371851, EQ371852, EQ371853, EQ371854, EQ371855, EQ371856
    Molecule: LSD D-Lysergic acid N,N-diethylamide N,N-diethyl-7-methyl-6,6a,8,9-tetrahydro-4H-indolo[4,3-fg]quinoline-9-carboxamide
    Entries: EQ368801, EQ368802, EQ368803, EQ368804, EQ368805, EQ368806, EQ368807, EQ368808, EQ368809
    Molecule: Forchlorfenuron N-(2-chloro-4-pyridyl)-N`-phenylurea 1-(2-chloropyridin-4-yl)-3-phenylurea
    Entries: EQ360101, EQ360102, EQ360103, EQ360104, EQ360105, EQ360106, EQ360107, EQ360108, EQ360109, EQ360151, EQ360152, EQ360153, EQ360154, EQ360155, EQ360156, EQ360157, EQ360158, EQ360159
    Molecule: Sarafloxacin 6-fluoro-1-(4-fluorophenyl)-4-oxo-7-piperazin-1-ylquinoline-3-carboxylic acid
    Entries: EQ362001, EQ362002, EQ362003, EQ362004, EQ362005, EQ362006, EQ362007, EQ362008, EQ362009
    Molecule: Bromazepam 7-bromo-5-pyridin-2-yl-1,3-dihydro-1,4-benzodiazepin-2-one
    Entries: EQ362801, EQ362802, EQ362803, EQ362804, EQ362805, EQ362806, EQ362807, EQ362808, EQ362809, EQ362851, EQ362852, EQ362853, EQ362854, EQ362855, EQ362856, EQ362857, EQ362858, EQ362859
    Molecule: 3,4-Methylenedioxy-N-methylamphetamine (MDMA) N-Methyl-3,4-methylenedioxyamphetamine 1-(1,3-benzodioxol-5-yl)-N-methylpropan-2-amine
    Entries: EQ361301, EQ361302, EQ361303, EQ361304, EQ361305, EQ361306, EQ361307, EQ361308, EQ361309
    Molecule: Fentanyl N-phenyl-N-[1-(2-phenylethyl)piperidin-4-yl]propanamide
    Entries: EQ364001, EQ364002, EQ364003, EQ364004, EQ364005, EQ364006, EQ364007, EQ364008, EQ364009
    Molecule: Chlordiazepoxide Zetran 7-chloro-4-hydroxy-N-methyl-5-phenyl-3H-1,4-benzodiazepin-2-imine
    Entries: EQ361701, EQ361702, EQ361703, EQ361704, EQ361705, EQ361706, EQ361707, EQ361708, EQ361709, EQ361751, EQ361752, EQ361753, EQ361754, EQ361755, EQ361756, EQ361757, EQ361758
    Molecule: Rifaximin (7S,9E,11S,12R,13S,14R,15R,16R,17S,18S,19E,21Z)-2,15,17,36-Tetrahydroxy-11-methoxy-3,7,12,14,16,18,22,30-octamethyl-6,23-dioxo-8,37-dioxa-24,27,33-triazahexacyclo[23.10.1.14,7.05,35.026,34.027, 32]heptatriaconta-1(35),2,4,9,19,21,25(36),26(34),...
    Entries: EQ369351, EQ369352, EQ369353, EQ369354, EQ369355, EQ369356, EQ369357, EQ369358, EQ369359
    Molecule: Dicloxacillin 6-[[3-(2,6-dichlorophenyl)-5-methyl-1,2-oxazole-4-carbonyl]amino]-3,3-dimethyl-7-oxo-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid
    Entries: EQ366501, EQ366502, EQ366503, EQ366504, EQ366505, EQ366506, EQ366507, EQ366508, EQ366509, EQ366551, EQ366552, EQ366553, EQ366554, EQ366555, EQ366556, EQ366557, EQ366558, EQ366559
    Molecule: Oxadiazon 5-tert-butyl-3-(2,4-dichloro-5-propan-2-yloxyphenyl)-1,3,4-oxadiazol-2-one
    Entries: EQ370401, EQ370402, EQ370403, EQ370404, EQ370405, EQ370406, EQ370407, EQ370408, EQ370409
    Molecule: Fenamiphos-sulfone Fenamiphos sulphone Ethyl 3-methyl-4-(methylsulfonyl)phenyl isopropylphosphoramidate
    Entries: EQ370301, EQ370302, EQ370303, EQ370304, EQ370305, EQ370306, EQ370307, EQ370308, EQ370309
    Molecule: 8-Hydroxy Mirtazapine 2-Methyl-1,2,3,4,10,14b-hexahydropyrazino[2,1-a]pyrido[2,3-c][2]benzazepin-8-ol
    Entries: EQ371201, EQ371202, EQ371203, EQ371204, EQ371205, EQ371206, EQ371207, EQ371208, EQ371209
    Molecule: Flumequine 9-Fluoro-5-methyl-1-oxo-6,7-dihydro-1H,5H-pyrido[3,2,1-ij]quinoline-2-carboxylic acid
    Entries: EQ364201, EQ364202, EQ364203, EQ364204, EQ364205, EQ364206, EQ364207, EQ364208, EQ364209
    Molecule: Marbofloxacin 9-Fluoro-3-methyl-10-(4-methyl-1-piperazinyl)-7-oxo-2,3-dihydro-7H-[1,3,4]oxadiazino[6,5,4-ij]quinoline-6-carboxylic acid
    Entries: EQ368901, EQ368902, EQ368903, EQ368904, EQ368905, EQ368906, EQ368907, EQ368908, EQ368909
    Molecule: Picolinafen N-(4-fluorophenyl)-6-[3-(trifluoromethyl)phenoxy]pyridine-2-carboxamide
    Entries: EQ360301, EQ360302, EQ360303, EQ360304, EQ360305, EQ360306, EQ360307, EQ360308, EQ360309
    Molecule: 2-(3-Hydroxycyclohexyl)-5-(2-methyl-2-octanyl)phenol CP47.497 2-(3-hydroxycyclohexyl)-5-(2-methyloctan-2-yl)phenol
    Entries: EQ369901, EQ369902, EQ369903, EQ369904, EQ369905, EQ369906, EQ369907, EQ369908, EQ369909, EQ369951, EQ369952, EQ369953, EQ369954, EQ369955, EQ369956, EQ369957, EQ369958, EQ369959
    Molecule: Phenobarbital 5-ethyl-5-phenyl-1,3-diazinane-2,4,6-trione
    Entries: EQ362401, EQ362402, EQ362403, EQ362404, EQ362405, EQ362406, EQ362407, EQ362408, EQ362409, EQ362451, EQ362452, EQ362453, EQ362454, EQ362455
    Molecule: Nordiazepam Nordazepam (1-Demethyldiazepam) 7-chloro-5-phenyl-1,3-dihydro-1,4-benzodiazepin-2-one
    Entries: EQ360801, EQ360802, EQ360803, EQ360804, EQ360805, EQ360806, EQ360807, EQ360808, EQ360809, EQ360851, EQ360852, EQ360853, EQ360854, EQ360855, EQ360856, EQ360857, EQ360858
    Molecule: Diclazuril 2-(4-chlorophenyl)-2-[2,6-dichloro-4-(3,5-dioxo-1,2,4-triazin-2-yl)phenyl]acetonitrile
    Entries: EQ365851, EQ365852, EQ365853, EQ365854, EQ365855, EQ365856
    Molecule: 7-amino-flunitrazepam 7-aminoflunitrazepam 7-amino-5-(2-fluorophenyl)-1-methyl-3H-1,4-benzodiazepin-2-one
    Entries: EQ371901, EQ371902, EQ371903, EQ371904, EQ371905, EQ371906, EQ371907, EQ371908, EQ371909
    Molecule: 1-Benzylpiperazine N-Benzylpiperazine 1-(phenylmethyl)piperazine
    Entries: EQ282001, EQ282002, EQ282003, EQ282004, EQ282005, EQ282006, EQ282007, EQ282008, EQ282009
    Molecule: Cimaterol 2-amino-5-[1-hydroxy-2-(propan-2-ylamino)ethyl]benzonitrile
    Entries: EQ365701, EQ365702, EQ365703, EQ365704, EQ365705, EQ365706, EQ365707, EQ365708, EQ365709
    Molecule: 3,4-Methylenedioxy-N-ethylamphetamine (MDEA) 3,4-methylenedioxyethamphetamine 1-(1,3-benzodioxol-5-yl)-N-ethylpropan-2-amine
    Entries: EQ360601, EQ360602, EQ360603, EQ360604, EQ360605, EQ360606, EQ360607, EQ360608, EQ360609
    Molecule: Ecgonine-methyl-ester (EME) Ecgonine methyl ester Methyl 3-hydroxy-8-methyl-8-azabicyclo[3.2.1]octane-4-carboxylate
    Entries: EQ361901, EQ361902, EQ361903, EQ361904, EQ361905, EQ361906, EQ361907, EQ361908, EQ361909
    Molecule: Indapamide 4-chloro-N-(2-methyl-2,3-dihydroindol-1-yl)-3-sulfamoylbenzamide
    Entries: EQ368601, EQ368602, EQ368603, EQ368604, EQ368605, EQ368606, EQ368607, EQ368608, EQ368609, EQ368651, EQ368652, EQ368653, EQ368654, EQ368655, EQ368656, EQ368657, EQ368658, EQ368659
    Molecule: Doxycycline 8-carbamoyl-10-(dimethylazaniumyl)-4,6a,7,11-tetrahydroxy-12-methyl-6,9-dioxo-10a,11,11a,12-tetrahydro-10H-tetracen-5-olate
    Entries: EQ367801, EQ367802, EQ367803, EQ367804, EQ367805, EQ367806, EQ367807, EQ367808, EQ367809, EQ367851, EQ367852, EQ367853, EQ367854, EQ367855, EQ367856, EQ367857, EQ367858, EQ367859
    Molecule: Vigabatrin 4-aminohex-5-enoic acid 4-Amino-5-hexenoic acid
    Entries: EQ362601, EQ362602, EQ362603, EQ362604, EQ362605, EQ362606, EQ362607, EQ362608, EQ362609
    Molecule: Nigericin 2-{6-[(2-{5`-[6-hydroxy-6-(hydroxymethyl)-3,5-dimethyltetrahydro-2h-pyran-2-yl]-2,3`-dimethyloctahydro-2,2`-bifuran-5-yl}-9-methoxy-2,4,10-trimethyl-1,6-dioxaspiro[4.5]dec-7-yl)methyl]-3-methyltetrahydro-2h-pyran-2-yl}propanoic acid 2-[6-[[2-[5...
    Entries: EQ368251, EQ368252, EQ368253, EQ368254, EQ368255, EQ368256, EQ368257, EQ368258, EQ368259
    Molecule: Oxytetracycline (4S,4aR,5S,5aR,6S,12aS)-4-(Dimethylamino)-3,5,6,10,12,12a-hexahydroxy-6-methyl-1,11-dioxo-1,4,4a,5,5a,6,11,12a-octahydro-2-tetracenecarboxamide 2-carbamoyl-4-(dimethylazaniumyl)-5,6,10,11,12a-pentahydroxy-6-methyl-3,12-dioxo-4,4a,5,5a-tet...
    Entries: EQ361001, EQ361002, EQ361003, EQ361004, EQ361005, EQ361006, EQ361007, EQ361008, EQ361009, EQ361051, EQ361052, EQ361053, EQ361054, EQ361055, EQ361056, EQ361057, EQ361058, EQ361059
    Molecule: Oxolinic acid 5-ethyl-8-oxo-[1,3]dioxolo[4,5-g]quinoline-7-carboxylic acid
    Entries: EQ360901, EQ360902, EQ360903, EQ360904, EQ360905, EQ360906, EQ360907, EQ360908, EQ360909
    Molecule: Mirtazapine 2-Methyl-1,2,3,4,10,14b-hexahydropyrazino[2,1-a]pyrido[2,3-c][2]benzazepine
    Entries: EQ362201, EQ362202, EQ362203, EQ362204, EQ362205, EQ362206, EQ362207, EQ362208, EQ362209
    Molecule: 1-Chlorobenzotriazole 1-chloro-benzotriazole Chlorobenzotriazole
    Entries: EQ369401, EQ369402, EQ369403, EQ369404, EQ369405, EQ369406, EQ369407, EQ369408, EQ369409, EQ369451, EQ369452, EQ369453, EQ369454, EQ369455, EQ369456, EQ369457, EQ369458, EQ369459
    Molecule: Lacosamide N2-Acetyl-N-benzyl-O-methyl-D-serinamide 2-acetamido-N-benzyl-3-methoxypropanamide
    Entries: EQ362701, EQ362702, EQ362703, EQ362704, EQ362705, EQ362706, EQ362707, EQ362708, EQ362709, EQ362751, EQ362752, EQ362753, EQ362754, EQ362755, EQ362756, EQ362757, EQ362758, EQ362759
    Molecule: Lincomycin N-[2-hydroxy-1-(3,4,5-trihydroxy-6-methylsulfanyloxan-2-yl)propyl]-1-methyl-4-propylpyrrolidine-2-carboxamide
    Entries: EQ368701, EQ368702, EQ368703, EQ368704, EQ368705, EQ368706, EQ368707, EQ368708, EQ368709, EQ368751, EQ368752, EQ368753, EQ368754, EQ368755, EQ368756, EQ368757, EQ368758, EQ368759
    Molecule: Norephedrine Phenylpropanolamine 2-amino-1-phenylpropan-1-ol
    Entries: EQ368401, EQ368402, EQ368403, EQ368404, EQ368405, EQ368406, EQ368407, EQ368408, EQ368409
    Molecule: Dichlorophen 4-chloro-2-[(5-chloro-2-hydroxyphenyl)methyl]phenol
    Entries: EQ364851, EQ364852, EQ364853, EQ364854, EQ364855, EQ364856, EQ364857, EQ364858, EQ364859
    Molecule: Clobazam 7-chloro-1-methyl-5-phenyl-1,5-benzodiazepine-2,4-dione
    Entries: EQ359301, EQ359302, EQ359303, EQ359304, EQ359305, EQ359306, EQ359307, EQ359308, EQ359309
    Molecule: Flunitrazepam 5-(2-fluorophenyl)-1-methyl-7-nitro-3H-1,4-benzodiazepin-2-one
    Entries: EQ368101, EQ368102, EQ368103, EQ368104, EQ368105, EQ368106, EQ368107, EQ368108, EQ368109
    Molecule: Chlortetracycline 2-carbamoyl-7-chloro-4-(dimethylazaniumyl)-6,10,11,12a-tetrahydroxy-6-methyl-3,12-dioxo-4,4a,5,5a-tetrahydrotetracen-1-olate
    Entries: EQ368051, EQ368052, EQ368053, EQ368054, EQ368055, EQ368056, EQ368057, EQ368058, EQ368059
    Molecule: Mebendazole methyl N-(6-benzoyl-1H-benzimidazol-2-yl)carbamate
    Entries: EQ366801, EQ366802, EQ366803, EQ366804, EQ366805, EQ366806, EQ366807, EQ366808, EQ366809, EQ366851, EQ366852, EQ366853, EQ366854, EQ366855, EQ366856, EQ366857, EQ366858, EQ366859
    Molecule: Flunixine Banamine 2-[2-methyl-3-(trifluoromethyl)anilino]pyridine-3-carboxylic acid
    Entries: EQ368501, EQ368502, EQ368503, EQ368504, EQ368505, EQ368506, EQ368507, EQ368508, EQ368509, EQ368551, EQ368552, EQ368553, EQ368554, EQ368555, EQ368556, EQ368557, EQ368558, EQ368559
    Molecule: Ronidazole (1-methyl-5-nitroimidazol-2-yl)methyl carbamate
    Entries: EQ369501, EQ369502, EQ369503, EQ369504, EQ369505, EQ369506, EQ369507, EQ369508, EQ369509
    Molecule: Bromuconazole 1-[[4-bromo-2-(2,4-dichlorophenyl)oxolan-2-yl]methyl]-1,2,4-triazole
    Entries: EQ370001, EQ370002, EQ370003, EQ370004, EQ370005, EQ370006, EQ370007, EQ370008, EQ370009
    Molecule: Paroxetine 3-[(1,3-Benzodioxol-5-yloxy)methyl]-4-(4-fluorophenyl)piperidine 3-(1,3-benzodioxol-5-yloxymethyl)-4-(4-fluorophenyl)piperidine
    Entries: EQ361101, EQ361102, EQ361103, EQ361104, EQ361105, EQ361106, EQ361107, EQ361108, EQ361109

Authors: Beck B, Stravs M, Schymanski E, Singer H, Department of Environmental Chemistry, Eawag   
CopyRight: Copyright (C) 2015 Eawag, Duebendorf, Switzerland   

    Molecule: Pretilachlor 2-chloranyl-N-(2,6-diethylphenyl)-N-(2-propoxyethyl)ethanamide
    Entries: EQ311701, EQ311702, EQ311703, EQ311704, EQ311705, EQ311706
    Molecule: Hexazinone 3-cyclohexyl-6-(dimethylamino)-1-methyl-1,3,5-triazine-2,4-dione
    Entries: EQ025601, EQ025602, EQ025603, EQ025604, EQ025605, EQ025606, EQ025607, EQ025608, EQ025609
    Molecule: Rivastigmine N-ethyl-N-methyl-carbamic acid [3-[(1S)-1-(dimethylamino)ethyl]phenyl] ester
    Entries: EQ284401, EQ284402, EQ284403, EQ284404, EQ284405, EQ284406
    Molecule: Hydrocortisone Cortisol (8S,9S,10R,11S,13S,14S,17R)-11,17-dihydroxy-17-(2-hydroxyacetyl)-10,13-dimethyl-2,6,7,8,9,11,12,14,15,16-decahydro-1H-cyclopenta[a]phenanthren-3-one
    Entries: EQ320101, EQ320102, EQ320103, EQ320104, EQ320105, EQ320106, EQ320107, EQ320108, EQ320109
    Molecule: Fenhexamid N-(2,3-dichloro-4-hydroxy-phenyl)-1-methyl-cyclohexanecarboxamide
    Entries: EQ305701, EQ305702, EQ305703, EQ305704, EQ305705, EQ305706, EQ305751, EQ305752, EQ305753, EQ305754, EQ305755, EQ305756
    Molecule: Norlidocaine Monoethylglycinexylidide N-(2,6-dimethylphenyl)-2-(ethylamino)acetamide
    Entries: EQ347101, EQ347102, EQ347103, EQ347104, EQ347105, EQ347106, EQ347107, EQ347108, EQ347109
    Molecule: Dinotefuran 2-methyl-1-nitro-3-(oxolan-3-ylmethyl)guanidine
    Entries: EQ310801, EQ310802, EQ310803, EQ310804, EQ310805, EQ310806, EQ310851, EQ310852, EQ310853, EQ310854, EQ310855, EQ310856
    Molecule: Flufenamic acid Arlef 2-[3-(trifluoromethyl)anilino]benzoic acid
    Entries: EQ302101, EQ302102, EQ302103, EQ302104, EQ302105, EQ302106, EQ302151, EQ302152, EQ302153, EQ302154, EQ302155, EQ302156
    Molecule: Prilocaine Propanamide, N-(2-methylphenyl)-2-(propylamino)- N-(2-methylphenyl)-2-(propylamino)propanamide
    Entries: EQ314101, EQ314102, EQ314103, EQ314104, EQ314105, EQ314106
    Molecule: Betamethasone (8S,9R,10S,11S,13S,14S,16S,17R)-9-fluoro-11,17-dihydroxy-17-(2-hydroxyacetyl)-10,13,16-trimethyl-6,7,8,11,12,14,15,16-octahydrocyclopenta[a]phenanthren-3-one
    Entries: EQ324201, EQ324205
    Molecule: 1-Methyl-1,2,3-benzotriazole 1-methylbenzotriazole
    Entries: EQ279801, EQ279802, EQ279803, EQ279804, EQ279805, EQ279806, EQ279807, EQ279808, EQ279809
    Molecule: Mepivacaine N-(2,6-dimethylphenyl)-1-methyl-2-piperidinecarboxamide
    Entries: EQ312601, EQ312602, EQ312603, EQ312604, EQ312605, EQ312606
    Molecule: N-Bisdesmethyl Tramadol 2-(aminomethyl)-1-(3-methoxyphenyl)cyclohexan-1-ol
    Entries: EQ339801, EQ339802, EQ339803, EQ339804, EQ339805, EQ339806, EQ339807, EQ339808, EQ339809
    Molecule: Genistein 5,7-dihydroxy-3-(4-hydroxyphenyl)chromen-4-one
    Entries: EQ326501, EQ326502, EQ326503, EQ326504, EQ326505, EQ326506, EQ326507, EQ326508, EQ326509, EQ326551, EQ326552, EQ326553, EQ326554, EQ326555, EQ326556, EQ326557, EQ326558, EQ326559
    Molecule: Trospium spiro[8-azoniabicyclo[3.2.1]octane-8,1`-azolidin-1-ium]-3-yl 2-hydroxy-2,2-diphenylacetate
    Entries: EQ303017, EQ303018, EQ303019, EQ303020, EQ303021, EQ303022
    Molecule: Nordeprenyl Desmethylselegiline 1-phenyl-N-prop-2-ynylpropan-2-amine
    Entries: EQ328101, EQ328102, EQ328103, EQ328104, EQ328105, EQ328106, EQ328107, EQ328108, EQ328109
    Molecule: Risperidone 3-[2-[4-(6-fluoro-1,2-benzoxazol-3-yl)piperidin-1-yl]ethyl]-2-methyl-6,7,8,9-tetrahydropyrido[1,2-a]pyrimidin-4-one
    Entries: EQ335301, EQ335302, EQ335303, EQ335304, EQ335305, EQ335306, EQ335307, EQ335308, EQ335309
    Molecule: Emtricitabine 4-amino-5-fluoro-1-[(2R,5S)-2-(hydroxymethyl)-1,3-oxathiolan-5-yl]-2-pyrimidinone
    Entries: EQ310601, EQ310602, EQ310603, EQ310604, EQ310605, EQ310606, EQ310651, EQ310652, EQ310653, EQ310654, EQ310655, EQ310656
    Molecule: Gemfibrozil 5-[(2,5-dimethylphenyl)oxy]-2,2-dimethylpentanoic acid 5-(2,5-dimethylphenoxy)-2,2-dimethyl-pentanoic acid
    Entries: EQ307101, EQ307102, EQ307103, EQ307104, EQ307105, EQ307106, EQ307151, EQ307152, EQ307153, EQ307154, EQ307155, EQ307156
    Molecule: Risperidone Pargyline N-oxide N-benzyl-N-methylprop-2-yn-1-amine oxide
    Entries: EQ327101, EQ327102, EQ327103, EQ327104, EQ327105, EQ327106, EQ327107, EQ327108, EQ327109
    Molecule: Glycyrrhetinic Acid 18-beta-Glycyrrhetin acid (2S,4aS,6aR,6aS,6bR,8aR,10S,12aS,14bR)-10-hydroxy-2,4a,6a,6b,9,9,12a-heptamethyl-13-oxo-3,4,5,6,6a,7,8,8a,10,11,12,14b-dodecahydro-1H-picene-2-carboxylic acid
    Entries: EQ326301, EQ326302, EQ326303, EQ326304, EQ326305, EQ326306, EQ326307, EQ326308, EQ326309, EQ326351, EQ326352, EQ326353, EQ326354, EQ326355, EQ326356, EQ326357, EQ326358, EQ326359
    Molecule: Coumafuryl Fumarin 3-[1-(furan-2-yl)-3-oxobutyl]-4-hydroxychromen-2-one
    Entries: EQ309101, EQ309102, EQ309103, EQ309104, EQ309105, EQ309106, EQ309151, EQ309152, EQ309153, EQ309154, EQ309155, EQ309156
    Molecule: Ofloxacin DL-8280
    Entries: EQ307301, EQ307302, EQ307303, EQ307304, EQ307305, EQ307306
    Molecule: Fonofos ethoxy-ethyl-(phenylthio)-sulfanylidenephosphorane
    Entries: EQ311201, EQ311202, EQ311203, EQ311204, EQ311205, EQ311206
    Molecule: 2-Aminoheptane 1-methylhexylamine 1-Methylhexylamine
    Entries: EQ315301, EQ315302, EQ315303, EQ315304, EQ315305, EQ315306
    Molecule: Diltiazem 1,5-Benzothiazepin-4(5H)-one, 3-(acetyloxy)-5-[2-(dimethylamino)ethyl]-2,3-dihydro-2-(4-methoxyphenyl)-, (2S-cis)- [5-(2-dimethylaminoethyl)-2-(4-methoxyphenyl)-4-oxidanylidene-2,3-dihydro-1,5-benzothiazepin-3-yl] ethanoate
    Entries: EQ301701, EQ301702, EQ301703, EQ301704, EQ301705, EQ301706
    Molecule: R-Deprenyl N-Oxide Methyl(1-phenyl-2-propanyl)2-propyn-1-ylamine oxide N-methyl-1-phenyl-N-prop-2-ynylpropan-2-amine oxide
    Entries: EQ346101, EQ346102, EQ346103, EQ346104, EQ346105, EQ346106, EQ346107, EQ346108, EQ346109
    Molecule: Dehydroepiandrosterone (DHEA) trans-Dehydroandosterone (3S,8R,9S,10R,13S,14S)-10,13-dimethyl-3-oxidanyl-1,2,3,4,7,8,9,11,12,14,15,16-dodecahydrocyclopenta[a]phenanthren-17-one
    Entries: EQ308501, EQ308502, EQ308503, EQ308504, EQ308505, EQ308506
    Molecule: Oxprenolol 1-(2-allyloxyphenoxy)-3-(isopropylamino)propan-2-ol
    Entries: EQ312801, EQ312802, EQ312803, EQ312804, EQ312805, EQ312806
    Molecule: Benalaxyl 2-(2,6-dimethyl-N-(1-oxo-2-phenylethyl)anilino)propanoic acid methyl ester
    Entries: EQ303501, EQ303502, EQ303503, EQ303504, EQ303505, EQ303506
    Molecule: 1,3-Dimethylpteridine-2,4-dione 1,3-Dimethyllumazine
    Entries: EQ310401, EQ310402, EQ310403, EQ310404, EQ310405, EQ310406
    Molecule: Fludrocortisone [2-[(8S,9R,10S,11S,13S,14S,17R)-9-fluoro-11,17-dihydroxy-10,13-dimethyl-3-oxo-1,2,6,7,8,11,12,14,15,16-decahydrocyclopenta[a]phenanthren-17-yl]-2-oxoethyl] acetate
    Entries: EQ324001, EQ324002, EQ324003, EQ324004, EQ324005, EQ324006, EQ324051, EQ324052, EQ324053, EQ324054, EQ324055
    Molecule: Coumachlor Ratilan 3-[1-(4-chlorophenyl)-3-oxobutyl]-4-hydroxychromen-2-one
    Entries: EQ309001, EQ309002, EQ309003, EQ309004, EQ309005, EQ309006, EQ309051, EQ309052, EQ309053, EQ309054, EQ309055, EQ309056
    Molecule: Amantadine 1-adamantanamine
    Entries: EQ312401, EQ312402, EQ312403, EQ312404, EQ312405, EQ312406
    Molecule: Cycloxydim 2-[(1E)-N-Ethoxybutanimidoyl]-3-hydroxy-5-(tetrahydro-2H-thiopyran-3-yl)-2-cyclohexen-1-one
    Entries: EQ304501, EQ304502, EQ304503, EQ304504, EQ304505, EQ304506, EQ304551, EQ304552, EQ304553, EQ304554, EQ304555, EQ304556
    Molecule: Nitrofen 2,4-bis(chloranyl)-1-(4-nitrophenoxy)benzene
    Entries: EQ309801, EQ309802, EQ309803, EQ309804, EQ309805, EQ309806
    Molecule: Terbinafine Terbinex N,6,6-trimethyl-N-(naphthalen-1-ylmethyl)hept-2-en-4-yn-1-amine
    Entries: EQ358601, EQ358602, EQ358603, EQ358604, EQ358605, EQ358606, EQ358607, EQ358608, EQ358609
    Molecule: Valsartan acid 4-[2-(2H-1,2,3,4-tetrazol-5-yl)phenyl]benzoic acid
    Entries: EQ279501, EQ279502, EQ279503, EQ279504, EQ279505, EQ279506, EQ279507, EQ279508, EQ279509, EQ279551, EQ279552, EQ279553, EQ279554, EQ279555, EQ279556, EQ279557, EQ279558, EQ279559
    Molecule: Letrozole 4-[(4-cyanophenyl)-(1,2,4-triazol-1-yl)methyl]benzonitrile
    Entries: EQ358501, EQ358502, EQ358503, EQ358504, EQ358505, EQ358506, EQ358507, EQ358508, EQ358509, EQ358551, EQ358552, EQ358553, EQ358554, EQ358555, EQ358556, EQ358557, EQ358558, EQ358559
    Molecule: Mianserin-N-Oxide 1,2,3,4,10,14b-Hexahydro-2-methyldibenzo(c,f)pyrazino(1,2-a)azepine 2-oxide
    Entries: EQ327801, EQ327802, EQ327803, EQ327804, EQ327805, EQ327806, EQ327807, EQ327808, EQ327809, EQ327851, EQ327852, EQ327853, EQ327854, EQ327855, EQ327856, EQ327857, EQ327858, EQ327859
    Molecule: Primaquine (4-amino-1-methyl-butyl)-(6-methoxy-8-quinolyl)amine
    Entries: EQ300901, EQ300902, EQ300903, EQ300904, EQ300905, EQ300906
    Molecule: N-Desmethyltramadol N-demethyltramadol 1-(3-methoxyphenyl)-2-(methylaminomethyl)cyclohexan-1-ol
    Entries: EQ356501, EQ356502, EQ356503, EQ356504, EQ356505, EQ356506, EQ356507, EQ356508, EQ356509
    Molecule: Lorazepam 2H-1,4-Benzodiazepin-2-one, 7-chloro-5-(2-chlorophenyl)-1,3-dihydro-3-hydroxy- 7-chloranyl-5-(2-chlorophenyl)-3-oxidanyl-1,3-dihydro-1,4-benzodiazepin-2-one
    Entries: EQ313801, EQ313802, EQ313803, EQ313804, EQ313805, EQ313806, EQ313851, EQ313852, EQ313853, EQ313854, EQ313855, EQ313856
    Molecule: Propanil N-(3,4-dichlorophenyl)propanamide
    Entries: EQ305101, EQ305102, EQ305103, EQ305104, EQ305105, EQ305106, EQ305151, EQ305152, EQ305153, EQ305154, EQ305155, EQ305156
    Molecule: Torsemide SMR000466313 1-[4-(3-methylanilino)pyridin-3-yl]sulfonyl-3-propan-2-ylurea
    Entries: EQ314501, EQ314502, EQ314503, EQ314504, EQ314505, EQ314506, EQ314551, EQ314552, EQ314553, EQ314554, EQ314555, EQ314556
    Molecule: 2-Aminobenzothiazole benzo(d)thiazol-2-amine 1,3-benzothiazol-2-amine
    Entries: EQ337401, EQ337402, EQ337403, EQ337404, EQ337405, EQ337406, EQ337407, EQ337408, EQ337409
    Molecule: Tiapride N-(2-diethylaminoethyl)-2-methoxy-5-methylsulfonyl-benzamide
    Entries: EQ302801, EQ302802, EQ302803, EQ302804, EQ302805, EQ302806, EQ302851, EQ302852, EQ302853, EQ302854, EQ302855, EQ302856
    Molecule: Corticosterone 11b,21-Dihydroxyprogesterone (8S,9S,10R,11S,13S,14S,17S)-11-hydroxy-17-(2-hydroxyacetyl)-10,13-dimethyl-1,2,6,7,8,9,11,12,14,15,16,17-dodecahydrocyclopenta[a]phenanthren-3-one
    Entries: EQ319801, EQ319802, EQ319803, EQ319804, EQ319805, EQ319806, EQ319807, EQ319808, EQ319809
    Molecule: N,N-Dimethyldipropylenetriamine 3-(3-aminopropylamino)propyl-dimethyl-amine
    Entries: EQ300701, EQ300702, EQ300703, EQ300704, EQ300705, EQ300706
    Molecule: Butylisopropylamine N-propan-2-ylbutan-1-amine N-Isopropyl-1-butanamine
    Entries: EQ359101, EQ359102, EQ359103, EQ359104, EQ359105, EQ359106, EQ359107, EQ359108, EQ359109
    Molecule: 4-Chlorobenzophenone 4-Chlorobenzophenon (4-chlorophenyl)-phenylmethanone
    Entries: EQ338001, EQ338002, EQ338003, EQ338004, EQ338005, EQ338006, EQ338007, EQ338008, EQ338009
    Molecule: 1-[(4-Chlorophenyl)phenylmethyl]piperazine Norchlorcyclizine 1-(4-Chlorobenzhydryl) piperazine
    Entries: EQ338301, EQ338302, EQ338303, EQ338304, EQ338305, EQ338306, EQ338307, EQ338308, EQ338309
    Molecule: Celiprolol 3-[3-acetyl-4-[3-(tert-butylamino)-2-hydroxy-propoxy]phenyl]-1,1-diethyl-urea
    Entries: EQ301501, EQ301502, EQ301503, EQ301504, EQ301505, EQ301506, EQ301551, EQ301552, EQ301553, EQ301554, EQ301555, EQ301556
    Molecule: Dinoterb 2-tert-butyl-4,6-dinitro-phenol Phenol, 2-(1,1-dimethylethyl)-4,6-dinitro-
    Entries: EQ310951, EQ310952, EQ310953, EQ310954, EQ310955, EQ310956
    Molecule: N-(2,4-dimethylphenyl)formamide
    Entries: EQ033001, EQ033002, EQ033003, EQ033004, EQ033005, EQ033006, EQ033007, EQ033008, EQ033009
    Molecule: N-Nitrosopiperazine (NPAZ) 1-nitrosopiperazine
    Entries: EQ335601, EQ335602, EQ335603, EQ335604, EQ335605, EQ335606
    Molecule: Nilotinib 4-methyl-N-(3-(4-methylimidazol-1-yl)-5-(trifluoromethyl)phenyl)-3-((4-pyridin-3-ylpyrimidin-2-yl)amino)benzamide
    Entries: EQ358801, EQ358802, EQ358803, EQ358804, EQ358805, EQ358806, EQ358807, EQ358808, EQ358809, EQ358851, EQ358852, EQ358853, EQ358854, EQ358855, EQ358856, EQ358857, EQ358858, EQ358859
    Molecule: Metribuzin-desamino Deaminometribuzin 6-tert-butyl-3-(methylthio)-2H-1,2,4-triazin-5-one
    Entries: EQ009101, EQ009102, EQ009103, EQ009104, EQ009105, EQ009106, EQ009107, EQ009108, EQ009109, EQ009151, EQ009152, EQ009153, EQ009154, EQ009155, EQ009156, EQ009157, EQ009158, EQ009159
    Molecule: Cetirizine N-Oxide (R)-Cetirizine N-Oxide 2-[2-[4-[(4-chlorophenyl)-phenylmethyl]-1-oxidopiperazin-1-ium-1-yl]ethoxy]acetic acid
    Entries: EQ338101, EQ338102, EQ338103, EQ338104, EQ338105, EQ338106, EQ338107, EQ338108, EQ338109, EQ338151, EQ338152, EQ338153, EQ338154, EQ338155, EQ338156, EQ338157, EQ338158, EQ338159
    Molecule: Pargyline N-benzyl-N-methylprop-2-yn-1-amine N-Methyl-N-propargylbenzylamine
    Entries: EQ300401, EQ300402, EQ300403, EQ300404, EQ300405, EQ300406
    Molecule: N-Nitrosopyrrolidine (NPYR) 1-nitrosopyrrolidine
    Entries: EQ345001, EQ345002, EQ345003, EQ345004, EQ345005, EQ345006
    Molecule: Methamphetamine (2S)-N-methyl-1-phenylpropan-2-amine
    Entries: EQ282901, EQ282902, EQ282903, EQ282904, EQ282905, EQ282906
    Molecule: Chloroxynil 3,5-Dichloro-4-hydroxybenzonitrile 3,5-bis(chloranyl)-4-oxidanyl-benzenecarbonitrile
    Entries: EQ304251, EQ304252, EQ304253, EQ304254, EQ304255, EQ304256
    Molecule: Pencycuron 1-(4-chlorobenzyl)-1-cyclopentyl-3-phenyl-urea
    Entries: EQ306401, EQ306402, EQ306403, EQ306404, EQ306405, EQ306406, EQ306451, EQ306452, EQ306453, EQ306454, EQ306455, EQ306456
    Molecule: Oxacillin (2S,5R,6R)-3,3-dimethyl-6-[(5-methyl-3-phenyl-1,2-oxazole-4-carbonyl)amino]-7-oxo-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid
    Entries: EQ320701, EQ320702, EQ320703, EQ320704, EQ320705, EQ320706, EQ320707, EQ320708, EQ320709, EQ320751, EQ320752, EQ320753, EQ320754, EQ320755, EQ320756, EQ320757, EQ320758, EQ320759
    Molecule: Pregabalin (3S)-3-(aminomethyl)-5-methyl-hexanoic acid
    Entries: EQ312901, EQ312902, EQ312903, EQ312904, EQ312905, EQ312906, EQ312951, EQ312952, EQ312953, EQ312954, EQ312955
    Molecule: 4-chloro-N-methylaniline
    Entries: EQ347801, EQ347802, EQ347803, EQ347804, EQ347805, EQ347806, EQ347807, EQ347808, EQ347809
    Molecule: Mexiletine 1-(2,6-dimethylphenoxy)-2-propanamine
    Entries: EQ301001, EQ301002, EQ301003, EQ301004, EQ301005, EQ301006
    Molecule: Azobenzol Fentoxan [(Z)-Phenyl-NNO-azoxy]benzene
    Entries: EQ308801, EQ308802, EQ308803, EQ308804, EQ308805, EQ308806
    Molecule: Betamethasone 21 [2-[(8S,9R,10S,11S,13S,14S,16S,17R)-9-fluoro-11,17-dihydroxy-10,13,16-trimethyl-3-oxo-6,7,8,11,12,14,15,16-octahydrocyclopenta[a]phenanthren-17-yl]-2-oxoethyl] acetate
    Entries: EQ324401, EQ324402, EQ324404, EQ324405, EQ324406, EQ324451, EQ324452, EQ324453, EQ324454, EQ324455, EQ324456
    Molecule: Tenofovir [(1R)-2-adenin-9-yl-1-methyl-ethoxy]methylphosphonic acid
    Entries: EQ310501, EQ310502, EQ310503, EQ310504, EQ310505, EQ310506, EQ310551, EQ310552, EQ310553, EQ310554, EQ310555, EQ310556
    Molecule: Minocycline Minocin (4S,4aS,5aR,12aS)-4,7-Bis(dimethylamino)-3,10,12,12a-tetrahydroxy-1,11-dioxo-1,4,4a,5,5a,6,11,12a-octahydro-2-tetracenecarboxamide
    Entries: EQ320501, EQ320502, EQ320503, EQ320504, EQ320505, EQ320506, EQ320507, EQ320508, EQ320509, EQ320551, EQ320552, EQ320553, EQ320554, EQ320555, EQ320556, EQ320557, EQ320558, EQ320559
    Molecule: 4-chloro-N,N-dimethylaniline 4-Chloro-N,N-dimethylaniline
    Entries: EQ347701, EQ347702, EQ347703, EQ347704, EQ347705, EQ347706, EQ347707, EQ347708, EQ347709
    Molecule: Tributylphosphine oxide Tri-n-butylphosphine oxide
    Entries: EQ358901, EQ358902, EQ358903, EQ358904, EQ358905, EQ358906, EQ358907, EQ358908, EQ358909
    Molecule: N-Nitrosodibutylamine (NDBA) N,N-dibutylnitrous amide
    Entries: EQ346001, EQ346002, EQ346003, EQ346004, EQ346005, EQ346006
    Molecule: N-Nitrosomorpholine (NMOR) 4-nitrosomorpholine
    Entries: EQ345401, EQ345402, EQ345403, EQ345404, EQ345405, EQ345406
    Molecule: ortho-Chlorophenylpiperazine 1-(2-Chlorophenyl)piperazine 1-(2-chlorophenyl)piperazine
    Entries: EQ347001, EQ347002, EQ347003, EQ347004, EQ347005, EQ347006, EQ347007, EQ347008, EQ347009
    Molecule: N-desmethylpheniramine N-methyl-3-phenyl-3-pyridin-2-ylpropan-1-amine
    Entries: EQ342001, EQ342002, EQ342003, EQ342004, EQ342005, EQ342006, EQ342007, EQ342008, EQ342009
    Molecule: Sulpiride Sulpirid N-[(1-ethyl-2-pyrrolidinyl)methyl]-2-methoxy-5-sulfamoylbenzamide
    Entries: EQ302701, EQ302702, EQ302703, EQ302704, EQ302705, EQ302706, EQ302751, EQ302752, EQ302753, EQ302754, EQ302755, EQ302756
    Molecule: Nafcillin (2S,5R,6R)-6-[(2-ethoxynaphthalene-1-carbonyl)amino]-3,3-dimethyl-7-oxo-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid
    Entries: EQ320601, EQ320602, EQ320603, EQ320604, EQ320605, EQ320606, EQ320607, EQ320608, EQ320609, EQ320651, EQ320652, EQ320653, EQ320654, EQ320655, EQ320656, EQ320657, EQ320658, EQ320659
    Molecule: Dichlofluanid 1,1-dichloro-N-(dimethylsulfamoyl)-1-fluoro-N-phenylmethanesulfenamide
    Entries: EQ029701, EQ029702, EQ029703, EQ029704, EQ029705, EQ029706
    Molecule: 8phiC8SPC Octacarboxy sulfophenyl carboxylic acid 8-(4-sulfophenyl)octanoic acid
    Entries: EQ357101, EQ357102, EQ357103, EQ357104, EQ357105, EQ357106, EQ357107, EQ357108, EQ357109, EQ357151, EQ357152, EQ357153, EQ357154, EQ357155, EQ357156, EQ357157, EQ357158, EQ357159
    Molecule: Vildagliptin (2S)-1-[2-[(3-hydroxy-1-adamantyl)amino]-1-oxoethyl]-2-pyrrolidinecarbonitrile Equa
    Entries: EQ314601, EQ314602, EQ314603, EQ314604, EQ314605, EQ314606
    Molecule: Levofloxacin (3S)-9-Fluoro-3-methyl-10-(4-methyl-1-piperazinyl)-7-oxo-2,3-dihydro-7H-[1,4]oxazino[2,3,4-ij]quinoline-6-carboxylic acid
    Entries: EQ323901, EQ323902, EQ323904, EQ323905, EQ323906
    Molecule: Imiprothrin (2,5-dioxo-3-prop-2-ynylimidazolidin-1-yl)methyl 2,2-dimethyl-3-(2-methylprop-1-enyl)cyclopropane-1-carboxylate
    Entries: EQ293201, EQ293202, EQ293203, EQ293204, EQ293205, EQ293206
    Molecule: Crotamiton N-ethyl-N-(2-methylphenyl)but-2-enamide
    Entries: EQ320001, EQ320002, EQ320003, EQ320004, EQ320005, EQ320006, EQ320007, EQ320008, EQ320009
    Molecule: Alfuzosin N-{3-[(4-Amino-6,7-dimethoxy-2-quinazolinyl)(methyl)amino]propyl}tetrahydro-2-furancarboxamide N-[3-[(4-amino-6,7-dimethoxyquinazolin-2-yl)-methylamino]propyl]oxolane-2-carboxamide
    Entries: EQ335001, EQ335002, EQ335003, EQ335004, EQ335005, EQ335006, EQ335007, EQ335008, EQ335009, EQ335051, EQ335052, EQ335053, EQ335054, EQ335055, EQ335056, EQ335057, EQ335058, EQ335059
    Molecule: Flecainide N-(2-piperidinylmethyl)-2,5-bis(2,2,2-trifluoroethoxy)benzamide
    Entries: EQ302001, EQ302002, EQ302003, EQ302004, EQ302005, EQ302006, EQ302051, EQ302052, EQ302053, EQ302054, EQ302055, EQ302056
    Molecule: Amisulpride N-Oxide 4-amino-N-[(1-ethyl-1-oxidopyrrolidin-1-ium-2-yl)methyl]-5-ethylsulfonyl-2-methoxybenzamide
    Entries: EQ339901, EQ339902, EQ339903, EQ339904, EQ339905, EQ339906, EQ339907, EQ339908, EQ339909, EQ339951, EQ339952, EQ339953, EQ339954, EQ339955, EQ339956, EQ339957, EQ339958, EQ339959
    Molecule: Triclabendazole Triclabendazol 6-chloro-5-(2,3-dichlorophenoxy)-2-methylsulfanyl-1H-benzimidazole
    Entries: EQ321101, EQ321102, EQ321103, EQ321104, EQ321105, EQ321106, EQ321107, EQ321108, EQ321109, EQ321151, EQ321152, EQ321153, EQ321154, EQ321155, EQ321156, EQ321157, EQ321158, EQ321159
    Molecule: Prednisone (8S,9S,10R,13S,14S,17R)-17-hydroxy-17-(2-hydroxyacetyl)-10,13-dimethyl-6,7,8,9,12,14,15,16-octahydrocyclopenta[a]phenanthrene-3,11-dione
    Entries: EQ324301, EQ324305, EQ324351, EQ324355, EQ324356
    Molecule: Cortisone 17,21-Dihydroxypregn-4-ene-3,11,20-trione (8S,9S,10R,13S,14S,17R)-17-hydroxy-17-(2-hydroxyacetyl)-10,13-dimethyl-1,2,6,7,8,9,12,14,15,16-decahydrocyclopenta[a]phenanthrene-3,11-dione
    Entries: EQ319901, EQ319902, EQ319903, EQ319904, EQ319905, EQ319906, EQ319907, EQ319908, EQ319909, EQ319951, EQ319952, EQ319953, EQ319954, EQ319955, EQ319956, EQ319957, EQ319958, EQ319959
    Molecule: CGA62826 (2-[2,6-dimethylphenyl)-methoxyacetylamino]propionic acid N-(2,6-Dimethylphenyl)-N-(methoxyacetyl)-L-alanine (2S)-2-(N-(2-methoxyacetyl)-2,6-dimethylanilino)propanoic acid
    Entries: EQ323301, EQ323302, EQ323303, EQ323304, EQ323305, EQ323306, EQ323307, EQ323308, EQ323309, EQ323351, EQ323352, EQ323353, EQ323354, EQ323355, EQ323356, EQ323357, EQ323358, EQ323359
    Molecule: 2-Imidazolidinethione Ethylenethiourea
    Entries: EQ027101, EQ027102, EQ027103, EQ027104, EQ027105, EQ027106, EQ027107, EQ027108, EQ027109
    Molecule: 2-Mercaptobenzothiazole 1,3-Benzothiazole-2(3H)-thione 3H-1,3-benzothiazole-2-thione
    Entries: EQ319001, EQ319002, EQ319003, EQ319004, EQ319005, EQ319006, EQ319007, EQ319008, EQ319009, EQ319051, EQ319052, EQ319053, EQ319054, EQ319055, EQ319056, EQ319057, EQ319058, EQ319059
    Molecule: Medazepam 7-chloranyl-1-methyl-5-phenyl-2,3-dihydro-1,4-benzodiazepine
    Entries: EQ302201, EQ302202, EQ302203, EQ302204, EQ302205, EQ302206
    Molecule: 1-(Methoxymethyl)-1H-benzotriazole 1-(methoxymethyl)benzotriazole
    Entries: EQ318301, EQ318302, EQ318303, EQ318304, EQ318305, EQ318306, EQ318307, EQ318308, EQ318309
    Molecule: Didanosine 9-[(2R,5S)-5-(hydroxymethyl)-2-oxolanyl]-3H-purin-6-one
    Entries: EQ313501, EQ313502, EQ313503, EQ313504, EQ313505, EQ313506, EQ313551, EQ313552, EQ313553, EQ313554, EQ313555, EQ313556
    Molecule: Deferasirox
    Entries: EQ284201, EQ284202, EQ284203, EQ284205, EQ284206, EQ284251, EQ284252, EQ284253, EQ284254, EQ284255, EQ284256
    Molecule: Chlorthalidone Chlortalidone 2-chloranyl-5-(1-oxidanyl-3-oxidanylidene-2H-isoindol-1-yl)benzenesulfonamide
    Entries: EQ313401, EQ313402, EQ313403, EQ313404, EQ313405, EQ313406, EQ313451, EQ313452, EQ313453, EQ313454, EQ313455, EQ313456
    Molecule: Venlafaxine 1-[2-(dimethylamino)-1-(4-methoxyphenyl)ethyl]-1-cyclohexanol
    Entries: EQ064501, EQ064502, EQ064503, EQ064504, EQ064505, EQ064506
    Molecule: Procymidone 3-(3,5-dichlorophenyl)-1,5-dimethyl-3-azabicyclo[3.1.0]hexane-2,4-dione
    Entries: EQ310201, EQ310202, EQ310203, EQ310204, EQ310205, EQ310206
    Molecule: Venlafaxine N-Oxide 2-(1-hydroxycyclohexyl)-2-(4-methoxyphenyl)-N,N-dimethylethanamine oxide
    Entries: EQ339601, EQ339602, EQ339603, EQ339604, EQ339605, EQ339606, EQ339607, EQ339608, EQ339609
    Molecule: Chlorthiazide 6-chloro-1,1-dioxo-4H-1$l^{6},2,4-benzothiadiazine-7-sulfonamide
    Entries: EQ325601, EQ325602, EQ325603, EQ325604, EQ325605, EQ325606, EQ325607, EQ325608, EQ325609, EQ325651, EQ325652, EQ325653, EQ325654, EQ325655, EQ325656, EQ325657, EQ325658, EQ325659
    Molecule: Oxybutynin 2-cyclohexyl-2-hydroxy-2-phenyl-acetic acid 4-(diethylamino)but-2-ynyl ester
    Entries: EQ302501, EQ302502, EQ302503, EQ302504, EQ302505, EQ302506
    Molecule: Trimipramine N-Oxide Trimipramine N-oxide, (+/-)- 3-(5,6-dihydrobenzo[b][1]benzazepin-11-yl)-N,N,2-trimethylpropan-1-amine oxide
    Entries: EQ356101, EQ356102, EQ356103, EQ356104, EQ356105, EQ356106, EQ356107, EQ356108, EQ356109
    Molecule: Ametryn 4-N-ethyl-6-methylsulfanyl-2-N-propan-2-yl-1,3,5-triazine-2,4-diamine
    Entries: EQ303401, EQ303402, EQ303403, EQ303404, EQ303405, EQ303406
    Molecule: Bisoprolol 1-(propan-2-ylamino)-3-[4-(2-propan-2-yloxyethoxymethyl)phenoxy]-2-propanol
    Entries: EQ301301, EQ301302, EQ301303, EQ301304, EQ301305, EQ301306
    Molecule: Rufinamide 1-(2,6-difluorobenzyl)triazole-4-carboxamide
    Entries: EQ314301, EQ314302, EQ314303, EQ314304, EQ314305, EQ314306
    Molecule: N-tert-Butylisopropylamine 2-methyl-N-propan-2-ylpropan-2-amine N-Isopropyl-2-methyl-2-propanamine
    Entries: EQ363101, EQ363102, EQ363103, EQ363104, EQ363105, EQ363106, EQ363107, EQ363108, EQ363109
    Molecule: Fluazifop-butyl 2-[4-[[5-(trifluoromethyl)-2-pyridinyl]oxy]phenoxy]propanoic acid butyl ester
    Entries: EQ309301, EQ309302, EQ309303, EQ309304, EQ309305, EQ309306
    Molecule: Imidacloprid-guanidine Desnitro-imidacloprid 1-[(6-Chloro-3-pyridinyl)methyl]-4,5-dihydro-1H-imidazol-2-amine
    Entries: EQ305901, EQ305902, EQ305903, EQ305904, EQ305905, EQ305906
    Molecule: Norfenfluramine 1-[3-(trifluoromethyl)phenyl]propan-2-amine
    Entries: EQ346801, EQ346802, EQ346803, EQ346804, EQ346805, EQ346806, EQ346807, EQ346808, EQ346809
    Molecule: Cotinine (5S)-1-methyl-5-pyridin-3-ylpyrrolidin-2-one
    Entries: EQ328201, EQ328202, EQ328203, EQ328204, EQ328205, EQ328206, EQ328207, EQ328208, EQ328209
    Molecule: Atazanavir N-[(1S)-1-[[(1S,2S)-1-benzyl-3-[[[(2S)-2-(carbomethoxyamino)-3,3-dimethyl-butanoyl]amino]-[4-(2-pyridyl)benzyl]amino]-2-hydroxy-propyl]carbamoyl]-2,2-dimethyl-propyl]carbamic acid methyl ester
    Entries: EQ301201, EQ301202, EQ301203, EQ301204, EQ301205, EQ301206, EQ301251, EQ301252, EQ301253, EQ301254, EQ301255, EQ301256
    Molecule: Pheniramine-N-oxide N,N-dimethyl-3-phenyl-3-pyridin-2-ylpropan-1-amine oxide
    Entries: EQ325401, EQ325402, EQ325403, EQ325404, EQ325405, EQ325406, EQ325407, EQ325408, EQ325409
    Molecule: 5-Hydroxy Diclofenac 5-hydroxydiclofenac 2-[2-(2,6-dichloroanilino)-5-hydroxyphenyl]acetic acid
    Entries: EQ340301, EQ340302, EQ340303, EQ340304, EQ340305, EQ340306, EQ340307, EQ340308, EQ340309, EQ340351, EQ340352, EQ340353, EQ340354, EQ340355, EQ340356, EQ340357, EQ340358, EQ340359
    Molecule: 2-Hydroxybenzothiazole 2-Benzothiazolol 3H-1,3-benzothiazol-2-one
    Entries: EQ337501, EQ337502, EQ337503, EQ337504, EQ337505, EQ337506, EQ337507, EQ337508, EQ337509, EQ337551, EQ337552, EQ337553, EQ337554, EQ337555, EQ337556, EQ337557, EQ337558, EQ337559
    Molecule: Naftifine Suadian N-methyl-N-(naphthalen-1-ylmethyl)-3-phenylprop-2-en-1-amine
    Entries: EQ358101, EQ358102, EQ358103, EQ358104, EQ358105, EQ358106, EQ358107, EQ358108, EQ358109
    Molecule: 10phiC10SPC Decacarboxy sulfophenyl carboxylic acid 10-(4-sulfophenyl)octanoic acid
    Entries: EQ357001, EQ357002, EQ357003, EQ357004, EQ357005, EQ357006, EQ357007, EQ357008, EQ357009, EQ357051, EQ357052, EQ357053, EQ357054, EQ357055, EQ357056, EQ357057, EQ357058, EQ357059
    Molecule: 1H-Benzotriazole-5-carboxylic acid 1H-benzotriazole-5-carboxylic acid 2H-benzotriazole-5-carboxylic acid
    Entries: EQ318401, EQ318402, EQ318403, EQ318404, EQ318405, EQ318406, EQ318407, EQ318408, EQ318409, EQ318451, EQ318452, EQ318453, EQ318454, EQ318455, EQ318456, EQ318457, EQ318458, EQ318459
    Molecule: Carbaryl Carbaril N-methylcarbamic acid 1-naphthalenyl ester
    Entries: EQ303801, EQ303802, EQ303803, EQ303804, EQ303805, EQ303806, EQ303851, EQ303852, EQ303853, EQ303854, EQ303855, EQ303856
    Molecule: N-Nitrosodimethylamine (NDMA) N,N-dimethylnitrous amide
    Entries: EQ344701, EQ344702, EQ344703, EQ344704, EQ344705, EQ344706
    Molecule: Nadolol NADOLOL (2R,3S)-5-[3-(tert-butylamino)-2-hydroxy-propoxy]tetralin-2,3-diol
    Entries: EQ307201, EQ307202, EQ307203, EQ307204, EQ307205, EQ307206
    Molecule: Timolol 1-(tert-butylamino)-3-[(4-morpholin-4-yl-1,2,5-thiadiazol-3-yl)oxy]propan-2-ol
    Entries: EQ321001, EQ321002, EQ321003, EQ321004, EQ321005, EQ321006, EQ321007, EQ321008, EQ321009
    Molecule: 4-Hydroxybenzotriazole
    Entries: EQ290001, EQ290002, EQ290003, EQ290004, EQ290005, EQ290006, EQ290007, EQ290008, EQ290009, EQ290051, EQ290052, EQ290053, EQ290054, EQ290055, EQ290056, EQ290057, EQ290058, EQ290059
    Molecule: Triflumuron 2-chloranyl-N-[[4-(trifluoromethyloxy)phenyl]carbamoyl]benzamide
    Entries: EQ311801, EQ311802, EQ311803, EQ311804, EQ311805, EQ311806, EQ311851, EQ311852, EQ311853, EQ311854, EQ311855, EQ311856
    Molecule: Normianserin N-desmethylmianserin 1,2,3,4,10,14b-Hexahydrodibenzo[c,f]pyrazino[1,2-a]azepine
    Entries: EQ327701, EQ327702, EQ327703, EQ327704, EQ327705, EQ327706, EQ327707, EQ327708, EQ327709
    Molecule: Mianserin 2-Methyl-1,2,3,4,10,14b-hexahydrodibenzo[c,f]pyrazino[1,2-a]azepine
    Entries: EQ327401, EQ327402, EQ327403, EQ327404, EQ327405, EQ327406, EQ327407, EQ327408, EQ327409
    Molecule: N-Nitrosomethylethylamine (NMEA) N-ethyl-N-methylnitrous amide
    Entries: EQ344901, EQ344902, EQ344903, EQ344904, EQ344905, EQ344906
    Molecule: Enoxacin 1-Ethyl-6-fluoro-1,4-dihydro-4-oxo-7-(1-piperazinyl)-1,8-naphthyridine-3-carboxylic acid 1-ethyl-6-fluoranyl-4-oxidanylidene-7-piperazin-1-yl-1,8-naphthyridine-3-carboxylic acid
    Entries: EQ307801, EQ307802, EQ307803, EQ307804, EQ307805, EQ307806
    Molecule: Ticlopidine 5-(2-chlorobenzyl)-6,7-dihydro-4H-thieno[3,2-c]pyridine
    Entries: EQ302901, EQ302902, EQ302903, EQ302904, EQ302905, EQ302906, EQ302951, EQ302952, EQ302953
    Molecule: Fosthiazate 3-[(butan-2-ylthio)-ethoxyphosphoryl]-2-thiazolidinone
    Entries: EQ305001, EQ305002, EQ305003, EQ305004, EQ305005
    Molecule: Vinclozolin 3-(3,5-dichlorophenyl)-5-ethenyl-5-methyl-1,3-oxazolidine-2,4-dione
    Entries: EQ311901, EQ311902, EQ311903, EQ311904, EQ311905, EQ311906
    Molecule: Ketoconazole 1-[4-[4-[[(2R,4S)-2-(2,4-dichlorophenyl)-2-(imidazol-1-ylmethyl)-1,3-dioxolan-4-yl]methoxy]phenyl]piperazin-1-yl]ethanone
    Entries: EQ326201, EQ326202, EQ326203, EQ326204, EQ326205, EQ326206, EQ326207, EQ326208, EQ326209
    Molecule: Varenicline 7,8,9,10-tetrahydro-6h-6,10-methanoazepino[4,5-g]quinoxaline 7,8,9,10-tetrahydro-6,10-Methano-6H-pyrazino[2,3-h][3]benzazepine
    Entries: EQ363601, EQ363602, EQ363603, EQ363604, EQ363605, EQ363606, EQ363607, EQ363608, EQ363609
    Molecule: 2-Phenylpiperidine-2-acetamide 2-phenyl-2-piperidin-2-ylacetamide 2-Phenyl-2-(2-piperidinyl)acetamide
    Entries: EQ358401, EQ358402, EQ358403, EQ358404, EQ358405, EQ358406, EQ358407, EQ358408, EQ358409
    Molecule: Heptenophos (6-chloranyl-7-bicyclo[3.2.0]hepta-3,6-dienyl) dimethyl phosphate Ragadan
    Entries: EQ311401, EQ311402, EQ311403, EQ311404, EQ311405, EQ311406, EQ311451, EQ311452, EQ311453, EQ311454, EQ311455, EQ311456
    Molecule: Capecitabine (1-(5-Deoxy-beta-D-ribofuranosyl)-5-fluoro-1,2-dihydro-2-oxo-4-pyrimidinyl)-carbamic acid pentyl ester N-[1-[(2R,3R,4S,5R)-3,4-dihydroxy-5-methyl-2-oxolanyl]-5-fluoro-2-oxo-4-pyrimidinyl]carbamic acid pentyl ester
    Entries: EQ284501, EQ284505, EQ284506, EQ284551, EQ284552, EQ284553, EQ284554, EQ284555, EQ284556
    Molecule: Pheniramine N,N-dimethyl-3-phenyl-3-pyridin-2-ylpropan-1-amine
    Entries: EQ300501, EQ300502, EQ300503, EQ300504, EQ300505, EQ300506, EQ300507, EQ300508, EQ300509
    Molecule: Zonisamide AD-810 1,2-benzoxazol-3-ylmethanesulfonamide
    Entries: EQ313201, EQ313202, EQ313203, EQ313204, EQ313205, EQ313206, EQ313251, EQ313252, EQ313253, EQ313254, EQ313255, EQ313256
    Molecule: Fluvastatin (3R,5S)-7-[3-(4-fluorophenyl)-1-isopropyl-2-indolyl]-3,5-dihydroxy-6-heptenoic acid
    Entries: EQ313601, EQ313602, EQ313603, EQ313604, EQ313605, EQ313606, EQ313651, EQ313652, EQ313653, EQ313654, EQ313655, EQ313656
    Molecule: Fluoxastrobin (E)-[[2-[6-(2-chlorophenoxy)-5-fluoro-pyrimidin-4-yl]oxyphenyl]-(5,6-dihydro-1,4,2-dioxazin-3-yl)methylene]-methoxy-amine
    Entries: EQ305801, EQ305802, EQ305803, EQ305804, EQ305805, EQ305806
    Molecule: 1-Hydroxybenzotriazole
    Entries: EQ289901, EQ289902, EQ289903, EQ289904, EQ289905, EQ289906, EQ289907, EQ289908, EQ289909, EQ289951, EQ289952, EQ289953, EQ289954, EQ289955, EQ289956, EQ289957, EQ289958, EQ289959
    Molecule: Naptalam 2-(1-naphthylcarbamoyl)benzoic acid
    Entries: EQ309701, EQ309702, EQ309703, EQ309704, EQ309705, EQ309706, EQ309751, EQ309752, EQ309753, EQ309754, EQ309755, EQ309756
    Molecule: Memantine 3,5-dimethyladamantan-1-amine
    Entries: EQ335101, EQ335102, EQ335103, EQ335104, EQ335105, EQ335106, EQ335107, EQ335108, EQ335109
    Molecule: Penconazole 1-[2-(2,4-dichlorophenyl)pentyl]-1,2,4-triazole
    Entries: EQ310701, EQ310702, EQ310703, EQ310704, EQ310705, EQ310706, EQ310707, EQ310708, EQ310709
    Molecule: Flucytosine (5-FC) 6-Amino-5-fluoro-2(1H)-pyrimidinone 6-amino-5-fluoro-1H-pyrimidin-2-one
    Entries: EQ312101, EQ312102, EQ312103, EQ312104, EQ312105, EQ312106, EQ312151, EQ312152, EQ312153, EQ312154, EQ312155, EQ312156
    Molecule: Benthiavalicarb-isopropyl Isopropyl (1-{[1-(6-fluoro-1,3-benzothiazol-2-yl)ethyl]amino}-3-methyl-1-oxo-2-butanyl)carbamate propan-2-yl N-[(2S)-1-[[(1R)-1-(6-fluoro-1, 3-benzothiazol-2-yl)ethyl]amino]-3-methyl-1-oxobutan-2-yl]carbamate
    Entries: EQ305501, EQ305502, EQ305503, EQ305504, EQ305505, EQ305506, EQ305551, EQ305552, EQ305553, EQ305554, EQ305555, EQ305556
    Molecule: Tributylamine 1-Butanamine, N,N-dibutyl- N,N-dibutyl-1-butanamine
    Entries: EQ315001, EQ315002, EQ315003, EQ315004, EQ315005, EQ315006
    Molecule: Nornicotine 3-pyrrolidin-2-ylpyridine
    Entries: EQ328001, EQ328002, EQ328003, EQ328004, EQ328005, EQ328006, EQ328007, EQ328008, EQ328009
    Molecule: Dexamethasone acetate Dexamethason-21-acetate [2-(9-fluoro-11,17-dihydroxy-10,13,16-trimethyl-3-oxo-6,7,8,11,12,14,15,16-octahydrocyclopenta[a]phenanthren-17-yl)-2-oxoethyl] acetate
    Entries: EQ326401, EQ326402, EQ326403, EQ326404, EQ326405, EQ326406, EQ326407, EQ326408, EQ326409, EQ326451, EQ326452, EQ326453, EQ326454, EQ326455, EQ326456, EQ326457, EQ326458, EQ326459
    Molecule: Parathion diethoxy-(4-nitrophenoxy)-sulfanylidene-$l^{5}-phosphane
    Entries: EQ311501, EQ311502, EQ311503, EQ311504, EQ311505, EQ311506
    Molecule: 1H-Benzotriazole, 4(or 5)-methyl- 4-methyl-2H-benzotriazole
    Entries: EQ279101, EQ279102, EQ279103, EQ279104, EQ279105, EQ279106, EQ279107, EQ279108, EQ279109, EQ279151, EQ279152, EQ279153, EQ279154, EQ279155, EQ279156, EQ279157, EQ279158, EQ279159
    Molecule: Melamine 1,3,5-Triazine-2,4,6-triamine 1,3,5-triazine-2,4,6-triamine
    Entries: EQ315101, EQ315102, EQ315103, EQ315104, EQ315105, EQ315106
    Molecule: NOA407475 3-[(2-chloro-1,3-thiazol-5-yl)methyl]-5-methyl-1,3,5-oxadiazinan-4-imine 3-(2-chloro-thiazol-5-ylmethyl)-5-methyl-[1,3,5]oxadiazinan-4-ylideneamine
    Entries: EQ323401, EQ323402, EQ323403, EQ323404, EQ323405, EQ323406, EQ323407, EQ323408, EQ323409
    Molecule: Propazine-2-hydroxy 1,3,5-Triazin-2(1H)-one, 4,6-bis((1-methylethyl)amino)- 2,6-bis(isopropylamino)-1H-s-triazin-4-one
    Entries: EQ014201, EQ014202, EQ014203, EQ014204, EQ014205, EQ014206, EQ014207, EQ014208, EQ014209, EQ014251, EQ014252, EQ014253, EQ014254, EQ014255, EQ014256, EQ014257, EQ014258, EQ014259
    Molecule: Metoxuron 3-(3-chloro-4-methoxyphenyl)-1,1-dimethylurea
    Entries: EQ320401, EQ320402, EQ320403, EQ320404, EQ320405, EQ320406, EQ320407, EQ320408, EQ320409, EQ320451, EQ320452, EQ320453, EQ320454, EQ320455, EQ320456, EQ320457, EQ320458, EQ320459
    Molecule: Clomazone 2-(2-chlorobenzyl)-4,4-dimethyl-isoxazolidin-3-one
    Entries: EQ013901, EQ013902, EQ013903, EQ013904, EQ013905, EQ013906
    Molecule: Valganciclovir (2S)-2-amino-3-methyl-butyric acid [2-[(2-amino-6-keto-3H-purin-9-yl)methoxy]-3-hydroxy-propyl] ester
    Entries: EQ284351, EQ284352, EQ284353, EQ284354, EQ284355, EQ284356
    Molecule: Losartan [2-butyl-5-chloranyl-3-[[4-[2-(2H-1,2,3,4-tetrazol-5-yl)phenyl]phenyl]methyl]imidazol-4-yl]methanol
    Entries: EQ279401, EQ279402, EQ279403, EQ279404, EQ279405, EQ279406, EQ279407, EQ279408, EQ279409, EQ279451, EQ279452, EQ279453, EQ279454, EQ279455, EQ279456, EQ279457, EQ279458, EQ279459
    Molecule: Fexofenadine 2-[4-[1-hydroxy-4-[4-[hydroxy(diphenyl)methyl]-1-piperidinyl]butyl]phenyl]-2-methylpropanoic acid
    Entries: EQ301901, EQ301902, EQ301903, EQ301904, EQ301905, EQ301906, EQ301951, EQ301952, EQ301953, EQ301954, EQ301955, EQ301956
    Molecule: Tetracycline (4S,4aS,5aS,6S,12aS)-4-(Dimethylamino)-3,6,10,12,12a-pentahydroxy-6-methyl-1,11-dioxo-1,4,4a,5,5a,6,11,12a-octahydro-2-tetracenecarboxamide tetrex
    Entries: EQ323501, EQ323502, EQ323552, EQ323553, EQ323555
    Molecule: Simvastatin (1S,3R,7S,8S,8aR)-8-{2-[(2R,4R)-4-hydroxy-6-oxotetrahydro-2H-pyran-2-yl]ethyl}-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl 2,2-dimethylbutanoate 2,2-dimethylbutanoic acid [(1S,3R,7S,8S,8aR)-8-[2-[(2R,4R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]...
    Entries: EQ314401, EQ314402, EQ314403, EQ314404, EQ314405, EQ314406
    Molecule: 4,6-dinitro-o-cresol DNOC 2-methyl-4,6-dinitro-phenol
    Entries: EQ312051, EQ312052, EQ312053, EQ312054, EQ312055, EQ312056
    Molecule: Aliskiren (2S,4S,5S,7S)-5-amino-N-(3-amino-2,2-dimethyl-3-oxopropyl)-4-hydroxy-7-[[4-methoxy-3-(3-methoxypropoxy)phenyl]methyl]-8-methyl-2-propan-2-ylnonanamide
    Entries: EQ312301, EQ312302, EQ312303, EQ312304, EQ312305, EQ312306, EQ312351, EQ312352, EQ312353, EQ312354, EQ312355, EQ312356
    Molecule: Haloperidol 4-[4-(4-chlorophenyl)-4-hydroxypiperidin-1-yl]-1-(4-fluorophenyl)butan-1-one
    Entries: EQ356601, EQ356602, EQ356603, EQ356604, EQ356605, EQ356606, EQ356607, EQ356608, EQ356609
    Molecule: Perindopril Coversyl (2S,3aS, 7aS)-1-[(2S)-2-[[(2S)-1-ethoxy-1-oxopentan-2-yl]amino]propanoyl]-2,3,3a, 4,5,6,7,7a-octahydroindole-2-carboxylic acid
    Entries: EQ302601, EQ302602, EQ302603, EQ302604, EQ302605, EQ302606, EQ302651, EQ302652, EQ302653, EQ302654, EQ302655, EQ302656
    Molecule: Thiazopyr 2-(difluoromethyl)-4-isobutyl-5-(2-thiazolin-2-yl)-6-(trifluoromethyl)nicotinic acid methyl ester
    Entries: EQ305202, EQ305203, EQ305204, EQ305205, EQ305206
    Molecule: N-Nitrosodipropylamine (NDPA) N,N-dipropylnitrous amide
    Entries: EQ345601, EQ345602, EQ345603, EQ345604, EQ345605, EQ345606
    Molecule: Benserazide 2-amino-3-hydroxy-N`-(2,3,4-trihydroxybenzyl)propionohydrazide
    Entries: EQ285301, EQ285303, EQ285304, EQ285306, EQ285351, EQ285352, EQ285353, EQ285354, EQ285355, EQ285356
    Molecule: Pramoxine 4-[3-(4-butoxyphenoxy)propyl]morpholine
    Entries: EQ327301, EQ327302, EQ327303, EQ327304, EQ327305, EQ327306, EQ327307, EQ327308, EQ327309
    Molecule: Mepanipyrim (4-methyl-6-prop-1-ynyl-pyrimidin-2-yl)-phenyl-amine
    Entries: EQ306201, EQ306202, EQ306203, EQ306204, EQ306205, EQ306206
    Molecule: N-Nitrosofenfluramine N-Ethyl-N-nitroso-1-[3-(trifluoromethyl)phenyl]-2-propanamine
    Entries: EQ347201, EQ347202, EQ347203, EQ347204, EQ347205, EQ347206, EQ347207, EQ347208, EQ347209
    Molecule: Atenolol-desisopropyl 4-(3-Amino-2-hydroxypropoxy)phenylacetamide
    Entries: EQ267001, EQ267002, EQ267003, EQ267004, EQ267005, EQ267006, EQ267007, EQ267008, EQ267009
    Molecule: Rimantadine 1-(1-adamantyl)ethanamine
    Entries: EQ314901, EQ314902, EQ314903, EQ314904, EQ314905, EQ314906
    Molecule: 3-[(4-chlorobenzoyl)amino]propanoic acid N-(4-Chlorobenzoyl)-Beta-alanine
    Entries: EQ338401, EQ338402, EQ338403, EQ338404, EQ338405, EQ338406, EQ338407, EQ338408, EQ338409, EQ338451, EQ338452, EQ338453, EQ338454, EQ338455, EQ338456, EQ338457
    Molecule: Carbofuran (2,2-dimethyl-3H-1-benzofuran-7-yl) N-methylcarbamate
    Entries: EQ304001, EQ304002, EQ304003, EQ304004, EQ304005, EQ304006
    Molecule: Pyrilamine N`-[(4-methoxyphenyl)methyl]-N,N-dimethyl-N`-(2-pyridinyl)ethane-1,2-diamine
    Entries: EQ300601, EQ300602, EQ300603, EQ300604, EQ300605, EQ300606
    Molecule: Tramadol N-Oxide 1-[(1R,2R)-2-hydroxy-2-(3-methoxyphenyl)cyclohexyl]-N,N-dimethylmethanamine oxide
    Entries: EQ331001, EQ331002, EQ331003, EQ331004, EQ331005, EQ331006, EQ331007, EQ331008, EQ331009
    Molecule: Repaglinide 2-ethoxy-4-[2-[[3-methyl-1-(2-piperidin-1-ylphenyl)butyl]amino]-2-oxoethyl]benzoic acid
    Entries: EQ334901, EQ334902, EQ334903, EQ334904, EQ334905, EQ334906, EQ334907, EQ334908, EQ334909, EQ334951, EQ334952, EQ334953, EQ334954, EQ334955, EQ334956, EQ334957, EQ334958, EQ334959
    Molecule: Xylometazoline 2-(4-tert-butyl-2,6-dimethyl-benzyl)-2-imidazoline
    Entries: EQ303101, EQ303102, EQ303103, EQ303104, EQ303105, EQ303106
    Molecule: Lidocaine-N-Oxide Lignocaine N-oxide 2-(2,6-dimethylanilino)-N,N-diethyl-2-oxoethanamine oxide
    Entries: EQ327901, EQ327902, EQ327903, EQ327904, EQ327905, EQ327906, EQ327907, EQ327908, EQ327909
    Molecule: Buflomedil 4-(1-pyrrolidinyl)-1-(2,4,6-trimethoxyphenyl)-1-butanone
    Entries: EQ306801, EQ306802, EQ306803, EQ306804, EQ306805, EQ306806
    Molecule: Warfarin 2-hydroxy-3-(3-oxo-1-phenylbutyl)chromen-4-one Coumafen
    Entries: EQ310301, EQ310302, EQ310303, EQ310304, EQ310305, EQ310306, EQ310351, EQ310352, EQ310353, EQ310354, EQ310355, EQ310356
    Molecule: Acetazolamide AZM N-(5-sulfamoyl-1,3,4-thiadiazol-2-yl)acetamide
    Entries: EQ301101, EQ301102, EQ301103, EQ301104, EQ301105, EQ301106, EQ301151, EQ301152, EQ301153, EQ301154, EQ301155, EQ301156
    Molecule: Piracetam 1-Pyrrolidineacetamide, 2-oxo- 2-(2-ketopyrrolidino)acetamide
    Entries: EQ307601, EQ307602, EQ307603, EQ307604, EQ307605, EQ307606
    Molecule: Fenuron 1,1-dimethyl-3-phenyl-urea
    Entries: EQ311101, EQ311102, EQ311103, EQ311104, EQ311105, EQ311106
    Molecule: Uniconazole 1-(4-chlorophenyl)-4,4-dimethyl-2-(1,2,4-triazol-1-yl)pent-1-en-3-ol
    Entries: EQ358201, EQ358202, EQ358203, EQ358204, EQ358205, EQ358206, EQ358207, EQ358208, EQ358209
    Molecule: Midazolam 8-chloranyl-6-(2-fluorophenyl)-1-methyl-4H-imidazo[1,5-a][1,4]benzodiazepine
    Entries: EQ302301, EQ302302, EQ302303, EQ302304, EQ302305, EQ302306
    Molecule: Sulfentrazon N-(2,4-dichloro-5-(4-(difluoromethyl)-4,5-dihydro-3-methyl-5-oxo-1H-1,2,4-triazol-1-yl)phenyl)methanesulfonamide N-[2,4-dichloro-5-[4-(difluoromethyl)-3-methyl-5-oxo-1,2,4-triazol-1-yl]phenyl]methanesulfonamide
    Entries: EQ320901, EQ320902, EQ320903, EQ320904, EQ320905, EQ320906, EQ320907, EQ320908, EQ320909, EQ320951, EQ320952, EQ320953, EQ320954, EQ320955, EQ320956, EQ320957, EQ320958, EQ320959
    Molecule: Propranolol 1-(1-naphthalenyloxy)-3-(propan-2-ylamino)-2-propanol
    Entries: EQ017101, EQ017102, EQ017103, EQ017104, EQ017105, EQ017106
    Molecule: Fenamidone (5S)-3-anilino-5-methyl-2-(methylthio)-5-phenyl-2-imidazolin-4-one
    Entries: EQ305601, EQ305602, EQ305603, EQ305604, EQ305605, EQ305606
    Molecule: Furalaxyl 2-(N-[2-furanyl(oxo)methyl]-2,6-dimethylanilino)propanoic acid methyl ester
    Entries: EQ311301, EQ311302, EQ311303, EQ311304, EQ311305, EQ311306
    Molecule: 4-Piperidinecarboxamide 4-piperidinecarboxamide
    Entries: EQ315501, EQ315502, EQ315503, EQ315504, EQ315505, EQ315506
    Molecule: Deprenyl Selegiline (2R)-N-methyl-1-phenyl-N-prop-2-ynylpropan-2-amine
    Entries: EQ327501, EQ327502, EQ327503, EQ327504, EQ327505, EQ327506, EQ327507, EQ327508, EQ327509
    Molecule: Pipemidic acid 5,8-Dihydro-8-ethyl-5-oxo-2-(1-piperazinyl)pyrido(2,3-d)pyrimidine-6-carboxylic acid 8-ethyl-5-keto-2-piperazin-4-ium-1-yl-pyrido[2,3-d]pyrimidine-6-carboxylate
    Entries: EQ307501, EQ307502, EQ307503, EQ307504, EQ307505, EQ307506
    Molecule: Amiodarone (2-butyl-1-benzofuran-3-yl)-[4-(2-diethylaminoethyloxy)-3,5-bis(iodanyl)phenyl]methanone
    Entries: EQ306701, EQ306702, EQ306703, EQ306704, EQ306705, EQ306706
    Molecule: Oryzalin 4-(dipropylamino)-3,5-dinitro-benzenesulfonamide
    Entries: EQ309901, EQ309902, EQ309903, EQ309904, EQ309905, EQ309906, EQ309951, EQ309952, EQ309953, EQ309954, EQ309955, EQ309956
    Molecule: Imatinib 4-[(4-methylpiperazin-1-yl)methyl]-N-[4-methyl-3-[(4-pyridin-3-ylpyrimidin-2-yl)amino]phenyl]benzamide
    Entries: EQ363301, EQ363302, EQ363303, EQ363304, EQ363305, EQ363306, EQ363307, EQ363308, EQ363309, EQ363351, EQ363352, EQ363353, EQ363354, EQ363355, EQ363356, EQ363357, EQ363358, EQ363359
    Molecule: N-Nitrosomethylaniline N-methyl-N-phenylnitrous amide
    Entries: EQ335701, EQ335702, EQ335703, EQ335704, EQ335705, EQ335706
    Molecule: N-Nitrosodiphenylamine N,N-diphenylnitrous amide
    Entries: EQ335501, EQ335502, EQ335503, EQ335504, EQ335505, EQ335506
    Molecule: 4`-Hydroxy Diclofenac 4`-hydroxydiclofenac 2-[2-(2,6-dichloro-4-hydroxyanilino)phenyl]acetic acid
    Entries: EQ339701, EQ339702, EQ339703, EQ339704, EQ339705, EQ339706, EQ339707, EQ339708, EQ339709, EQ339751, EQ339752, EQ339753, EQ339754, EQ339755, EQ339756, EQ339757, EQ339758, EQ339759
    Molecule: 2,6-Di-tert-butylpyridine 2,6-ditert-butylpyridine 2,6-Bis(2-methyl-2-propanyl)pyridine
    Entries: EQ359001, EQ359002, EQ359003, EQ359004, EQ359005, EQ359006, EQ359007, EQ359008, EQ359009
    Molecule: Gabapentin Related Compound E 1-(carboxymethyl)cyclohexanecarboxylic acid 1-(carboxymethyl)cyclohexane-1-carboxylic acid
    Entries: EQ331201, EQ331202, EQ331203, EQ331204, EQ331205, EQ331206, EQ331207, EQ331208, EQ331209, EQ331251, EQ331252, EQ331253, EQ331254, EQ331255, EQ331256, EQ331257
    Molecule: Triphenylphosphine oxide Triphenylphosphineoxide (diphenylphosphoroso)benzene
    Entries: EQ358701, EQ358702, EQ358703, EQ358704, EQ358705, EQ358706, EQ358707, EQ358708, EQ358709
    Molecule: Doxylamine N,N-dimethyl-2-(1-phenyl-1-pyridin-2-yl-ethoxy)ethanamine
    Entries: EQ301801, EQ301802, EQ301803, EQ301804, EQ301805, EQ301806

Authors: Vogler B, Stravs M, Schymanski E, Singer H, Department of Environmental Chemistry, Eawag   
CopyRight: Copyright (C) 2015 Eawag, Duebendorf, Switzerland   

    Molecule: 2-Chlorobenzenesulfonamide 2-chlorobenzenesulfonamide
    Entries: EQ326001, EQ326002, EQ326003, EQ326004, EQ326005, EQ326006
    Molecule: 4-Amino-6-chloro-1,3-benzenedisulfonamide 4-amino-6-chlorobenzene-1,3-disulfonamide
    Entries: EQ325701, EQ325704, EQ325705, EQ325706, EQ325751, EQ325752
    Molecule: Progesterone (8S,9S,10R,13S,14S,17S)-17-acetyl-10,13-dimethyl-1,2,6,7,8,9,11,12,14,15,16,17-dodecahydrocyclopenta[a]phenanthren-3-one
    Entries: EQ325502, EQ325503, EQ325504, EQ325505, EQ325506
    Molecule: Guanylurea diaminomethylideneurea
    Entries: EQ325301, EQ325302, EQ325303, EQ325304
    Molecule: 4-Chlorobenzenesulfonamide 4-chlorobenzenesulfonamide
    Entries: EQ325951, EQ325954, EQ325955, EQ325956

Authors: Huntscha S, Stravs M, Schymanski E, Singer H, Department of Environmental Chemistry, Eawag   
CopyRight: Copyright (C) 2015 Eawag, Duebendorf, Switzerland   

    Molecule: Dimefuron 3-[4-(5-tert-butyl-2-keto-1,3,4-oxadiazol-3-yl)-3-chloro-phenyl]-1,1-dimethyl-urea
    Entries: EQ316601, EQ316602, EQ316603, EQ316604, EQ316605, EQ316606, EQ316651, EQ316652, EQ316653, EQ316654, EQ316655, EQ316656
    Molecule: Oxyfluorfen 2-chloranyl-1-(3-ethoxy-4-nitro-phenoxy)-4-(trifluoromethyl)benzene
    Entries: EQ317401, EQ317402, EQ317403, EQ317404, EQ317405, EQ317406
    Molecule: Monolinuron 3-(4-chlorophenyl)-1-methoxy-1-methyl-urea
    Entries: EQ317301, EQ317302, EQ317303, EQ317304, EQ317305, EQ317306
    Molecule: Tepraloxydim E-Tepraloxydim
    Entries: EQ317801, EQ317802, EQ317803, EQ317804, EQ317805, EQ317806, EQ317851, EQ317852, EQ317853, EQ317854, EQ317855, EQ317856
    Molecule: Prothioconazole-desethio 1H-1,2,4-Triazole-1-ethanol, alpha-(1-chlorocyclopropyl)-alpha-((2-chlorophenyl)methyl)- 2-(1-chloranylcyclopropyl)-1-(2-chlorophenyl)-3-(1,2,4-triazol-1-yl)propan-2-ol
    Entries: EQ317601, EQ317602, EQ317603, EQ317604, EQ317605, EQ317606, EQ317651, EQ317652, EQ317653, EQ317654, EQ317655, EQ317656
    Molecule: Dichlorvos DDVP 2,2-bis(chloranyl)ethenyl dimethyl phosphate
    Entries: EQ304701, EQ304702, EQ304703, EQ304704, EQ304705, EQ304706, EQ304751, EQ304752, EQ304753, EQ304754, EQ304755, EQ304756
    Molecule: Terbacil 3-tert-butyl-5-chloranyl-6-methyl-1H-pyrimidine-2,4-dione
    Entries: EQ317951, EQ317952, EQ317953, EQ317954, EQ317955, EQ317956
    Molecule: Flumioxazin 2-(7-fluoranyl-3-oxidanylidene-4-prop-2-ynyl-1,4-benzoxazin-6-yl)-4,5,6,7-tetrahydroisoindole-1,3-dione
    Entries: EQ316701, EQ316702, EQ316703, EQ316704, EQ316705, EQ316706
    Molecule: Metosulam N-(2,6-dichloro-3-methyl-phenyl)-5,7-dimethoxy-[1,2,4]triazolo[1,5-a]pyrimidine-2-sulfonamide
    Entries: EQ317201, EQ317202, EQ317203, EQ317204, EQ317205, EQ317206, EQ317251, EQ317252, EQ317253, EQ317254, EQ317255, EQ317256
    Molecule: 1H-1,2,3-triazole-5-OH 1H-1,2,3-Triazol-4-ol
    Entries: EQ318601, EQ318602, EQ318603, EQ318604, EQ318605, EQ318606, EQ318651, EQ318652, EQ318653, EQ318654, EQ318655, EQ318656
    Molecule: Metconazole 5-(4-chlorobenzyl)-2,2-dimethyl-1-(1,2,4-triazol-1-ylmethyl)cyclopentanol
    Entries: EQ317101, EQ317102, EQ317103, EQ317104, EQ317105, EQ317106
    Molecule: Mesosulfuron-methyl 2-[(4,6-dimethoxypyrimidin-2-yl)carbamoylsulfamoyl]-4-(methanesulfonamidomethyl)benzoic acid methyl ester
    Entries: EQ317001, EQ317002, EQ317003, EQ317004, EQ317005, EQ317006, EQ317051, EQ317052, EQ317053, EQ317054, EQ317055, EQ317056
    Molecule: 1H-1,2,3-Triazole 2H-1,2,3-triazole
    Entries: EQ318501, EQ318502, EQ318503, EQ318504, EQ318505, EQ318506
    Molecule: 4,5-dihydro-5,5-diphenyl-1,2-oxazole-3-carboxylic acid, ethyl ester 5,5-diphenyl-2-isoxazoline-3-carboxylic acid ethyl ester
    Entries: EQ316801, EQ316802, EQ316803, EQ316804, EQ316805, EQ316806
    Molecule: Thiabendazole 2-(4-Thiazolyl)benzimidazole 4-(1H-benzimidazol-2-yl)-1,3-thiazole
    Entries: EQ318001, EQ318002, EQ318003, EQ318004, EQ318005, EQ318006

Authors: Lege S, Zwiener C, Environmental Analytical Chemistry (EAC), University of Tuebingen   
CopyRight: Copyright (C) 2015, Environmental Analytical Chemistry (EAC), University of Tuebingen, Germany   

    Molecule: Benzophenone-4 Sulisobenzone
    Entries: TUE00145, TUE00146, TUE00147
    Molecule: Atenolol
    Entries: TUE00081, TUE00082, TUE00083
    Molecule: 2-Aminobenzimidazol
    Entries: TUE00011, TUE00012, TUE00013, TUE00015, TUE00016, TUE00017
    Molecule: Diethyltoluamide DEET
    Entries: TUE00311, TUE00312, TUE00313
    Molecule: Sulfamethoxazole
    Entries: TUE00671, TUE00672, TUE00673, TUE00675, TUE00676, TUE00677
    Molecule: Melamine
    Entries: TUE00491, TUE00492, TUE00493
    Molecule: Bisoprolol
    Entries: TUE00191, TUE00192, TUE00193
    Molecule: Phenylbenzimidazole sulfonic acid
    Entries: TUE00581, TUE00582, TUE00583, TUE00585, TUE00586, TUE00587
    Molecule: Sotalol
    Entries: TUE00651, TUE00652, TUE00653
    Molecule: Benzophenone-2 2,2',4,4'-tetrahydroxybenzophenone
    Entries: TUE00121, TUE00122, TUE00123, TUE00125, TUE00126, TUE00127
    Molecule: 5-Methylbenzotriazole
    Entries: TUE00031, TUE00032, TUE00033, TUE00035, TUE00036, TUE00037
    Molecule: O-desmethylvenlafaxine
    Entries: TUE00821, TUE00822, TUE00823, TUE00825, TUE00826, TUE00827
    Molecule: Venlafaxine
    Entries: TUE00811, TUE00812, TUE00813
    Molecule: 4-Methylbenzotriazole
    Entries: TUE00041, TUE00042, TUE00043, TUE00045, TUE00046, TUE00047
    Molecule: Guanylurea
    Entries: TUE00361, TUE00362, TUE00363
    Molecule: Hydrochlorothiazide
    Entries: TUE00375, TUE00376, TUE00377
    Molecule: Metoprolol
    Entries: TUE00531, TUE00532, TUE00533
    Molecule: Gabapentin
    Entries: TUE00341, TUE00342, TUE00343, TUE00345, TUE00346, TUE00347
    Molecule: Sucralose
    Entries: TUE00661, TUE00662, TUE00665, TUE00666, TUE00667
    Molecule: Tramadol
    Entries: TUE00731, TUE00732, TUE00733
    Molecule: Lidocaine Lignocaine
    Entries: TUE00461, TUE00462, TUE00463
    Molecule: Benzotriazole
    Entries: TUE00161, TUE00162, TUE00163, TUE00165, TUE00166, TUE00167

Authors: Stravs M, Schymanski E, Singer H, Department of Environmental Chemistry, Eawag   
CopyRight: Copyright (C) 2012 Eawag, Duebendorf, Switzerland   

    Molecule: Acetochlor ESA 2-[(Ethoxymethyl)(2-ethyl-6-methylphenyl)amino]-2-oxoethanesulfonic acid
    Entries: EA010351, EA010352, EA010353, EA010354, EA010355, EA010356, EA010357, EA010358, EA010359, EA010360, EA010361, EA010362, EA010363, EA010364
    Molecule: Chlortoluron Chlorotoluron 3-(3-chloranyl-4-methyl-phenyl)-1,1-dimethyl-urea
    Entries: EA027301, EA027302, EA027303, EA027304, EA027305, EA027306, EA027307, EA027308, EA027309, EA027310, EA027311, EA027312, EA027313, EA027314
    Molecule: Terbutylazine-desethyl-2-hydroxy 1,3,5-Triazin-2(1H)-one, 4-amino-6-((1,1-dimethylethyl)amino)- 2-amino-6-(tert-butylamino)-1H-1,3,5-triazin-4-one
    Entries: EA067201, EA067202, EA067203, EA067204, EA067205, EA067206, EA067207, EA067208, EA067209, EA067210, EA067211, EA067212, EA067213, EA067214, EA067251, EA067252, EA067253, EA067254, EA067255, EA067256, EA067257, EA067258, EA067259, EA067260, EA067261, EA067262
    Molecule: Dimethenamid OXA [(2,4-dimethylthiophen-3-yl)(1-methoxypropan-2-yl)amino](oxo)acetic acid
    Entries: EA025901, EA025902, EA025903, EA025904, EA025905, EA025906, EA025907, EA025908, EA025909, EA025910, EA025911, EA025912, EA025913, EA025914, EA025951, EA025952, EA025953, EA025954, EA025955, EA025956, EA025957, EA025958, EA025959, EA025960, EA025961, EA025962, EA025963, EA025964
    Molecule: Fipronil-sulfone 5-amino-1-[2,6-dichloro-4-(trifluoromethyl)phenyl]-4-(trifluoromethylsulfonyl)-3-pyrazolecarbonitrile
    Entries: EA266751, EA266752, EA266753, EA266754, EA266755, EA266756, EA266757, EA266758, EA266759, EA266760, EA266761, EA266762, EA266763, EA266764
    Molecule: Perfluorohexanoic acid 2,2,3,3,4,4,5,5,6,6,6-undecafluorohexanoic acid
    Entries: EA271551, EA271552, EA271553, EA271554, EA271555, EA271557, EA271558, EA271559
    Molecule: Ioxynil 4-Hydroxy-3,5-diiodobenzonitrile 3,5-bis(iodanyl)-4-oxidanyl-benzenecarbonitrile
    Entries: EA008751, EA008752, EA008753, EA008754, EA008755, EA008756, EA008757, EA008758, EA008759, EA008760, EA008761, EA008762, EA008763, EA008764
    Molecule: Dextromethorphan DXM (14alpha)-3-Methoxy-17-methylmorphinan
    Entries: EA282401, EA282402, EA282403, EA282404, EA282405, EA282406, EA282407, EA282408, EA282409, EA282410, EA282411, EA282412, EA282413, EA282414
    Molecule: Pyrimidinol 2-isopropyl-6-methyl-1H-pyrimidin-4-one
    Entries: EA270601, EA270602, EA270603, EA270604, EA270605, EA270606, EA270607, EA270608, EA270609, EA270610, EA270611, EA270612, EA270613, EA270614
    Molecule: Acetochlor OXA [(Ethoxymethyl)(2-ethyl-6-methylphenyl)amino](oxo)acetic acid
    Entries: EA010251, EA010252, EA010253, EA010254, EA010255, EA010256, EA010257, EA010258, EA010259, EA010260, EA010261, EA010262, EA010263, EA010264
    Molecule: Fenpropidin 1-[3-(4-tert-butylphenyl)-2-methyl-propyl]piperidine
    Entries: EA295801, EA295802, EA295803, EA295804, EA295805, EA295806, EA295807, EA295808, EA295809, EA295810, EA295811, EA295812, EA295813, EA295814
    Molecule: Cyanazine 2-[[4-chloranyl-6-(ethylamino)-1,3,5-triazin-2-yl]amino]-2-methyl-propanenitrile
    Entries: EA275901, EA275902, EA275903, EA275904, EA275905, EA275906, EA275907, EA275908, EA275909, EA275910, EA275911, EA275912, EA275913, EA275914
    Molecule: Gemcitabine 2'-deoxy-2',2'-difluorocytidine 4-amino-1-[(2R,4R,5R)-3,3-difluoro-4-hydroxy-5-(hydroxymethyl)-2-oxolanyl]-2-pyrimidinone
    Entries: EA260301, EA260302, EA260303, EA260304, EA260305, EA260306, EA260307, EA260308, EA260309, EA260310, EA260311, EA260312, EA260313, EA260314, EA260351, EA260352, EA260353, EA260354, EA260355, EA260356, EA260357, EA260358, EA260359, EA260360, EA260361, EA260362, EA260363, EA260364
    Molecule: Thiacloprid-amide [3-(6-Chloro-3-pyridylmethyl)thiazolidin-2-ylidene]urea
    Entries: EA295501, EA295502, EA295503, EA295504, EA295505, EA295506, EA295507, EA295508, EA295509, EA295510, EA295511, EA295512, EA295513, EA295514
    Molecule: Carbamazepine 5H-Dibenz[b,f]azepine-5-carboxamide 11-benzo[b][1]benzazepinecarboxamide
    Entries: EA019401, EA019402, EA019403, EA019404, EA019405, EA019406, EA019407, EA019408, EA019409, EA019410, EA019411, EA019412, EA019413, EA019414
    Molecule: Ciprofloxacin 1-cyclopropyl-6-fluoro-4-keto-7-piperazin-4-ium-1-yl-quinoline-3-carboxylate
    Entries: EA027601, EA027602, EA027603, EA027604, EA027605, EA027606, EA027607, EA027608, EA027609, EA027610, EA027611, EA027612, EA027613, EA027614
    Molecule: Teflubenzuron N-[(3,5-dichloro-2,4-difluoro-phenyl)carbamoyl]-2,6-difluoro-benzamide
    Entries: EA293101, EA293102, EA293103, EA293104, EA293105, EA293106, EA293107, EA293108, EA293109, EA293110, EA293111, EA293112, EA293113, EA293114, EA293153, EA293155, EA293156
    Molecule: 4-Isopropylaniline 4-propan-2-ylaniline
    Entries: EA005201, EA005202, EA005203, EA005204, EA005205, EA005206, EA005207, EA005208, EA005209, EA005210, EA005211, EA005212, EA005213, EA005214
    Molecule: Bicalutamide (2R)-N-[4-cyano-3-(trifluoromethyl)phenyl]-3-[(4-fluorophenyl)sulfonyl]-2-hydroxy-2-methylpropanamide
    Entries: EA280901, EA280902, EA280903, EA280904, EA280905, EA280906, EA280907, EA280908, EA280909, EA280910, EA280911, EA280912, EA280913, EA280914, EA280951, EA280952, EA280953, EA280954, EA280955, EA280956, EA280957, EA280958, EA280959, EA280960, EA280961, EA280962, EA280963, EA280964
    Molecule: Mefenamic acid 2-(2,3-dimethylanilino)benzoic acid
    Entries: EA020801, EA020802, EA020803, EA020804, EA020805, EA020806, EA020807, EA020808, EA020809, EA020810, EA020811, EA020812, EA020813, EA020814, EA020851, EA020852, EA020853, EA020854, EA020855, EA020856, EA020857, EA020858, EA020859, EA020860, EA020861, EA020862, EA020863, EA020864
    Molecule: Albuterol Salbutamol 4-[2-(tert-butylamino)-1-hydroxy-ethyl]-2-methylol-phenol
    Entries: EA285101, EA285102, EA285103, EA285104, EA285105, EA285106, EA285107, EA285108, EA285109, EA285110, EA285111, EA285112, EA285113, EA285114, EA285151, EA285152, EA285153, EA285154, EA285155, EA285156, EA285158, EA285159, EA285160, EA285161, EA285162, EA285164
    Molecule: Cetirizine 2-[2-[4-[(4-chlorophenyl)-phenyl-methyl]piperazin-1-yl]ethoxy]ethanoic acid
    Entries: EA277201, EA277202, EA277203, EA277204, EA277205, EA277206, EA277207, EA277208, EA277209, EA277210, EA277211, EA277212, EA277213, EA277214, EA277251, EA277252, EA277253, EA277254, EA277258, EA277259, EA277260, EA277264
    Molecule: Fipronil-sulfide 5-Amino-1-[2,6-dichloro-4-(trifluoromethyl)phenyl]-4-[(trifluoromethyl)sulfanyl]-1H-pyrazole-3-carbonitrile
    Entries: EA266801, EA266802, EA266803, EA266804, EA266805, EA266806, EA266807, EA266808, EA266809, EA266810, EA266811, EA266812, EA266813, EA266814, EA266853, EA266854, EA266859, EA266860
    Molecule: N,N-Didesmethylvenlafaxine 1-[2-Amino-1-(4-methoxyphenyl)ethyl]cyclohexanol
    Entries: EA265801, EA265802, EA265803, EA265804, EA265805, EA265806, EA265807, EA265808, EA265809, EA265810, EA265811, EA265812, EA265813, EA265814
    Molecule: Benzisothiazolone (BIT) 1,2-Benzisothiazol-3(2H)-one 1,2-benzothiazol-3-one
    Entries: EA030101, EA030102, EA030103, EA030104, EA030105, EA030106, EA030107, EA030108, EA030109, EA030110, EA030111, EA030112, EA030113, EA030114, EA030152, EA030153, EA030154, EA030158, EA030159
    Molecule: Metsulfuron-methyl 2-[(4-methoxy-6-methyl-s-triazin-2-yl)carbamoylsulfamoyl]benzoic acid methyl ester
    Entries: EA012801, EA012802, EA012803, EA012804, EA012805, EA012806, EA012807, EA012808, EA012809, EA012810, EA012811, EA012812, EA012813, EA012814, EA012851, EA012852, EA012853, EA012854, EA012855, EA012856, EA012857, EA012858, EA012859, EA012860, EA012861, EA012862, EA012863, EA012864
    Molecule: Isoproturon 1,1-dimethyl-3-(4-propan-2-ylphenyl)urea
    Entries: EA028601, EA028602, EA028603, EA028604, EA028605, EA028606, EA028607, EA028608, EA028609, EA028610, EA028611, EA028612, EA028613, EA028614
    Molecule: Triclopyr 2-(3,5,6-trichloropyridin-2-yl)oxyacetic acid
    Entries: EA294601, EA294602, EA294603, EA294604, EA294605, EA294606, EA294607, EA294608, EA294609, EA294610, EA294611, EA294612, EA294613, EA294614
    Molecule: Perfluorodecyl phosphate 8:2PAP 1H,1H,2H,2H-perfluorodecylphosphate
    Entries: EA292502, EA292503, EA292504, EA292505, EA292506, EA292507, EA292508, EA292509, EA292510, EA292511, EA292512, EA292513
    Molecule: Oxazepam 7-chloranyl-3-oxidanyl-5-phenyl-1,3-dihydro-1,4-benzodiazepin-2-one
    Entries: EA274301, EA274302, EA274303, EA274304, EA274305, EA274306, EA274307, EA274308, EA274309, EA274310, EA274311, EA274312, EA274313, EA274314, EA274351, EA274352, EA274353, EA274354, EA274356, EA274358, EA274359, EA274360, EA274361, EA274362, EA274364
    Molecule: Acetamiprid (1E)-N-[(6-Chloro-3-pyridinyl)methyl]-N'-cyano-N-methylethanimidamide
    Entries: EA298601, EA298602, EA298603, EA298604, EA298605, EA298606, EA298607, EA298608, EA298609, EA298610, EA298611, EA298612, EA298613, EA298614, EA298651, EA298652, EA298653, EA298658
    Molecule: 2-Naphthalenesulfonic acid 2-naphthalenesulfonic acid
    Entries: EA065301, EA065302, EA065303, EA065304, EA065305, EA065306, EA065307, EA065308, EA065309, EA065310, EA065311, EA065312, EA065313, EA065314, EA065351, EA065352, EA065353, EA065354, EA065355, EA065356, EA065357, EA065358, EA065359, EA065360, EA065361, EA065362, EA065363
    Molecule: Dimethachlor ESA 2-[(2,6-dimethylphenyl)(2-methoxyethyl)amino]-2-oxoethanesulfonic acid
    Entries: EA253501, EA253502, EA253503, EA253504, EA253505, EA253506, EA253507, EA253508, EA253509, EA253510, EA253511, EA253512, EA253513, EA253514, EA253551, EA253552, EA253553, EA253554, EA253555, EA253556, EA253557, EA253558, EA253559, EA253560, EA253561, EA253562, EA253563, EA253564
    Molecule: N'-(2,4-Dimethylphenyl)-N-methylformamidine Methanimidamide, N-(2,4-dimethylphenyl)-N'-methyl- N-(2,4-dimethylphenyl)-N'-methyl-formamidine
    Entries: EA033101, EA033102, EA033103, EA033104, EA033105, EA033106, EA033107, EA033108, EA033109, EA033110, EA033111, EA033112, EA033113, EA033114
    Molecule: Methamphetamine (2S)-N-methyl-1-phenylpropan-2-amine
    Entries: EA282901, EA282902, EA282903, EA282904, EA282905, EA282906, EA282907, EA282908, EA282909, EA282910, EA282911, EA282912, EA282913, EA282914
    Molecule: Difenoconazole 1-[[2-[2-chloranyl-4-(4-chloranylphenoxy)phenyl]-4-methyl-1,3-dioxolan-2-yl]methyl]-1,2,4-triazole
    Entries: EA293401, EA293402, EA293403, EA293404, EA293405, EA293406, EA293407, EA293408, EA293409, EA293410, EA293411, EA293412, EA293413, EA293414
    Molecule: N4-Acetylsulfadimethoxine N-[4-[(4,6-dimethoxy-2-pyrimidinyl)sulfamoyl]phenyl]acetamide
    Entries: EA024501, EA024502, EA024503, EA024504, EA024505, EA024506, EA024507, EA024508, EA024509, EA024510, EA024511, EA024512, EA024513, EA024514, EA024551, EA024553, EA024554, EA024555, EA024556, EA024557, EA024558, EA024559, EA024560, EA024561, EA024562, EA024563
    Molecule: Fluoxetine N-methyl-3-phenyl-3-[4-(trifluoromethyl)phenoxy]-1-propanamine
    Entries: EA033401, EA033402, EA033407, EA033408, EA033412, EA033413, EA033414
    Molecule: Atenolol 2-[4-[2-hydroxy-3-(isopropylamino)propoxy]phenyl]acetamide
    Entries: EA016901, EA016902, EA016903, EA016904, EA016905, EA016906, EA016907, EA016908, EA016909, EA016910, EA016911, EA016912, EA016913, EA016914
    Molecule: Diuron-desdimethyl Urea, (3,4-dichlorophenyl)- (3,4-dichlorophenyl)urea
    Entries: EA028901, EA028902, EA028903, EA028904, EA028905, EA028906, EA028907, EA028908, EA028909, EA028910, EA028911, EA028912, EA028913, EA028914, EA028951, EA028952, EA028953, EA028954, EA028955, EA028956, EA028958, EA028959, EA028960, EA028961, EA028962
    Molecule: Trinexapac-ethyl 4-[cyclopropyl(hydroxy)methylene]-3,5-diketo-cyclohexanecarboxylic acid ethyl ester
    Entries: EA015101, EA015102, EA015103, EA015104, EA015105, EA015106, EA015107, EA015108, EA015109, EA015110, EA015111, EA015112, EA015113, EA015114
    Molecule: Ioxitalamic acid 3-acetamido-5-(2-hydroxyethylcarbamoyl)-2,4,6-triiodo-benzoic acid
    Entries: EA034201, EA034202, EA034203, EA034204, EA034205, EA034206, EA034207, EA034208, EA034209, EA034210, EA034211, EA034212, EA034213, EA034214
    Molecule: Prosulfocarb N,N-dipropylcarbamothioic acid S-(phenylmethyl) ester
    Entries: EA007101, EA007102, EA007103, EA007104, EA007105, EA007106, EA007107, EA007108, EA007109, EA007110, EA007111, EA007112, EA007113, EA007114
    Molecule: Flufenacet N-(4-fluorophenyl)-N-isopropyl-2-[[5-(trifluoromethyl)-1,3,4-thiadiazol-2-yl]oxy]acetamide
    Entries: EA070901, EA070902, EA070903, EA070904, EA070905, EA070906, EA070907, EA070908, EA070909, EA070910, EA070911, EA070912, EA070913, EA070914
    Molecule: Acetochlor 2-chloranyl-N-(ethoxymethyl)-N-(2-ethyl-6-methyl-phenyl)ethanamide
    Entries: EA010401, EA010402, EA010403, EA010404, EA010405, EA010406, EA010407, EA010408, EA010409, EA010410, EA010411, EA010412, EA010413, EA010414
    Molecule: Sulcotrione 2-(2-chloranyl-4-methylsulfonyl-phenyl)carbonylcyclohexane-1,3-dione
    Entries: EA029001, EA029002, EA029003, EA029004, EA029005, EA029006, EA029007, EA029008, EA029009, EA029010, EA029011, EA029012, EA029013, EA029014, EA029051, EA029052, EA029053, EA029054, EA029055, EA029056, EA029058, EA029059, EA029060, EA029061, EA029062, EA029063, EA029064
    Molecule: Carbamazepine-10,11-epoxide 1a,10b-Dihydro-6H-dibenz(b,f)oxiren(d)azepine-6-carboxamide
    Entries: EA091601, EA091602, EA091603, EA091604, EA091605, EA091606, EA091607, EA091608, EA091609, EA091610, EA091611, EA091612, EA091613, EA091614
    Molecule: Carbetamide N-phenylcarbamic acid [1-(ethylamino)-1-oxopropan-2-yl] ester
    Entries: EA013601, EA013602, EA013603, EA013604, EA013605, EA013606, EA013607, EA013608, EA013609, EA013610, EA013611, EA013612, EA013613, EA013614, EA013651, EA013652, EA013653, EA013654, EA013655, EA013658, EA013659, EA013660, EA013661, EA013664
    Molecule: Moclobemide 4-chloranyl-N-(2-morpholin-4-ylethyl)benzamide
    Entries: EA267501, EA267502, EA267503, EA267504, EA267505, EA267506, EA267507, EA267508, EA267509, EA267510, EA267511, EA267512, EA267513, EA267514
    Molecule: Metolachlor OXA 2-(2-ethyl-N-(2-methoxy-1-methyl-ethyl)-6-methyl-anilino)-2-keto-acetic acid
    Entries: EA026501, EA026502, EA026503, EA026504, EA026505, EA026506, EA026507, EA026508, EA026509, EA026510, EA026511, EA026512, EA026513, EA026514, EA026551, EA026552, EA026553, EA026554, EA026555, EA026556, EA026557, EA026558, EA026559, EA026560, EA026561, EA026562, EA026563, EA026564
    Molecule: Mycophenolic acid (E)-6-(4-hydroxy-3-keto-6-methoxy-7-methyl-phthalan-5-yl)-4-methyl-hex-4-enoic acid
    Entries: EA280801, EA280802, EA280803, EA280804, EA280805, EA280806, EA280807, EA280808, EA280809, EA280810, EA280811, EA280812, EA280813, EA280814, EA280851, EA280852, EA280853, EA280854, EA280855, EA280856, EA280857, EA280858, EA280859, EA280860, EA280861, EA280862, EA280863, EA280864
    Molecule: Iobitridol 1,3-Benzenedicarboxamide, N,N'-bis(2,3-dihydroxypropyl)-5-((5-hydroxy-2-(hydroxymethyl)-1-oxopropyl)amino)-2,4,6-triiodo-N,N'-dimethyl- 1-N,3-N-bis(2,3-dihydroxypropyl)-5-[[3-hydroxy-2-(hydroxymethyl)propanoyl]amino]-2,4,6-triiodo-1-N,3-N-dime...
    Entries: EA289801, EA289802, EA289803, EA289804, EA289805, EA289806, EA289807, EA289808, EA289809, EA289810, EA289811, EA289812, EA289813, EA289814
    Molecule: Tritosulfuron N-{[4-methoxy-6-(trifluoromethyl)-1,3,5-triazin-2-yl]carbamoyl}-2-(trifluoromethyl)benzenesulfonamide
    Entries: EA274001, EA274002, EA274003, EA274004, EA274005, EA274006, EA274007, EA274008, EA274009, EA274010, EA274011, EA274012, EA274013, EA274014, EA274051, EA274052, EA274053, EA274054, EA274055, EA274056, EA274057, EA274058, EA274059, EA274060, EA274061, EA274062, EA274063, EA274064
    Molecule: 1-(3-Chlorophenyl)piperazine
    Entries: EA281801, EA281802, EA281803, EA281804, EA281805, EA281806, EA281807, EA281808, EA281809, EA281810, EA281811, EA281812, EA281813, EA281814
    Molecule: Pinoxaden 8-(2,6-diethyl-4-methylphenyl)-1,2,4,5-tetrahydro-7-oxo-7H-pyrazolo(1,2-d)(1,4,5)oxadiazepin-9-yl 2,2-dimethylpropanoate
    Entries: EA070401, EA070402, EA070403, EA070404, EA070405, EA070406, EA070407, EA070408, EA070409, EA070410, EA070411, EA070412, EA070413, EA070414
    Molecule: N,O-Didesmethylvenlafaxine 4-[1-(1-Hydroxycyclohexyl)-2-(methylamino)ethyl]phenol
    Entries: EA265701, EA265702, EA265703, EA265704, EA265705, EA265706, EA265707, EA265708, EA265709, EA265710, EA265711, EA265712, EA265713, EA265714
    Molecule: Ranitidine 1,1-Ethenediamine, N-(2-(((5-((dimethylamino)methyl)-2-furanyl)methyl)thio)ethyl)-N'-methyl-2-nitro- 1-N'-[2-[[5-(dimethylaminomethyl)furan-2-yl]methylsulfanyl]ethyl]-1-N-methyl-2-nitroethene-1,1-diamine
    Entries: EA019601, EA019602, EA019603, EA019604, EA019605, EA019606, EA019607, EA019608, EA019609, EA019610, EA019611, EA019612, EA019613, EA019614, EA019651, EA019652, EA019653, EA019654, EA019655, EA019656, EA019657, EA019658, EA019659, EA019660, EA019661, EA019662, EA019664
    Molecule: MCPA 2-(4-chloranyl-2-methyl-phenoxy)ethanoic acid
    Entries: EA026151, EA026152, EA026153, EA026154, EA026155, EA026156, EA026157, EA026158, EA026159, EA026160, EA026161, EA026162, EA026163, EA026164
    Molecule: Diazoxon diethyl (6-methyl-2-propan-2-yl-pyrimidin-4-yl) phosphate
    Entries: EA270501, EA270502, EA270503, EA270504, EA270505, EA270506, EA270507, EA270508, EA270509, EA270510, EA270511, EA270512, EA270513, EA270514
    Molecule: Thiamethoxam 3-((2-Chloro-5-thiazolyl)methyl)tetrahydro-5-methyl-N-nitro-4H-1,3,5-oxadiazin-4-imine
    Entries: EA294101, EA294102, EA294103, EA294104, EA294105, EA294106, EA294107, EA294108, EA294109, EA294110, EA294111, EA294112, EA294113, EA294114
    Molecule: Metribuzin 1,2,4-Triazin-5(4H)-one, 4-amino-6-(1,1-dimethylethyl)-3-(methylthio)- 4-amino-6-tert-butyl-3-(methylthio)-1,2,4-triazin-5-one
    Entries: EA009001, EA009002, EA009003, EA009004, EA009005, EA009006, EA009007, EA009008, EA009009, EA009010, EA009011, EA009012, EA009013, EA009014
    Molecule: Fluazinam 3-chloranyl-N-[3-chloranyl-2,6-dinitro-4-(trifluoromethyl)phenyl]-5-(trifluoromethyl)pyridin-2-amine
    Entries: EA011951, EA011952, EA011953, EA011954, EA011955, EA011956, EA011958, EA011959, EA011960, EA011961, EA011962, EA011963, EA011964
    Molecule: Clofibric acid 2-(p-Chlorophenoxy)isobutyric acid 2-(4-chloranylphenoxy)-2-methyl-propanoic acid
    Entries: EA020401, EA020403, EA020404, EA020409, EA020411, EA020413, EA020451, EA020452, EA020453, EA020454, EA020455, EA020456, EA020457, EA020458, EA020459, EA020460, EA020461, EA020462, EA020463
    Molecule: Benzoylecgonine 3-benzoyloxy-8-methyl-8-azabicyclo[3.2.1]octane-4-carboxylic acid
    Entries: EA282301, EA282302, EA282303, EA282304, EA282305, EA282306, EA282307, EA282308, EA282309, EA282310, EA282311, EA282312, EA282313, EA282314
    Molecule: Fluazifop 2-[4-[5-(trifluoromethyl)pyridin-2-yl]oxyphenoxy]propanoic acid
    Entries: EA014701, EA014702, EA014703, EA014704, EA014705, EA014706, EA014707, EA014708, EA014709, EA014710, EA014711, EA014712, EA014713, EA014714, EA014751, EA014752, EA014753, EA014754, EA014755, EA014756, EA014757, EA014758, EA014759, EA014760, EA014761, EA014762, EA014763, EA014764
    Molecule: 5-Chloro-2-methyl-4-isothiazolin-3-one (CMI) 3(2H)-Isothiazolone, 5-chloro-2-methyl- 5-chloranyl-2-methyl-1,2-thiazol-3-one
    Entries: EA028001, EA028002, EA028003, EA028004, EA028005, EA028006, EA028007, EA028008, EA028009, EA028010, EA028011, EA028012, EA028013, EA028014
    Molecule: Trimipramine 3-(5,6-dihydrobenzo[b][1]benzazepin-11-yl)-N,N,2-trimethyl-1-propanamine
    Entries: EA286101, EA286102, EA286103, EA286104, EA286105, EA286106, EA286107, EA286108, EA286109, EA286110, EA286111, EA286112, EA286113, EA286114
    Molecule: Ranitidine-S-oxide 1,1-Ethenediamine, N-(2-(((5-((dimethylamino)methyl)-2-furanyl)methyl)sulfinyl)ethyl)-N'-methyl-2-nitro-
    Entries: EA273701, EA273702, EA273703, EA273704, EA273705, EA273706, EA273707, EA273708, EA273709, EA273710, EA273711, EA273712, EA273713, EA273714, EA273751, EA273755, EA273758, EA273759
    Molecule: Mephedrone 2-(Methylamino)-1-(4-methylphenyl)-1-propanone 4-methylmethcathinone
    Entries: EA282701, EA282702, EA282703, EA282704, EA282705, EA282706, EA282707, EA282708, EA282709, EA282710, EA282711, EA282712, EA282713, EA282714
    Molecule: Galaxolidone 1,3,4,6,7,8-Hexahydro-4,6,6,7,8,8-hexamethylcyclopenta[g]-2-benzopyran-1-on
    Entries: EA066001, EA066002, EA066003, EA066004, EA066005, EA066006, EA066007, EA066008, EA066009, EA066010, EA066011, EA066012, EA066013, EA066014
    Molecule: Exemestane (8R,9S,10R,13S,14S)-10,13-dimethyl-6-methylene-7,8,9,11,12,14,15,16-octahydrocyclopenta[a]phenanthrene-3,17-dione
    Entries: EA066101, EA066102, EA066103, EA066104, EA066105, EA066106, EA066107, EA066108, EA066109, EA066110, EA066111, EA066112, EA066113, EA066114
    Molecule: Pymetrozine 6-methyl-4-(3-pyridinylmethylideneamino)-2,5-dihydro-1,2,4-triazin-3-one
    Entries: EA294701, EA294702, EA294703, EA294704, EA294705, EA294706, EA294707, EA294708, EA294709, EA294710, EA294711, EA294712, EA294713, EA294714
    Molecule: Iopromide 1-N,3-N-bis(2,3-dihydroxypropyl)-2,4,6-triiodo-5-[(2-methoxyacetyl)amino]-3-N-methylbenzene-1,3-dicarboxamide
    Entries: EA024201, EA024202, EA024203, EA024204, EA024205, EA024206, EA024207, EA024208, EA024209, EA024210, EA024211, EA024212, EA024214
    Molecule: Iprovalicarb Isopropyl (3-methyl-1-{[(1S)-1-(4-methylphenyl)ethyl]amino}-1-oxo-2-butanyl)carbamate
    Entries: EA293601, EA293602, EA293603, EA293604, EA293605, EA293606, EA293607, EA293608, EA293609, EA293610, EA293611, EA293612, EA293613, EA293614
    Molecule: Pyraclostrobin methyl [2-({[1-(4-chlorophenyl)-1H-pyrazol-3-yl]oxy}methyl)phenyl]methoxycarbamate
    Entries: EA277901, EA277902, EA277903, EA277904, EA277905, EA277906, EA277907, EA277908, EA277909, EA277910, EA277911, EA277912, EA277913, EA277914
    Molecule: Venlafaxine 1-[2-(dimethylamino)-1-(4-methoxyphenyl)ethyl]-1-cyclohexanol
    Entries: EA064501, EA064502, EA064503, EA064504, EA064505, EA064506, EA064507, EA064508, EA064509, EA064510, EA064511, EA064512, EA064513, EA064514
    Molecule: Spironolactone S-[(7R,8R,9S,10R,13S,14S,17R)-10,13-dimethyl-3,5'-bis(oxidanylidene)spiro[2,6,7,8,9,11,12,14,15,16-decahydro-1H-cyclopenta[a]phenanthrene-17,2'-oxolane]-7-yl] ethanethioate
    Entries: EA290202, EA290208
    Molecule: Sulfadimethoxine 4-amino-N-(2,6-dimethoxy-4-pyrimidinyl)benzenesulfonamide
    Entries: EA018301, EA018302, EA018303, EA018304, EA018305, EA018306, EA018307, EA018308, EA018309, EA018310, EA018311, EA018312, EA018313, EA018314, EA018351, EA018352, EA018353, EA018354, EA018355, EA018356, EA018357, EA018358, EA018359, EA018360, EA018361, EA018362, EA018363, EA018364
    Molecule: 4-Formylaminoantipyrine N-(1,5-dimethyl-3-oxidanylidene-2-phenyl-pyrazol-4-yl)methanamide
    Entries: EA103801, EA103802, EA103803, EA103804, EA103805, EA103806, EA103807, EA103808, EA103809, EA103810, EA103811, EA103812, EA103813, EA103814
    Molecule: Irgarol 2-N-tert-butyl-4-N-cyclopropyl-6-methylsulfanyl-1,3,5-triazine-2,4-diamine
    Entries: EA030201, EA030202, EA030203, EA030204, EA030205, EA030206, EA030207, EA030208, EA030209, EA030210, EA030211, EA030212, EA030213, EA030214
    Molecule: Propachlor ESA 2-[Isopropyl(phenyl)amino]-2-oxoethanesulfonic acid
    Entries: EA066501, EA066502, EA066503, EA066504, EA066505, EA066506, EA066507, EA066508, EA066509, EA066510, EA066511, EA066512, EA066513, EA066514, EA066551, EA066552, EA066553, EA066554, EA066555, EA066556, EA066557, EA066558, EA066559, EA066560, EA066561, EA066562, EA066563, EA066564
    Molecule: Metolachlor morpholinone 4-(2-ethyl-6-methylphenyl)-5-methylmorpholin-3-one
    Entries: EA006701, EA006702, EA006703, EA006704, EA006705, EA006706, EA006707, EA006708, EA006709, EA006710, EA006711, EA006712, EA006713, EA006714
    Molecule: Bisperfluorooctyl phosphate 6:2diPAP Bis(1H,1H,2H,2H-perfluorooctyl)phosphate
    Entries: EA292201, EA292202, EA292203, EA292204, EA292205, EA292206, EA292207, EA292208, EA292209, EA292210, EA292211, EA292212, EA292213, EA292214
    Molecule: Bromacil Bromazil 5-bromanyl-3-butan-2-yl-6-methyl-1H-pyrimidine-2,4-dione
    Entries: EA026601, EA026602, EA026603, EA026604, EA026605, EA026606, EA026607, EA026608, EA026609, EA026610, EA026611, EA026612, EA026613, EA026614, EA026651, EA026652, EA026653, EA026654, EA026655, EA026656, EA026657, EA026658, EA026659, EA026660, EA026661, EA026662, EA026663, EA026664
    Molecule: Dronedarone N-[2-butyl-3-[4-[3-(dibutylamino)propoxy]benzoyl]-1-benzofuran-5-yl]methanesulfonamide
    Entries: EA285501, EA285502, EA285503, EA285504, EA285505, EA285506, EA285507, EA285508, EA285509, EA285510, EA285511, EA285512, EA285513, EA285514
    Molecule: Azoxystrobin (free acid)
    Entries: EA273401, EA273402, EA273403, EA273404, EA273405, EA273406, EA273407, EA273408, EA273409, EA273410, EA273411, EA273412, EA273413, EA273414, EA273451, EA273452, EA273453, EA273458, EA273459
    Molecule: Orbencarb N,N-diethylcarbamothioic acid S-[(2-chlorophenyl)methyl] ester
    Entries: EA006901, EA006902, EA006903, EA006904, EA006905, EA006906, EA006907, EA006908, EA006909, EA006910, EA006911, EA006912, EA006913, EA006914
    Molecule: Dimethenamid-P 2-chloranyl-N-(2,4-dimethylthiophen-3-yl)-N-[(2S)-1-methoxypropan-2-yl]ethanamide
    Entries: EA025401, EA025402, EA025403, EA025404, EA025405, EA025406, EA025407, EA025408, EA025409, EA025410, EA025411, EA025412, EA025413, EA025414
    Molecule: Dioxoaminopyrine AMDOPH 2-(N-[acetyl(methyl)amino]anilino)-2-keto-N,N-dimethyl-acetamide
    Entries: EA270801, EA270802, EA270803, EA270804, EA270805, EA270806, EA270807, EA270808, EA270809, EA270810, EA270811, EA270812, EA270813, EA270814
    Molecule: Terbutylazine 2-N-tert-butyl-6-chloro-4-N-ethyl-1,3,5-triazine-2,4-diamine
    Entries: EA028401, EA028402, EA028403, EA028404, EA028405, EA028406, EA028407, EA028408, EA028409, EA028410, EA028411, EA028412, EA028413, EA028414
    Molecule: Ethofumesate (2-ethoxy-3,3-dimethyl-2H-1-benzofuran-5-yl) methanesulfonate
    Entries: EA011701, EA011702, EA011703, EA011704, EA011705, EA011706, EA011707, EA011708, EA011709, EA011710, EA011711, EA011712, EA011713, EA011714
    Molecule: Flufenacet ESA 2-[(4-Fluorophenyl)(isopropyl)amino]-2-oxoethanesulfonic acid
    Entries: EA066301, EA066302, EA066303, EA066304, EA066305, EA066306, EA066307, EA066308, EA066309, EA066310, EA066311, EA066312, EA066313, EA066314, EA066351, EA066352, EA066353, EA066354, EA066355, EA066356, EA066357, EA066358, EA066359, EA066360, EA066361, EA066362, EA066363, EA066364
    Molecule: N4-Acetylsulfamethoxazole N-[4-[(5-methyl-1,2-oxazol-3-yl)sulfamoyl]phenyl]acetamide
    Entries: EA029901, EA029902, EA029903, EA029904, EA029905, EA029906, EA029907, EA029908, EA029909, EA029910, EA029911, EA029912, EA029913, EA029914, EA029951, EA029952, EA029953, EA029954, EA029955, EA029956, EA029957, EA029958, EA029959, EA029960, EA029961, EA029962
    Molecule: N,N-Dimethyl-N'-phenylsulfamide DMSA
    Entries: EA029601, EA029602, EA029603, EA029604, EA029605, EA029606, EA029607, EA029608, EA029609, EA029610, EA029611, EA029612, EA029613, EA029614, EA029651, EA029652, EA029653, EA029654, EA029655, EA029656, EA029657, EA029658, EA029659, EA029660, EA029661, EA029662, EA029664
    Molecule: Sulfadiazine 4-amino-N-(2-pyrimidinyl)benzenesulfonamide
    Entries: EA017901, EA017902, EA017903, EA017904, EA017905, EA017906, EA017907, EA017908, EA017909, EA017910, EA017911, EA017912, EA017913, EA017914
    Molecule: Furosemide Frusemide 4-chloranyl-2-(furan-2-ylmethylamino)-5-sulfamoyl-benzoic acid
    Entries: EA260001, EA260002, EA260003, EA260004, EA260005, EA260006, EA260007, EA260008, EA260009, EA260010, EA260011, EA260012, EA260013, EA260014, EA260051, EA260052, EA260053, EA260054, EA260055, EA260056, EA260057, EA260058, EA260059, EA260060, EA260061, EA260062, EA260063, EA260064
    Molecule: 4-Aminoantipyrine 4-amino-1,5-dimethyl-2-phenyl-3-pyrazolin-3-one
    Entries: EA084501, EA084502, EA084503, EA084504, EA084505, EA084506, EA084507, EA084508, EA084509, EA084510, EA084511, EA084512, EA084513, EA084514
    Molecule: Iomeprol 1-N,3-N-bis(2,3-dihydroxypropyl)-5-[(2-hydroxyacetyl)-methylamino]-2,4,6-triiodobenzene-1,3-dicarboxamide
    Entries: EA024001, EA024002, EA024003, EA024004, EA024005, EA024006, EA024007, EA024008, EA024009, EA024010, EA024011, EA024012, EA024013, EA024014
    Molecule: Chlorpyrifos-methyl dimethoxy-sulfanylidene-(3,5,6-trichloropyridin-2-yl)oxy-{5}-phosphane
    Entries: EA295101, EA295102, EA295103, EA295104, EA295105, EA295106, EA295107, EA295108, EA295109, EA295110, EA295111, EA295112, EA295113, EA295114
    Molecule: Chloridazone-desphenyl 3(2H)-Pyridazinone, 5-amino-4-chloro- 4-amino-5-chloro-1H-pyridazin-6-one
    Entries: EA069801, EA069802, EA069803, EA069804, EA069805, EA069806, EA069807, EA069808, EA069809, EA069810, EA069811, EA069812, EA069813, EA069814, EA069852, EA069853, EA069854, EA069855, EA069858, EA069859, EA069860, EA069861
    Molecule: Fluconazole 2-(2,4-difluorophenyl)-1,3-bis(1,2,4-triazol-1-yl)-2-propanol
    Entries: EA032801, EA032802, EA032803, EA032804, EA032805, EA032806, EA032807, EA032808, EA032809, EA032810, EA032811, EA032812, EA032813, EA032814, EA032851, EA032852, EA032853, EA032854, EA032855, EA032856, EA032857, EA032858, EA032859, EA032860, EA032861, EA032862, EA032863, EA032864
    Molecule: Foramsulfuron 2-[(4,6-dimethoxypyrimidin-2-yl)carbamoylsulfamoyl]-4-formamido-N,N-dimethyl-benzamide
    Entries: EA277801, EA277802, EA277803, EA277804, EA277805, EA277806, EA277807, EA277808, EA277809, EA277810, EA277811, EA277812, EA277813, EA277814, EA277851, EA277852, EA277853, EA277854, EA277855, EA277856, EA277857, EA277858, EA277859, EA277860, EA277861, EA277863, EA277864
    Molecule: Dinoseb 2,4-dinitro-6-sec-butyl-phenol
    Entries: EA025751, EA025752, EA025753, EA025754, EA025755, EA025756, EA025757, EA025758, EA025759, EA025760, EA025761, EA025762, EA025763, EA025764
    Molecule: Ethofumesate-2-keto 3,3-Dimethyl-2-oxo-2,3-dihydrobenzofurane-5-yl methane sulfonate
    Entries: EA011801, EA011802, EA011803, EA011804, EA011805, EA011806, EA011807, EA011808, EA011809, EA011810, EA011811, EA011812, EA011813, EA011814
    Molecule: Metazachlor OXA 2-[2,6-dimethyl-N-(pyrazol-1-ylmethyl)anilino]-2-oxo-acetic acid
    Entries: EA070601, EA070602, EA070603, EA070604, EA070605, EA070606, EA070607, EA070608, EA070609, EA070610, EA070611, EA070612, EA070613, EA070614, EA070652, EA070653, EA070654, EA070655, EA070656, EA070657, EA070658, EA070659, EA070660, EA070661, EA070662, EA070663
    Molecule: Metamitron-desamino 3-Methyl-6-phenyl-1,2,4-triazin-5-ol
    Entries: EA000401, EA000402, EA000403, EA000404, EA000405, EA000406, EA000407, EA000408, EA000409, EA000410, EA000411, EA000412, EA000413, EA000414, EA000451, EA000452, EA000453, EA000454, EA000455, EA000456, EA000457, EA000458, EA000459, EA000460, EA000461, EA000462, EA000463, EA000464
    Molecule: Prometon 6-methoxy-2-N,4-N-di(propan-2-yl)-1,3,5-triazine-2,4-diamine
    Entries: EA013201, EA013202, EA013203, EA013204, EA013205, EA013206, EA013207, EA013208, EA013209, EA013210, EA013211, EA013212, EA013213, EA013214
    Molecule: 3-Phenoxybenzoic acid
    Entries: EA032601, EA032609, EA032614, EA032651, EA032652, EA032653, EA032654, EA032655, EA032656, EA032657, EA032658, EA032659, EA032660, EA032661, EA032662, EA032663, EA032664
    Molecule: DEET N,N-Diethyl-m-toluamide N,N-diethyl-3-methyl-benzamide
    Entries: EA021301, EA021302, EA021303, EA021304, EA021305, EA021306, EA021307, EA021308, EA021309, EA021310, EA021311, EA021312, EA021313, EA021314
    Molecule: Primidone Primaclone 5-ethyl-5-phenyl-1,3-diazinane-4,6-dione
    Entries: EA019501, EA019502, EA019503, EA019504, EA019505, EA019506, EA019507, EA019508, EA019509, EA019510, EA019511, EA019512, EA019513, EA019514
    Molecule: 4-Trifluoromethylphenol 4-(trifluoromethyl)phenol
    Entries: EA067351, EA067352, EA067353, EA067354, EA067355, EA067356, EA067357, EA067358, EA067359, EA067360, EA067361, EA067362, EA067363
    Molecule: Sulfamethoxazole 4-Amino-N-(5-methyl-3-isoxazolyl)benzenesulfonamide 4-amino-N-(5-methyl-1,2-oxazol-3-yl)benzenesulfonamide
    Entries: EA029801, EA029802, EA029803, EA029804, EA029805, EA029806, EA029807, EA029808, EA029809, EA029810, EA029811, EA029812, EA029813, EA029814, EA029851, EA029852, EA029853, EA029854, EA029855, EA029856, EA029858, EA029859, EA029860, EA029861, EA029864
    Molecule: Myclobutanil MYC 2-(4-chlorophenyl)-2-(1,2,4-triazol-1-ylmethyl)hexanenitrile
    Entries: EA295701, EA295702, EA295703, EA295704, EA295705, EA295706, EA295707, EA295708, EA295709, EA295710, EA295711, EA295712, EA295713, EA295714
    Molecule: Sitagliptin (3R)-3-amino-1-[3-(trifluoromethyl)-6,8-dihydro-5H-[1,2,4]triazolo[4,3-a]pyrazin-7-yl]-4-(2,4,5-trifluorophenyl)-1-butanone
    Entries: EA290301, EA290302, EA290303, EA290304, EA290305, EA290306, EA290307, EA290308, EA290309, EA290310, EA290311, EA290312, EA290313, EA290314
    Molecule: Hydrochlorothiazide 6-Chloro-3,4-dihydro-2H-1,2,4-benzothiadiazine-7-sulfonamide 1,1-dioxide 6-chloranyl-1,1-bis(oxidanylidene)-3,4-dihydro-2H-1$l^{6},2,4-benzothiadiazine-7-sulfonamide
    Entries: EA261001, EA261002, EA261003, EA261004, EA261005, EA261006, EA261007, EA261008, EA261009, EA261010, EA261011, EA261012, EA261013, EA261014, EA261051, EA261052, EA261053, EA261054, EA261055, EA261056, EA261057, EA261058, EA261059, EA261060, EA261061, EA261062, EA261063, EA261064
    Molecule: Cyprodinil (4-cyclopropyl-6-methyl-pyrimidin-2-yl)-phenyl-amine
    Entries: EA014801, EA014802, EA014803, EA014804, EA014805, EA014806, EA014807, EA014808, EA014809, EA014810, EA014811, EA014812, EA014813, EA014814
    Molecule: Pantoprazole 6-(difluoromethoxy)-2-[(3,4-dimethoxy-2-pyridinyl)methylsulfinyl]-1H-benzimidazole
    Entries: EA064401, EA064402, EA064403, EA064404, EA064405, EA064406, EA064407, EA064408, EA064409, EA064410, EA064411, EA064412, EA064413, EA064414, EA064451, EA064452, EA064453, EA064454, EA064455, EA064456, EA064457, EA064458, EA064459, EA064460, EA064461, EA064462, EA064463, EA064464
    Molecule: Methiocarb-sulfoxide Mesurol sulfoxide (3,5-dimethyl-4-methylsulfinyl-phenyl) N-methylcarbamate
    Entries: EA295301, EA295302, EA295303, EA295304, EA295305, EA295306, EA295307, EA295308, EA295309, EA295310, EA295311, EA295312, EA295313, EA295314
    Molecule: Pirimicarb N,N-dimethylcarbamic acid [2-(dimethylamino)-5,6-dimethyl-4-pyrimidinyl] ester
    Entries: EA271101, EA271102, EA271103, EA271104, EA271105, EA271106, EA271107, EA271108, EA271109, EA271110, EA271111, EA271112, EA271113, EA271114
    Molecule: Chloridazone-methyl-desphenyl 1-Methyl-4-amino-5-chloro-6-oxo-(1H)-pyridazine 5-amino-4-chloro-2-methyl-3-pyridazinone
    Entries: EA069901, EA069902, EA069903, EA069904, EA069905, EA069906, EA069907, EA069908, EA069909, EA069910, EA069911, EA069912, EA069913, EA069914, EA069952, EA069953, EA069959
    Molecule: Valsartan (2S)-3-methyl-2-[1-oxopentyl-[[4-[2-(2H-tetrazol-5-yl)phenyl]phenyl]methyl]amino]butanoic acid
    Entries: EA258301, EA258302, EA258303, EA258304, EA258305, EA258306, EA258307, EA258308, EA258309, EA258310, EA258311, EA258312, EA258313, EA258314, EA258351, EA258352, EA258353, EA258354, EA258355, EA258356, EA258358, EA258359, EA258360, EA258361, EA258362, EA258364
    Molecule: Sulfamethazine 4-amino-N-(4,6-dimethylpyrimidin-2-yl)benzenesulfonamide 4-amino-N-(4,6-dimethyl-2-pyrimidinyl)benzenesulfonamide
    Entries: EA018101, EA018102, EA018103, EA018104, EA018105, EA018106, EA018107, EA018108, EA018109, EA018110, EA018111, EA018112, EA018113, EA018114
    Molecule: EDDP (2E)-2-ethylidene-1,5-dimethyl-3,3-diphenyl-pyrrolidine
    Entries: EA282501, EA282502, EA282503, EA282504, EA282505, EA282506, EA282507, EA282508, EA282509, EA282510, EA282511, EA282512, EA282513, EA282514
    Molecule: Ibuprofen (S)-(+)-Ibuprofen (2S)-2-(4-isobutylphenyl)propionic acid
    Entries: EA020301, EA020302, EA020303, EA020304, EA020305, EA020306, EA020307, EA020308, EA020309, EA020310, EA020311, EA020312, EA020314, EA020351, EA020352, EA020358, EA020364
    Molecule: Dimethenamide ESA 2-[(2,4-Dimethyl-3-thienyl)(1-methoxy-2-propanyl)amino]-2-oxoethanesulfonic acid
    Entries: EA026001, EA026002, EA026003, EA026004, EA026005, EA026006, EA026007, EA026008, EA026009, EA026010, EA026011, EA026012, EA026013, EA026014, EA026051, EA026052, EA026053, EA026054, EA026055, EA026056, EA026057, EA026058, EA026059, EA026060, EA026061, EA026062, EA026063, EA026064
    Molecule: Terbutylazine-2-hydroxy 4-(Ethylamino)-6-[(2-methyl-2-propanyl)amino]-1,3,5-triazin-2(5H)-one
    Entries: EA034701, EA034702, EA034703, EA034704, EA034705, EA034706, EA034707, EA034708, EA034709, EA034710, EA034711, EA034712, EA034713, EA034714, EA034751, EA034752, EA034753, EA034754, EA034755, EA034756, EA034757, EA034758, EA034759, EA034760, EA034761, EA034762, EA034763
    Molecule: 4,5-Dichloro-2-n-octyl-3(2H)-isothiazolone (DCOIT) 3(2H)-Isothiazolone, 4,5-dichloro-2-octyl- 4,5-bis(chloranyl)-2-octyl-1,2-thiazol-3-one
    Entries: EA034401, EA034402, EA034403, EA034404, EA034405, EA034406, EA034407, EA034408, EA034409, EA034410, EA034411, EA034412, EA034413, EA034414
    Molecule: Atorvastatin (3R,5R)-7-[3-(anilinocarbonyl)-5-(4-fluorophenyl)-2-isopropyl-4-phenyl-1H-pyrrol-1-yl]-3,5-dihydroxyheptanoic acid (3R,5R)-7-[2-(4-fluorophenyl)-3-phenyl-4-(phenylcarbamoyl)-5-propan-2-yl-pyrrol-1-yl]-3,5-bis(oxidanyl)heptanoic acid
    Entries: EA281001, EA281002, EA281003, EA281004, EA281005, EA281006, EA281007, EA281008, EA281009, EA281010, EA281011, EA281012, EA281013, EA281014, EA281051, EA281052, EA281053, EA281054, EA281055, EA281056, EA281057, EA281058, EA281059, EA281060, EA281061, EA281062, EA281063, EA281064
    Molecule: O-desmethylvenlafaxine Desvenlafaxine 4-[2-(dimethylamino)-1-(1-hydroxycyclohexyl)ethyl]phenol
    Entries: EA105301, EA105302, EA105303, EA105304, EA105305, EA105306, EA105307, EA105308, EA105309, EA105310, EA105311, EA105312, EA105313, EA105314
    Molecule: Triclocarban 1-(4-chlorophenyl)-3-(3,4-dichlorophenyl)urea
    Entries: EA298801, EA298802, EA298803, EA298804, EA298805, EA298806, EA298807, EA298808, EA298809, EA298810, EA298811, EA298812, EA298813, EA298814, EA298851, EA298859
    Molecule: Fenoxycarb N-[2-(4-phenoxyphenoxy)ethyl]carbamic acid ethyl ester
    Entries: EA294001, EA294002, EA294003, EA294004, EA294005, EA294006, EA294007, EA294008, EA294009, EA294010, EA294011, EA294012, EA294013, EA294014
    Molecule: Metolachlor ESA 2-[(2-Ethyl-6-methylphenyl)(1-methoxy-2-propanyl)amino]-2-oxoethanesulfonic acid
    Entries: EA050201, EA050202, EA050203, EA050204, EA050205, EA050206, EA050207, EA050208, EA050209, EA050210, EA050211, EA050212, EA050213, EA050214, EA050251, EA050252, EA050253, EA050254, EA050255, EA050256, EA050257, EA050258, EA050259, EA050260, EA050261, EA050262, EA050263, EA050264
    Molecule: Napropamid N,N-diethyl-2-(1-naphthalenyloxy)propanamide
    Entries: EA012001, EA012002, EA012003, EA012004, EA012005, EA012006, EA012007, EA012008, EA012009, EA012010, EA012011, EA012012, EA012013, EA012014
    Molecule: Iminostilbene 11H-benzo[b][1]benzazepine
    Entries: EA099201, EA099202, EA099203, EA099204, EA099205, EA099206, EA099207, EA099208, EA099209, EA099210, EA099211, EA099212, EA099213, EA099214
    Molecule: N-Desmethylvenlafaxine 1-[1-(4-methoxyphenyl)-2-(methylamino)ethyl]-1-cyclohexanol
    Entries: EA103401, EA103402, EA103403, EA103404, EA103405, EA103406, EA103407, EA103408, EA103409, EA103410, EA103411, EA103412, EA103413, EA103414
    Molecule: Candesartan 2-ethoxy-3-[4-[2-(2H-tetrazol-5-yl)phenyl]benzyl]benzimidazole-4-carboxylic acid
    Entries: EA280401, EA280402, EA280403, EA280404, EA280405, EA280406, EA280407, EA280408, EA280409, EA280410, EA280411, EA280412, EA280413, EA280414, EA280451, EA280452, EA280453, EA280454, EA280455, EA280456, EA280457, EA280458, EA280459, EA280460, EA280461, EA280462, EA280463, EA280464
    Molecule: Alachlor 2-chloranyl-N-(2,6-diethylphenyl)-N-(methoxymethyl)ethanamide
    Entries: EA027401, EA027402, EA027403, EA027404, EA027405, EA027406, EA027407, EA027408, EA027409, EA027410, EA027411, EA027412, EA027413, EA027414
    Molecule: 3,5-Dibromo-4-hydroxybenzoic acid 3,5-bis(bromanyl)-4-oxidanyl-benzoic acid
    Entries: EA080451, EA080452, EA080453, EA080454, EA080455, EA080456, EA080457, EA080458, EA080459, EA080460, EA080461, EA080462, EA080463, EA080464
    Molecule: 2-Octyl-3(2H)-isothiazolone 2-octyl-1,2-thiazol-3-one
    Entries: EA030001, EA030002, EA030003, EA030004, EA030005, EA030006, EA030007, EA030008, EA030009, EA030010, EA030011, EA030012, EA030013, EA030014
    Molecule: Pethoxamide 2-chloro-N-(2-ethoxyethyl)-N-(2-methyl-1-phenylprop-1-enyl)acetamide
    Entries: EA070001, EA070002, EA070003, EA070004, EA070005, EA070006, EA070007, EA070008, EA070009, EA070010, EA070011, EA070012, EA070013, EA070014
    Molecule: Isoproturon-monodemethyl 1-(4-Isoprophenyl)-3-methylurea
    Entries: EA030401, EA030402, EA030403, EA030404, EA030405, EA030406, EA030407, EA030408, EA030409, EA030410, EA030411, EA030412, EA030413, EA030414, EA030451, EA030452, EA030453, EA030454, EA030455, EA030459, EA030460, EA030464
    Molecule: Climbazol Climbazole 1-(4-chloranylphenoxy)-1-imidazol-1-yl-3,3-dimethyl-butan-2-one
    Entries: EA273901, EA273902, EA273903, EA273904, EA273905, EA273906, EA273907, EA273908, EA273909, EA273910, EA273911, EA273912, EA273913, EA273914, EA273951, EA273952, EA273953, EA273954, EA273955, EA273956, EA273957, EA273958, EA273959, EA273960, EA273961, EA273962, EA273964
    Molecule: Diuron 3-(3,4-dichlorophenyl)-1,1-dimethyl-urea
    Entries: EA029201, EA029202, EA029203, EA029204, EA029205, EA029206, EA029207, EA029208, EA029209, EA029210, EA029211, EA029212, EA029213, EA029214, EA029251, EA029252, EA029253, EA029254, EA029255, EA029256, EA029257, EA029258, EA029259, EA029260, EA029261, EA029262, EA029263, EA029264
    Molecule: 2-Methyl-4-amino-6-methoxy-s-triazine (4-methoxy-6-methyl-s-triazin-2-yl)amine
    Entries: EA012501, EA012502, EA012503, EA012504, EA012505, EA012506, EA012507, EA012508, EA012509, EA012510, EA012511, EA012512, EA012513, EA012514
    Molecule: Bromoxynil 3,5-bis(bromanyl)-4-oxidanyl-benzenecarbonitrile
    Entries: EA002451, EA002452, EA002453, EA002454, EA002455, EA002456, EA002457, EA002458, EA002459, EA002460, EA002461, EA002462, EA002463, EA002464
    Molecule: Sucralose 1,6-Dichloro-1,6-dideoxy-beta-D-fructofuranosyl 4-chloro-4-deoxy-alpha-D-galactopyranoside (2R,3R,4R,5R,6R)-2-[(2R,3S,4S,5S)-2,5-bis(chloromethyl)-3,4-bis(oxidanyl)oxolan-2-yl]oxy-5-chloranyl-6-(hydroxymethyl)oxane-3,4-diol
    Entries: EA070351, EA070352, EA070358, EA070364
    Molecule: Methadone 6-(dimethylamino)-4,4-diphenyl-3-heptanone
    Entries: EA282801, EA282802, EA282803, EA282804, EA282805, EA282806, EA282807, EA282808, EA282809, EA282810, EA282811, EA282812, EA282813, EA282814
    Molecule: Boscalid 2-chloranyl-N-[2-(4-chlorophenyl)phenyl]pyridine-3-carboxamide
    Entries: EA293801, EA293802, EA293803, EA293804, EA293805, EA293806, EA293807, EA293808, EA293809, EA293810, EA293811, EA293812, EA293813, EA293814, EA293852, EA293858
    Molecule: Perfluorobutyric acid Heptafluorobutyric acid 2,2,3,3,4,4,4-heptafluorobutanoic acid
    Entries: EA271351, EA271352, EA271353, EA271357, EA271358, EA271359, EA271360
    Molecule: Fenofibric acid 2-[4-(4-chlorobenzoyl)phenoxy]-2-methyl-propionic acid
    Entries: EA033701, EA033702, EA033703, EA033704, EA033705, EA033706, EA033707, EA033708, EA033709, EA033710, EA033711, EA033712, EA033713, EA033714, EA033751, EA033752, EA033753, EA033754, EA033755, EA033756, EA033758, EA033759, EA033760, EA033761, EA033762, EA033763
    Molecule: Mecoprop 2-(4-chloro-2-methyl-phenoxy)propionic acid
    Entries: EA030851, EA030852, EA030853, EA030854, EA030855, EA030856, EA030857, EA030858, EA030859, EA030860, EA030861, EA030862, EA030863, EA030864
    Molecule: Secbumeton 2-N-butan-2-yl-4-N-ethyl-6-methoxy-1,3,5-triazine-2,4-diamine
    Entries: EA066901, EA066902, EA066903, EA066904, EA066905, EA066906, EA066907, EA066908, EA066909, EA066910, EA066911, EA066912, EA066913, EA066914
    Molecule: Tebufenozide N-tert-butyl-N'-(4-ethylbenzoyl)-3,5-dimethyl-benzohydrazide
    Entries: EA295203, EA295204, EA295205, EA295206, EA295208, EA295209, EA295210, EA295211, EA295212, EA295213, EA295214, EA295251, EA295252, EA295253, EA295254, EA295255, EA295256, EA295258, EA295259, EA295260, EA295261
    Molecule: Thifensulfuron-methyl 3-[(4-methoxy-6-methyl-s-triazin-2-yl)carbamoylsulfamoyl]thiophene-2-carboxylic acid methyl ester
    Entries: EA012401, EA012402, EA012403, EA012404, EA012405, EA012406, EA012407, EA012408, EA012409, EA012410, EA012411, EA012412, EA012413, EA012414, EA012451, EA012452, EA012453, EA012454, EA012455, EA012456, EA012457, EA012458, EA012459, EA012460, EA012461, EA012462, EA012463, EA012464
    Molecule: Prednisolone 1-Dehydrohydrocortisone (8S,9S,10R,11S,13S,14S,17R)-10,13-dimethyl-11,17-bis(oxidanyl)-17-(2-oxidanylethanoyl)-7,8,9,11,12,14,15,16-octahydro-6H-cyclopenta[a]phenanthren-3-one
    Entries: EA278301, EA278302, EA278303, EA278304, EA278305, EA278306, EA278307, EA278308, EA278309, EA278310, EA278311, EA278312, EA278313, EA278314
    Molecule: Paracetamol Acetaminophen N-(4-hydroxyphenyl)acetamide
    Entries: EA024301, EA024302, EA024303, EA024304, EA024305, EA024306, EA024307, EA024308, EA024309, EA024310, EA024311, EA024312, EA024313, EA024314, EA024351, EA024352, EA024353, EA024354, EA024355, EA024356, EA024357, EA024358, EA024359, EA024360, EA024361, EA024362, EA024364
    Molecule: Azoxystrobin methyl (2E)-2-(2-{[6-(2-cyanophenoxy)pyrimidin-4-yl]oxy}phenyl)-3-methoxyacrylate
    Entries: EA008901, EA008902, EA008903, EA008904, EA008905, EA008906, EA008907, EA008908, EA008909, EA008910, EA008911, EA008912, EA008913, EA008914
    Molecule: 5-Fluorouracil 5-fluoranyl-1H-pyrimidine-2,4-dione
    Entries: EA256601, EA256602, EA256603, EA256604, EA256605, EA256606, EA256607, EA256608, EA256609, EA256610, EA256611, EA256612, EA256613, EA256614, EA256652, EA256653, EA256654, EA256655, EA256656, EA256657, EA256658, EA256659, EA256660, EA256661, EA256662, EA256663
    Molecule: Imidacloprid-urea 1-[(6-Chloropyridin-3-yl)methyl]imidazolidin-2-one
    Entries: EA295601, EA295602, EA295603, EA295604, EA295605, EA295606, EA295607, EA295608, EA295609, EA295610, EA295611, EA295612, EA295613, EA295614
    Molecule: Rimsulfuron 1-(4,6-dimethoxy-2-pyrimidinyl)-3-[(3-ethylsulfonyl-2-pyridinyl)sulfonyl]urea
    Entries: EA013001, EA013002, EA013003, EA013004, EA013005, EA013006, EA013007, EA013008, EA013009, EA013010, EA013011, EA013012, EA013013, EA013014, EA013051, EA013052, EA013053, EA013054, EA013055, EA013056, EA013057, EA013058, EA013059, EA013060, EA013061, EA013062, EA013063, EA013064
    Molecule: Ifosfamide Isophosphamide N,3-Bis(2-chloroethyl)-1,3,2-oxazaphosphinan-2-amine 2-oxide
    Entries: EA268301, EA268302, EA268303, EA268304, EA268305, EA268306, EA268307, EA268308, EA268309, EA268310, EA268311, EA268312, EA268313, EA268314
    Molecule: Monuron 3-(4-chlorophenyl)-1,1-dimethyl-urea
    Entries: EA016101, EA016102, EA016103, EA016104, EA016105, EA016106, EA016107, EA016108, EA016109, EA016110, EA016111, EA016112, EA016113, EA016114, EA016151, EA016152, EA016153, EA016154, EA016155, EA016158, EA016159, EA016160, EA016161, EA016164
    Molecule: Lenacil 3-cyclohexyl-1,5,6,7-tetrahydrocyclopenta[d]pyrimidine-2,4-dione
    Entries: EA294901, EA294902, EA294903, EA294904, EA294905, EA294906, EA294907, EA294908, EA294909, EA294910, EA294911, EA294912, EA294913, EA294914
    Molecule: Bisphenol A 4-[1-(4-hydroxyphenyl)-1-methyl-ethyl]phenol
    Entries: EA016301, EA016303, EA016308, EA016309, EA016310, EA016312, EA016314
    Molecule: Clopidogrel carboxylic acid 2-(2-chlorophenyl)-2-(6,7-dihydro-4H-thieno[3,2-c]pyridin-5-yl)acetic acid
    Entries: EA285401, EA285402, EA285403, EA285404, EA285405, EA285406, EA285407, EA285408, EA285409, EA285410, EA285411, EA285412, EA285413, EA285414, EA285451, EA285464
    Molecule: Metamitron 4-Amino-4,5-dihydro-3-methyl-6-phenyl-1,2,4-triazin-5-one 4-amino-3-methyl-6-phenyl-1,2,4-triazin-5-one
    Entries: EA005801, EA005802, EA005803, EA005804, EA005805, EA005806, EA005807, EA005808, EA005809, EA005810, EA005811, EA005812, EA005813, EA005814, EA005851, EA005852, EA005854, EA005855, EA005857, EA005858, EA005862, EA005864
    Molecule: Flufenacet OXA [(4-Fluorophenyl)(isopropyl)amino](oxo)acetic acid
    Entries: EA066401, EA066402, EA066403, EA066404, EA066405, EA066406, EA066407, EA066408, EA066409, EA066410, EA066411, EA066412, EA066413, EA066414, EA066451, EA066452, EA066453, EA066454, EA066455, EA066456, EA066457, EA066458, EA066459, EA066460, EA066461, EA066462, EA066463, EA066464
    Molecule: Simeton 2-N,4-N-diethyl-6-methoxy-1,3,5-triazine-2,4-diamine
    Entries: EA066701, EA066702, EA066703, EA066704, EA066705, EA066706, EA066707, EA066708, EA066709, EA066710, EA066711, EA066712, EA066713, EA066714
    Molecule: Fludioxonil 4-(2,2-difluoro-1,3-benzodioxol-4-yl)-1H-pyrrole-3-carbonitrile
    Entries: EA016251, EA016252, EA016253, EA016254, EA016255, EA016256, EA016257, EA016258, EA016259, EA016260, EA016261, EA016262, EA016263, EA016264
    Molecule: Cocaine (1S,3S,4R,5R)-3-benzoyloxy-8-methyl-8-azabicyclo[3.2.1]octane-4-carboxylic acid methyl ester
    Entries: EA281701, EA281702, EA281703, EA281704, EA281705, EA281706, EA281707, EA281708, EA281709, EA281710, EA281711, EA281712, EA281713, EA281714
    Molecule: Aspartame 3-amino-4-[(1-benzyl-2-keto-2-methoxy-ethyl)amino]-4-keto-butyric acid
    Entries: EA277001, EA277002, EA277003, EA277004, EA277005, EA277006, EA277007, EA277008, EA277009, EA277010, EA277011, EA277012, EA277013, EA277014, EA277051, EA277052, EA277053, EA277054, EA277055, EA277056, EA277057, EA277058, EA277059, EA277060, EA277061, EA277062, EA277064
    Molecule: 2-Aminosulfonyl-benzoic acid methyl ester Methyl 2-(aminosulfonyl)benzoate 2-sulfamoylbenzoic acid methyl ester
    Entries: EA012601, EA012602, EA012603, EA012604, EA012605, EA012606, EA012607, EA012608, EA012609, EA012610, EA012611, EA012612, EA012613, EA012614
    Molecule: 2',2'-Difluoro-2'-deoxyuridine
    Entries: EA264901, EA264902, EA264903, EA264904, EA264905, EA264906, EA264907, EA264908, EA264909, EA264910, EA264911, EA264912, EA264913, EA264914, EA264951, EA264955, EA264956, EA264959, EA264962, EA264963
    Molecule: Bezafibrate 2-[4-[2-[(4-chlorobenzoyl)amino]ethyl]phenoxy]-2-methyl-propionic acid
    Entries: EA020901, EA020902, EA020903, EA020904, EA020905, EA020906, EA020907, EA020908, EA020909, EA020910, EA020911, EA020912, EA020913, EA020914, EA020951, EA020952, EA020953, EA020954, EA020955, EA020956, EA020958, EA020959, EA020960, EA020961, EA020962, EA020964
    Molecule: 10,11-trans-Dihydroxy-10,11-dihydrocarbamazepine (5S,6S)-5,6-bis(oxidanyl)-5,6-dihydrobenzo[b][1]benzazepine-11-carboxamide
    Entries: EA270101, EA270102, EA270103, EA270104, EA270105, EA270106, EA270107, EA270108, EA270109, EA270110, EA270111, EA270112, EA270113, EA270114
    Molecule: Terbumeton 2-N-tert-butyl-4-N-ethyl-6-methoxy-1,3,5-triazine-2,4-diamine 1,3,5-Triazine-2,4-diamine, N-(1,1-dimethylethyl)-N'-ethyl-6-methoxy-
    Entries: EA034601, EA034602, EA034603, EA034604, EA034605, EA034606, EA034607, EA034608, EA034609, EA034610, EA034611, EA034612, EA034613, EA034614
    Molecule: Flonicamid N-(cyanomethyl)-4-(trifluoromethyl)-3-pyridinecarboxamide
    Entries: EA294301, EA294302, EA294303, EA294304, EA294305, EA294306, EA294307, EA294308, EA294309, EA294310, EA294311, EA294312, EA294313, EA294314, EA294352, EA294354, EA294358, EA294359, EA294361
    Molecule: DDAO N,N-dimethyl-1-decanamine oxide NN-Dimethyldicylamine N-oxide
    Entries: EA278401, EA278402, EA278403, EA278404, EA278405, EA278406, EA278407, EA278408, EA278409, EA278410, EA278411, EA278412, EA278413, EA278414
    Molecule: Ephedrine (1R,2S)-2-(methylamino)-1-phenyl-1-propanol
    Entries: EA275801, EA275802, EA275803, EA275804, EA275805, EA275806, EA275807, EA275808, EA275809, EA275810, EA275811, EA275812, EA275813, EA275814
    Molecule: Metribuzin-diketo 4-amino-6-tert-butyl-2H-1,2,4-triazine-3,5-dione
    Entries: EA009301, EA009302, EA009303, EA009304, EA009305, EA009306, EA009307, EA009308, EA009309, EA009310, EA009311, EA009312, EA009313, EA009314, EA009352, EA009353, EA009354, EA009358, EA009359
    Molecule: Epoxiconazole 1-[[(2S,3R)-3-(2-chlorophenyl)-2-(4-fluorophenyl)-2-oxiranyl]methyl]-1,2,4-triazole
    Entries: EA009501, EA009502, EA009503, EA009504, EA009505, EA009506, EA009507, EA009508, EA009509, EA009510, EA009511, EA009512, EA009513, EA009514
    Molecule: Tebutam 2,2-dimethyl-N-(phenylmethyl)-N-propan-2-yl-propanamide N-Benzyl-N-isopropylpivalamide
    Entries: EA025501, EA025502, EA025503, EA025504, EA025505, EA025506, EA025507, EA025508, EA025509, EA025510, EA025511, EA025512, EA025513, EA025514
    Molecule: Aminopyrine 4-(dimethylamino)-1,5-dimethyl-2-phenyl-3-pyrazolin-3-one
    Entries: EA070201, EA070202, EA070203, EA070204, EA070205, EA070206, EA070207, EA070208, EA070209, EA070210, EA070211, EA070212, EA070213, EA070214
    Molecule: Methiocarb (3,5-dimethyl-4-methylsulfanyl-phenyl) N-methylcarbamate
    Entries: EA293701, EA293702, EA293703, EA293704, EA293705, EA293706, EA293707, EA293708, EA293709, EA293710, EA293711, EA293712, EA293713, EA293714
    Molecule: Alachlor ESA 2-[(2,6-diethylphenyl)-(methoxymethyl)amino]-2-oxidanylidene-ethanesulfonic acid
    Entries: EA009951, EA009952, EA009953, EA009954, EA009955, EA009956, EA009957, EA009958, EA009959, EA009960, EA009961, EA009962, EA009963, EA009964
    Molecule: Fenofibrate isopropyl 2-[4-(4-chlorobenzoyl)phenoxy]-2-methylpropanoate 2-[4-(4-chlorobenzoyl)phenoxy]-2-methyl-propionic acid isopropyl ester
    Entries: EA033601, EA033602, EA033603, EA033604, EA033605, EA033606, EA033607, EA033608, EA033609, EA033610, EA033611, EA033612, EA033613, EA033614
    Molecule: Asulam N-(4-aminophenyl)sulfonylcarbamic acid methyl ester
    Entries: EA015601, EA015602, EA015603, EA015604, EA015605, EA015606, EA015607, EA015608, EA015609, EA015610, EA015611, EA015612, EA015613, EA015614, EA015651, EA015652, EA015653, EA015654, EA015655, EA015658, EA015659, EA015660, EA015661, EA015664
    Molecule: Levamisole (6S)-6-phenyl-2,3,5,6-tetrahydroimidazo[2,1-b][1,3]thiazole
    Entries: EA285701, EA285702, EA285703, EA285704, EA285705, EA285706, EA285707, EA285708, EA285709, EA285710, EA285711, EA285712, EA285713, EA285714
    Molecule: Perfluorooctane sulfonamidoacetic acid {[(heptadecafluorooctyl)sulfonyl]amino}acetic acid FOSAA
    Entries: EA291951, EA291954, EA291955, EA291958, EA291959, EA291960, EA291961, EA291962
    Molecule: Dimethoate 2-(dimethoxyphosphinothioylthio)-N-methylacetamide
    Entries: EA276101, EA276102, EA276103, EA276104, EA276105, EA276106, EA276107, EA276108, EA276109, EA276110, EA276111, EA276112, EA276113, EA276114
    Molecule: Ranitidine N-oxide 1,1-Ethenediamine, N-(2-(((5-((dimethylamino)methyl)-2-furanyl)methyl)thio)ethyl)-N'-methyl-2-nitro-, N-oxide
    Entries: EA273801, EA273802, EA273803, EA273804, EA273805, EA273806, EA273807, EA273808, EA273809, EA273810, EA273811, EA273812, EA273813, EA273814, EA273851, EA273853
    Molecule: Bupropion 2-(tert-butylamino)-1-(3-chlorophenyl)-1-propanone
    Entries: EA280301, EA280302, EA280303, EA280304, EA280305, EA280306, EA280307, EA280308, EA280309, EA280310, EA280311, EA280312, EA280313, EA280314
    Molecule: Amitriptyline 3-(10,11-Dihydro-5H-dibenzo[a,d][7]annulen-5-ylidene)-N,N-dimethyl-1-propanamine
    Entries: EA282101, EA282102, EA282103, EA282104, EA282105, EA282106, EA282107, EA282108, EA282109, EA282110, EA282111, EA282112, EA282113, EA282114
    Molecule: Telmisartan 4'-[(1,7'-dimethyl-2'-propyl-1H,3'H-2,5'-bibenzimidazol-3'-yl)methyl]biphenyl-2-carboxylic acid 2-[4-[[4-methyl-6-(1-methyl-2-benzimidazolyl)-2-propyl-1-benzimidazolyl]methyl]phenyl]benzoic acid
    Entries: EA280501, EA280502, EA280503, EA280504, EA280505, EA280506, EA280507, EA280508, EA280509, EA280510, EA280511, EA280512, EA280513, EA280514, EA280551, EA280552, EA280553, EA280554, EA280555, EA280556, EA280558, EA280559, EA280560, EA280561, EA280562, EA280563, EA280564
    Molecule: Iohexol 5-[acetyl(2,3-dihydroxypropyl)amino]-1-N,3-N-bis(2,3-dihydroxypropyl)-2,4,6-triiodobenzene-1,3-dicarboxamide
    Entries: EA023101, EA023102, EA023103, EA023104, EA023105, EA023106, EA023108, EA023109, EA023110, EA023111, EA023112, EA023114
    Molecule: Trifluralin 2,6-dinitro-N,N-dipropyl-4-(trifluoromethyl)aniline
    Entries: EA012302, EA012308
    Molecule: Propaquizafop 2-(propan-2-ylideneamino)oxyethyl 2-[4-(6-chloranylquinoxalin-2-yl)oxyphenoxy]propanoate
    Entries: EA012201, EA012202, EA012203, EA012204, EA012205, EA012206, EA012207, EA012208, EA012209, EA012210, EA012211, EA012212, EA012213, EA012214
    Molecule: Clozapine 8-chloro-11-(4-methylpiperazin-1-yl)-5H-dibenzo[b,e][1,4]diazepine
    Entries: EA284101, EA284102, EA284103, EA284104, EA284105, EA284106, EA284107, EA284108, EA284109, EA284110, EA284111, EA284112, EA284113, EA284114
    Molecule: Flusilazole bis(4-fluorophenyl)-methyl-(1,2,4-triazol-1-ylmethyl)silane
    Entries: EA009701, EA009702, EA009703, EA009704, EA009705, EA009706, EA009707, EA009708, EA009709, EA009710, EA009711, EA009712, EA009713, EA009714
    Molecule: Gabapentin 1-(Aminomethyl)cyclohexaneacetic acid 2-[1-(aminomethyl)cyclohexyl]acetic acid
    Entries: EA256101, EA256102, EA256103, EA256104, EA256105, EA256106, EA256107, EA256108, EA256109, EA256110, EA256111, EA256112, EA256113, EA256114
    Molecule: Cyproconazole 2-(4-chlorophenyl)-3-cyclopropyl-1-(1,2,4-triazol-1-yl)-2-butanol
    Entries: EA008601, EA008602, EA008603, EA008604, EA008605, EA008606, EA008607, EA008608, EA008609, EA008610, EA008611, EA008612, EA008613, EA008614
    Molecule: Metazachlor ESA 2-[(2,6-dimethylphenyl)(1H-pyrazol-1-ylmethyl)amino]-2-oxoethanesulfonic acid
    Entries: EA070551, EA070552, EA070553, EA070554, EA070555, EA070556, EA070557, EA070558, EA070559, EA070560, EA070561, EA070562, EA070563, EA070564
    Molecule: Diuron-desmethyl DCPMU 1-(3,4-dichlorophenyl)-3-methyl-urea
    Entries: EA029401, EA029402, EA029403, EA029404, EA029405, EA029406, EA029407, EA029408, EA029409, EA029410, EA029411, EA029412, EA029413, EA029414, EA029451, EA029452, EA029453, EA029454, EA029455, EA029456, EA029457, EA029458, EA029459, EA029460, EA029461, EA029462, EA029463, EA029464
    Molecule: Methoxyfenozide N'-tert-butyl-N'-(3,5-dimethylbenzoyl)-3-methoxy-2-methyl-benzohydrazide
    Entries: EA293501, EA293502, EA293503, EA293504, EA293505, EA293506, EA293507, EA293508, EA293509, EA293510, EA293511, EA293512, EA293513, EA293514, EA293551, EA293552, EA293553, EA293554, EA293555, EA293558, EA293559, EA293560, EA293561
    Molecule: Dichlorprop 2-(2,4-dichlorophenoxy)propanoic acid
    Entries: EA027051, EA027052, EA027053, EA027054, EA027055, EA027056, EA027057, EA027058, EA027059, EA027060, EA027061, EA027062, EA027063, EA027064
    Molecule: Lidocaine 2-(diethylamino)-N-(2,6-dimethylphenyl)acetamide
    Entries: EA257201, EA257202, EA257203, EA257204, EA257205, EA257206, EA257207, EA257208, EA257209, EA257210, EA257211, EA257212, EA257213, EA257214
    Molecule: Amphetamine Amfetamine (1-methyl-2-phenyl-ethyl)amine
    Entries: EA282201, EA282202, EA282203, EA282204, EA282205, EA282206, EA282207, EA282208, EA282209, EA282210, EA282211, EA282212, EA282213, EA282214
    Molecule: Irbesartan 8-butyl-7-[4-[2-(2H-tetrazol-5-yl)phenyl]benzyl]-7,9-diazaspiro[4.4]non-8-en-6-one
    Entries: EA277401, EA277402, EA277403, EA277404, EA277405, EA277406, EA277407, EA277408, EA277409, EA277410, EA277411, EA277412, EA277413, EA277414, EA277451, EA277452, EA277453, EA277454, EA277455, EA277456, EA277457, EA277458, EA277459, EA277460, EA277461, EA277462, EA277463, EA277464
    Molecule: Nicosulfuron 2-[(4,6-dimethoxypyrimidin-2-yl)carbamoylsulfamoyl]-N,N-dimethyl-nicotinamide
    Entries: EA012901, EA012902, EA012903, EA012904, EA012905, EA012906, EA012907, EA012908, EA012909, EA012910, EA012911, EA012912, EA012913, EA012914, EA012951, EA012952, EA012953, EA012954, EA012955, EA012956, EA012957, EA012958, EA012959, EA012960, EA012961, EA012962, EA012963, EA012964
    Molecule: Thiopental
    Entries: EA274251, EA274252, EA274253, EA274254, EA274255, EA274256, EA274257, EA274258, EA274259, EA274260, EA274261, EA274262, EA274263, EA274264
    Molecule: Cyclamate Cyclohexylsulfamic acid Cyclohexanesulfamic acid
    Entries: EA281301, EA281302, EA281303, EA281304, EA281305, EA281308, EA281309, EA281310, EA281314, EA281351, EA281352, EA281353, EA281354, EA281355, EA281356, EA281357, EA281358, EA281359, EA281360, EA281361, EA281362, EA281363, EA281364
    Molecule: Spiroxamine (8-tert-butyl-1,4-dioxaspiro[4.5]decan-3-yl)methyl-ethyl-propyl-amine
    Entries: EA278901, EA278902, EA278903, EA278904, EA278905, EA278906, EA278907, EA278908, EA278909, EA278910, EA278911, EA278912, EA278913, EA278914
    Molecule: Atraton 4-N-ethyl-6-methoxy-2-N-propan-2-yl-1,3,5-triazine-2,4-diamine
    Entries: EA015701, EA015702, EA015703, EA015704, EA015705, EA015706, EA015707, EA015708, EA015709, EA015710, EA015711, EA015712, EA015713, EA015714
    Molecule: Amisulpride 4-amino-5-esyl-N-[(1-ethylpyrrolidin-2-yl)methyl]-2-methoxy-benzamide
    Entries: EA285201, EA285202, EA285203, EA285204, EA285205, EA285206, EA285207, EA285208, EA285209, EA285210, EA285211, EA285212, EA285213, EA285214, EA285251, EA285252, EA285253, EA285254, EA285255, EA285256, EA285258, EA285259, EA285260, EA285264
    Molecule: Atrazine-2-hydroxy 2-Hydroxyatrazine 2-(ethylamino)-6-(isopropylamino)-1H-s-triazin-4-one
    Entries: EA027901, EA027902, EA027903, EA027904, EA027905, EA027906, EA027907, EA027908, EA027909, EA027910, EA027911, EA027912, EA027913, EA027914, EA027951, EA027952, EA027953, EA027954, EA027955, EA027956, EA027957, EA027958, EA027959, EA027960, EA027963
    Molecule: Mesotrione 2-[4-(methanesulfonyl)-2-nitrobenzoyl]cyclohexane-1,3-dione
    Entries: EA025801, EA025802, EA025803, EA025804, EA025805, EA025806, EA025807, EA025808, EA025809, EA025810, EA025811, EA025812, EA025813, EA025814, EA025851, EA025852, EA025853, EA025854, EA025855, EA025858, EA025859, EA025860, EA025864
    Molecule: Atrazine-desethyl-2-hydroxy 2-amino-6-(isopropylamino)-1H-s-triazin-4-one Prometon-Hydroxy-Desisopropyl
    Entries: EA027701, EA027702, EA027703, EA027704, EA027705, EA027706, EA027707, EA027708, EA027709, EA027710, EA027711, EA027712, EA027713, EA027714, EA027751, EA027752, EA027753, EA027754, EA027755, EA027756, EA027758, EA027759, EA027760, EA027761, EA027762
    Molecule: Dimethomorph (Z)-3-(4-chlorophenyl)-3-(3,4-dimethoxyphenyl)-1-(4-morpholinyl)-2-propen-1-one
    Entries: EA294401, EA294402, EA294403, EA294404, EA294405, EA294406, EA294407, EA294408, EA294409, EA294410, EA294411, EA294412, EA294413, EA294414
    Molecule: Pyrimethanil (4,6-dimethylpyrimidin-2-yl)-phenyl-amine
    Entries: EA271201, EA271202, EA271203, EA271204, EA271205, EA271206, EA271207, EA271208, EA271209, EA271210, EA271211, EA271212, EA271213, EA271214
    Molecule: Eprosartan 4-[[2-butyl-5-[(E)-3-hydroxy-3-keto-2-(2-thenyl)prop-1-enyl]imidazol-1-yl]methyl]benzoic acid
    Entries: EA277601, EA277602, EA277603, EA277604, EA277605, EA277606, EA277607, EA277608, EA277609, EA277610, EA277611, EA277612, EA277613, EA277614, EA277651, EA277652, EA277653, EA277654, EA277655, EA277656, EA277657, EA277658, EA277659, EA277660, EA277661, EA277662, EA277663, EA277664
    Molecule: Chloridazon 5-amino-4-chloro-2-phenyl-3-pyridazinone
    Entries: EA008801, EA008802, EA008803, EA008804, EA008805, EA008806, EA008807, EA008808, EA008809, EA008810, EA008811, EA008812, EA008813, EA008814, EA008851, EA008852, EA008853, EA008854, EA008855, EA008856, EA008857, EA008858, EA008859, EA008860, EA008861, EA008862, EA008863, EA008864
    Molecule: Terbutylazine-desethyl Desethylterbutylazine (4-amino-6-chloro-s-triazin-2-yl)-tert-butyl-amine
    Entries: EA067101, EA067102, EA067103, EA067104, EA067105, EA067106, EA067107, EA067108, EA067109, EA067110, EA067111, EA067112, EA067113, EA067114
    Molecule: Clomazone 2-(2-chlorobenzyl)-4,4-dimethyl-isoxazolidin-3-one
    Entries: EA013901, EA013902, EA013903, EA013904, EA013905, EA013906, EA013907, EA013908, EA013909, EA013910, EA013911, EA013912, EA013913, EA013914
    Molecule: Tramadol (1R,2R)-2-(dimethylaminomethyl)-1-(3-methoxyphenyl)-1-cyclohexanol
    Entries: EA256701, EA256702, EA256703, EA256704, EA256705, EA256706, EA256707, EA256708, EA256709, EA256710, EA256711, EA256712, EA256713, EA256714
    Molecule: Isoproturon-didemethyl [4-(propan-2-yl)phenyl]urea 1-(4-Isopropylphenyl)urea
    Entries: EA028501, EA028502, EA028503, EA028504, EA028505, EA028506, EA028507, EA028508, EA028509, EA028510, EA028511, EA028512, EA028513, EA028514
    Molecule: Propachlor 2-chloranyl-N-phenyl-N-propan-2-yl-ethanamide
    Entries: EA070801, EA070802, EA070803, EA070804, EA070805, EA070806, EA070807, EA070808, EA070809, EA070810, EA070811, EA070812, EA070813, EA070814
    Molecule: Cyclophosphamide N,N-bis(2-chloroethyl)-2-oxidanylidene-1,3,2$l^{5}-oxazaphosphinan-2-amine
    Entries: EA257901, EA257902, EA257903, EA257904, EA257905, EA257906, EA257907, EA257908, EA257909, EA257910, EA257911, EA257912, EA257913, EA257914
    Molecule: Terbutryn 2-N-tert-butyl-4-N-ethyl-6-methylsulfanyl-1,3,5-triazine-2,4-diamine
    Entries: EA030601, EA030602, EA030603, EA030604, EA030605, EA030606, EA030607, EA030608, EA030609, EA030610, EA030611, EA030612, EA030613, EA030614
    Molecule: Imazamox 2-(4-isopropyl-5-keto-4-methyl-2-imidazolin-2-yl)-5-(methoxymethyl)nicotinic acid
    Entries: EA293901, EA293902, EA293903, EA293904, EA293905, EA293906, EA293907, EA293908, EA293909, EA293910, EA293911, EA293912, EA293913, EA293914, EA293951, EA293952, EA293953, EA293954, EA293955, EA293956, EA293957, EA293958, EA293959, EA293960, EA293961, EA293962, EA293963
    Molecule: Prometryn 6-(methylthio)-N2,N4-di(propan-2-yl)-1,3,5-triazine-2,4-diamine
    Entries: EA013301, EA013302, EA013303, EA013304, EA013305, EA013306, EA013307, EA013308, EA013309, EA013310, EA013311, EA013312, EA013313, EA013314
    Molecule: 2-Aminobenzimidazole 1H-Benzimidazol-2-amine
    Entries: EA013801, EA013802, EA013803, EA013804, EA013805, EA013806, EA013807, EA013808, EA013809, EA013810, EA013811, EA013812, EA013813, EA013814, EA013851, EA013852, EA013853, EA013854, EA013855, EA013856, EA013857, EA013858, EA013859, EA013860, EA013861, EA013862, EA013863, EA013864
    Molecule: Thiacloprid (Z)-[3-[(6-Chloro-3-pyridinyl)methyl]-2-thiazolidinylidene]cyanamide Calypso
    Entries: EA295401, EA295402, EA295403, EA295404, EA295405, EA295406, EA295407, EA295408, EA295409, EA295410, EA295411, EA295412, EA295413, EA295414
    Molecule: Irgarol-descyclopropyl 2-N-tert-butyl-6-methylsulfanyl-1,3,5-triazine-2,4-diamine
    Entries: EA028301, EA028302, EA028303, EA028304, EA028305, EA028306, EA028307, EA028308, EA028309, EA028310, EA028311, EA028312, EA028313, EA028314
    Molecule: Sotalol N-[4-[1-hydroxy-2-(isopropylamino)ethyl]phenyl]methanesulfonamide
    Entries: EA017001, EA017002, EA017003, EA017004, EA017005, EA017006, EA017007, EA017008, EA017009, EA017010, EA017011, EA017012, EA017013, EA017014, EA017051, EA017052, EA017053, EA017054, EA017058, EA017059, EA017060, EA017064
    Molecule: 2,4-D 2-(2,4-dichlorophenoxy)acetic acid
    Entries: EA026751, EA026752, EA026753, EA026754, EA026755, EA026756, EA026757, EA026758, EA026759, EA026760, EA026761, EA026762, EA026763, EA026764
    Molecule: Imidacloprid N-[1-[(6-chloranylpyridin-3-yl)methyl]-4,5-dihydroimidazol-2-yl]nitramide
    Entries: EA270901, EA270902, EA270903, EA270904, EA270905, EA270906, EA270907, EA270908, EA270909, EA270910, EA270911, EA270912, EA270913, EA270914, EA270951, EA270952, EA270953, EA270954, EA270955, EA270956, EA270958, EA270959, EA270964
    Molecule: 1-Benzylpiperazine N-Benzylpiperazine 1-(phenylmethyl)piperazine
    Entries: EA282001, EA282002, EA282003, EA282004, EA282005, EA282006, EA282007, EA282008, EA282009, EA282010, EA282011, EA282012, EA282013, EA282014
    Molecule: Norfloxacin 1-ethyl-6-fluoro-4-oxo-7-piperazin-1-yl-1,4-dihydroquinoline-3-carboxylic acid
    Entries: EA023201, EA023202, EA023203, EA023204, EA023205, EA023206, EA023207, EA023208, EA023209, EA023210, EA023211, EA023212, EA023213, EA023214
    Molecule: Amitraz N'-(2,4-dimethylphenyl)-N-[(2,4-dimethylphenyl)iminomethyl]-N-methyl-formamidine
    Entries: EA010701, EA010702, EA010703, EA010704, EA010705, EA010706, EA010707, EA010708, EA010709, EA010710, EA010711, EA010712, EA010713, EA010714
    Molecule: 5-Methyl-1H-benzotriazole
    Entries: EA016701, EA016702, EA016703, EA016704, EA016705, EA016706, EA016707, EA016708, EA016709, EA016710, EA016711, EA016712, EA016713, EA016714, EA016751, EA016752, EA016753, EA016754, EA016755, EA016756, EA016757, EA016758, EA016759, EA016760, EA016761, EA016762, EA016763, EA016764
    Molecule: Trinexapac 4-[Cyclopropyl(hydroxy)methylene]-3,5-dioxocyclohexanecarboxylic acid
    Entries: EA015001, EA015002, EA015003, EA015004, EA015005, EA015006, EA015007, EA015008, EA015009, EA015010, EA015011, EA015012, EA015013, EA015014, EA015051, EA015052, EA015053, EA015054, EA015055, EA015057, EA015058, EA015059, EA015060, EA015061, EA015062
    Molecule: Atrazine-desethyl Deethylatrazine (4-amino-6-chloro-s-triazin-2-yl)-isopropyl-amine
    Entries: EA030901, EA030902, EA030903, EA030904, EA030905, EA030906, EA030907, EA030908, EA030909, EA030910, EA030911, EA030912, EA030913, EA030914
    Molecule: Naltrexone (5alpha)-17-(Cyclopropylmethyl)-3,14-dihydroxy-4,5-epoxymorphinan-6-one
    Entries: EA283001, EA283002, EA283003, EA283004, EA283005, EA283006, EA283007, EA283008, EA283009, EA283010, EA283011, EA283012, EA283013, EA283014
    Molecule: Trimethoprim 2,4-Diamino-5-(3,4,5-trimethoxybenzyl)pyrimidine 5-[(3,4,5-trimethoxyphenyl)methyl]pyrimidine-2,4-diamine
    Entries: EA019901, EA019902, EA019903, EA019904, EA019905, EA019906, EA019907, EA019908, EA019909, EA019910, EA019911, EA019912, EA019913, EA019914
    Molecule: Atenolol-desisopropyl 4-(3-Amino-2-hydroxypropoxy)phenylacetamide
    Entries: EA267001, EA267002, EA267003, EA267004, EA267005, EA267006, EA267007, EA267008, EA267009, EA267010, EA267011, EA267012, EA267013, EA267014
    Molecule: Acamprosate 3-acetamido-1-propanesulfonic acid
    Entries: EA284801, EA284802, EA284803, EA284804, EA284805, EA284806, EA284807, EA284808, EA284809, EA284810, EA284811, EA284812, EA284813, EA284814, EA284851, EA284852, EA284853, EA284854, EA284855, EA284856, EA284857, EA284858, EA284859, EA284860, EA284861, EA284862, EA284863, EA284864
    Molecule: Cilastatin
    Entries: EA255501, EA255502, EA255503, EA255504, EA255505, EA255506, EA255507, EA255508, EA255509, EA255510, EA255511, EA255512, EA255513, EA255514, EA255551, EA255552, EA255553, EA255554, EA255555, EA255556, EA255558, EA255559, EA255560, EA255561, EA255564
    Molecule: Neotame (3S)-3-(3,3-dimethylbutylamino)-4-[[(2S)-1-methoxy-1-oxidanylidene-3-phenyl-propan-2-yl]amino]-4-oxidanylidene-butanoic acid
    Entries: EA281501, EA281502, EA281503, EA281504, EA281505, EA281506, EA281507, EA281508, EA281509, EA281510, EA281511, EA281512, EA281513, EA281514, EA281551, EA281552, EA281553, EA281554, EA281555, EA281556, EA281557, EA281558, EA281559, EA281560, EA281561, EA281562, EA281563, EA281564
    Molecule: Ritalinic acid 2-phenyl-2-(2-piperidinyl)acetic acid
    Entries: EA070101, EA070102, EA070103, EA070104, EA070105, EA070106, EA070107, EA070108, EA070109, EA070110, EA070111, EA070112, EA070113, EA070114, EA070151, EA070152, EA070153, EA070158, EA070159
    Molecule: Iopamidol 1,3-Benzenedicarboxamide, N,N'-bis(2-hydroxy-1-(hydroxymethyl)ethyl)-5-((2-hydroxy-1-oxopropyl)amino)-2,4,6-triiodo- 1-N,3-N-bis(1,3-dihydroxypropan-2-yl)-5-(2-hydroxypropanoylamino)-2,4,6-triiodobenzene-1,3-dicarboxamide
    Entries: EA024101, EA024102, EA024103, EA024104, EA024105, EA024106, EA024107, EA024108, EA024109, EA024110, EA024111, EA024112, EA024113, EA024114
    Molecule: Fenpropimorph (2S,6R)-4-[3-(4-tert-butylphenyl)-2-methyl-propyl]-2,6-dimethyl-morpholine
    Entries: EA014601, EA014602, EA014603, EA014604, EA014605, EA014606, EA014607, EA014608, EA014609, EA014610, EA014611, EA014612, EA014613, EA014614
    Molecule: Perfluorooctyl phosphate 6:2PAP 1H,1H,2H,2H-perfluorooctyllphosphate
    Entries: EA292402, EA292403, EA292404, EA292405, EA292406, EA292407, EA292408, EA292409, EA292410, EA292411, EA292412, EA292413, EA292451, EA292452, EA292453, EA292454, EA292456, EA292458, EA292459, EA292460, EA292461, EA292462
    Molecule: Metazachlor 2-chloranyl-N-(2,6-dimethylphenyl)-N-(pyrazol-1-ylmethyl)ethanamide
    Entries: EA026901, EA026902, EA026903, EA026904, EA026905, EA026906, EA026907, EA026908, EA026909, EA026910, EA026911, EA026912, EA026913, EA026914
    Molecule: 1-(3-(Trifluoromethyl)phenyl)piperazine
    Entries: EA281901, EA281902, EA281903, EA281904, EA281905, EA281906, EA281907, EA281908, EA281909, EA281910, EA281911, EA281912, EA281913, EA281914
    Molecule: Citalopram 1-[3-(dimethylamino)propyl]-1-(4-fluorophenyl)-3H-2-benzofuran-5-carbonitrile
    Entries: EA290101, EA290102, EA290103, EA290104, EA290105, EA290106, EA290107, EA290108, EA290109, EA290110, EA290111, EA290112, EA290113, EA290114
    Molecule: Clindamycin (2S,4R)-N-[2-chloranyl-1-[(2R,3R,4S,5R,6R)-6-methylsulfanyl-3,4,5-tris(oxidanyl)oxan-2-yl]propyl]-1-methyl-4-propyl-pyrrolidine-2-carboxamide
    Entries: EA065701, EA065702, EA065703, EA065704, EA065705, EA065706, EA065707, EA065708, EA065709, EA065710, EA065711, EA065712, EA065713, EA065714
    Molecule: Atomoxetine (3R)-N-methyl-3-(2-methylphenoxy)-3-phenyl-1-propanamine
    Entries: EA284601, EA284602, EA284603, EA284604, EA284605, EA284606, EA284607, EA284608, EA284609, EA284610, EA284611, EA284612, EA284613, EA284614
    Molecule: Dexamethasone (8S,9R,10S,11S,13S,14S,16R,17R)-9-fluoranyl-10,13,16-trimethyl-11,17-bis(oxidanyl)-17-(2-oxidanylethanoyl)-6,7,8,11,12,14,15,16-octahydrocyclopenta[a]phenanthren-3-one
    Entries: EA261201, EA261202, EA261203, EA261204, EA261205, EA261206, EA261207, EA261208, EA261209, EA261210, EA261211, EA261212, EA261213, EA261214
    Molecule: Ethambutol (2S)-2-[2-[[(1S)-1-methylolpropyl]amino]ethylamino]butan-1-ol
    Entries: EA278201, EA278202, EA278203, EA278204, EA278205, EA278208, EA278209, EA278210, EA278211, EA278214
    Molecule: Levetiracetam (2R)-2-(2-ketopyrrolidino)butyramide
    Entries: EA256401, EA256402, EA256403, EA256404, EA256405, EA256406, EA256407, EA256408, EA256409, EA256410, EA256411, EA256412, EA256413, EA256414
    Molecule: Oseltamivir carboxylate (3R,4R,5S)-4-acetamido-5-amino-3-pentan-3-yloxycyclohexene-1-carboxylic acid
    Entries: EA065901, EA065902, EA065903, EA065904, EA065905, EA065906, EA065907, EA065908, EA065909, EA065910, EA065911, EA065912, EA065913, EA065914, EA065951, EA065952, EA065953, EA065954, EA065955, EA065956, EA065957, EA065958, EA065959, EA065960, EA065961, EA065962, EA065963
    Molecule: Tebuconazole 1-(4-chlorophenyl)-4,4-dimethyl-3-(1,2,4-triazol-1-ylmethyl)-3-pentanol
    Entries: EA032701, EA032702, EA032703, EA032704, EA032705, EA032706, EA032707, EA032708, EA032709, EA032710, EA032711, EA032712, EA032713, EA032714
    Molecule: Diflufenican N-(2,4-difluorophenyl)-2-[3-(trifluoromethyl)phenoxy]-3-pyridinecarboxamide
    Entries: EA011601, EA011602, EA011603, EA011604, EA011605, EA011606, EA011607, EA011608, EA011609, EA011610, EA011611, EA011612, EA011613, EA011614
    Molecule: Propachlor OXA N-(1-Methylethyl)-N-(phenyl)oxalamic acid 2-(N-isopropylanilino)-2-keto-acetic acid
    Entries: EA066601, EA066602, EA066603, EA066604, EA066605, EA066606, EA066607, EA066608, EA066609, EA066610, EA066611, EA066612, EA066613, EA066614
    Molecule: 2,6-Dichlorobenzamide 2,6-bis(chloranyl)benzamide
    Entries: EA008501, EA008502, EA008503, EA008504, EA008505, EA008506, EA008507, EA008508, EA008509, EA008510, EA008511, EA008512, EA008513, EA008514
    Molecule: Clothianidin 1-[(2-chloranyl-1,3-thiazol-5-yl)methyl]-2-methyl-3-nitro-guanidine 1-[(2-chloro-1,3-thiazol-5-yl)methyl]-2-methyl-3-nitroguanidine
    Entries: EA293301, EA293302, EA293303, EA293304, EA293305, EA293306, EA293307, EA293308, EA293309, EA293310, EA293311, EA293312, EA293313, EA293314, EA293351, EA293352, EA293353, EA293354, EA293355, EA293356, EA293357, EA293358, EA293359, EA293360, EA293361, EA293362, EA293363
    Molecule: Metoclopramide 4-amino-5-chloro-N-(2-diethylaminoethyl)-2-methoxy-benzamide
    Entries: EA278101, EA278102, EA278103, EA278104, EA278105, EA278106, EA278107, EA278108, EA278109, EA278110, EA278111, EA278112, EA278113, EA278114
    Molecule: Iodocarb Iodopropynyl butylcarbamate (IPBC) 3-iodanylprop-2-ynyl N-butylcarbamate
    Entries: EA030501, EA030502, EA030503, EA030504, EA030505, EA030506, EA030507, EA030508, EA030509, EA030510, EA030511, EA030512, EA030513, EA030514
    Molecule: Alachlor OXA 2-[(2,6-diethylphenyl)-(methoxymethyl)amino]-2-oxidanylidene-ethanoic acid
    Entries: EA010051, EA010052, EA010053, EA010054, EA010055, EA010056, EA010057, EA010058, EA010059, EA010060, EA010061, EA010062, EA010063, EA010064
    Molecule: 3,5,6-Trichloro-2-pyridinol 3,5,6-trichloro-1H-pyridin-2-one
    Entries: EA270451, EA270452, EA270453, EA270454, EA270455, EA270458, EA270459, EA270460, EA270461, EA270464
    Molecule: 1H-Benzotriazole
    Entries: EA016601, EA016602, EA016603, EA016604, EA016605, EA016606, EA016607, EA016608, EA016609, EA016610, EA016611, EA016612, EA016613, EA016614, EA016651, EA016652, EA016653, EA016654, EA016655, EA016656, EA016657, EA016658, EA016659, EA016660, EA016661, EA016662, EA016663, EA016664
    Molecule: 4-Acetamidoantipyrine N-(1,5-dimethyl-3-oxidanylidene-2-phenyl-pyrazol-4-yl)ethanamide
    Entries: EA023601, EA023602, EA023603, EA023604, EA023605, EA023606, EA023607, EA023608, EA023609, EA023610, EA023611, EA023612, EA023613, EA023614
    Molecule: N4-Acetylsulfadiazine N-[4-(2-pyrimidinylsulfamoyl)phenyl]acetamide
    Entries: EA024801, EA024802, EA024803, EA024804, EA024805, EA024806, EA024807, EA024808, EA024809, EA024810, EA024811, EA024812, EA024813, EA024814, EA024851, EA024852, EA024853, EA024854, EA024855, EA024856, EA024857, EA024858, EA024859, EA024860, EA024861, EA024862
    Molecule: Simazine-2-hydroxy 2-Hydroxysimazine 2,6-bis(ethylamino)-1H-1,3,5-triazin-4-one
    Entries: EA066801, EA066802, EA066803, EA066804, EA066805, EA066806, EA066807, EA066808, EA066809, EA066810, EA066811, EA066812, EA066813, EA066814, EA066851, EA066852, EA066853, EA066854, EA066855, EA066856, EA066857, EA066858, EA066859, EA066860, EA066861, EA066862, EA066863
    Molecule: Metoprolol 1-(isopropylamino)-3-[4-(2-methoxyethyl)phenoxy]propan-2-ol
    Entries: EA017201, EA017202, EA017203, EA017204, EA017205, EA017206, EA017207, EA017208, EA017209, EA017210, EA017211, EA017212, EA017213, EA017214
    Molecule: Atenolol acid Metoprolol acid 4-(2-Hydroxy-3-((1-methylethyl)amino)propoxy)benzeneacetic acid
    Entries: EA069701, EA069702, EA069703, EA069704, EA069705, EA069706, EA069707, EA069708, EA069709, EA069710, EA069711, EA069712, EA069713, EA069714, EA069751, EA069753
    Molecule: Metalaxyl 2-(N-(2-methoxy-1-oxoethyl)-2,6-dimethylanilino)propanoic acid methyl ester
    Entries: EA013501, EA013502, EA013503, EA013504, EA013505, EA013506, EA013507, EA013508, EA013509, EA013510, EA013511, EA013512, EA013513, EA013514
    Molecule: Lamotrigine 6-(2,3-dichlorophenyl)-1,2,4-triazine-3,5-diamine
    Entries: EA267601, EA267602, EA267603, EA267604, EA267605, EA267606, EA267607, EA267608, EA267609, EA267610, EA267611, EA267612, EA267613, EA267614
    Molecule: D617 Verapamil metabolite D617 2-(3,4-dimethoxyphenyl)-2-isopropyl-5-(methylamino)valeronitrile
    Entries: EA263701, EA263702, EA263703, EA263704, EA263705, EA263706, EA263707, EA263708, EA263709, EA263710, EA263711, EA263712, EA263713, EA263714
    Molecule: Propranolol 1-(1-naphthalenyloxy)-3-(propan-2-ylamino)-2-propanol
    Entries: EA017101, EA017102, EA017103, EA017104, EA017105, EA017106, EA017107, EA017108, EA017109, EA017110, EA017111, EA017112, EA017113, EA017114
    Molecule: N4-Acetylsulfamethazine N-[4-[(4,6-dimethyl-2-pyrimidinyl)sulfamoyl]phenyl]acetamide
    Entries: EA024701, EA024702, EA024703, EA024704, EA024705, EA024706, EA024707, EA024708, EA024709, EA024710, EA024711, EA024712, EA024713, EA024714, EA024751, EA024752, EA024753, EA024754, EA024755, EA024756, EA024757, EA024758, EA024759, EA024760, EA024761, EA024762, EA024763
    Molecule: Methomyl (1E)-N-(methylcarbamoyloxy)ethanimidothioic acid methyl ester
    Entries: EA294201, EA294202, EA294203, EA294204, EA294205, EA294206, EA294207, EA294208, EA294209, EA294210, EA294211, EA294212, EA294213, EA294214
    Molecule: Propazine 6-chloranyl-N2,N4-di(propan-2-yl)-1,3,5-triazine-2,4-diamine
    Entries: EA274101, EA274102, EA274103, EA274104, EA274105, EA274106, EA274107, EA274108, EA274109, EA274110, EA274111, EA274112, EA274113, EA274114
    Molecule: Prochloraz N-propyl-N-[2-(2,4,6-trichlorophenoxy)ethyl]-1-imidazolecarboxamide
    Entries: EA009601, EA009602, EA009603, EA009604, EA009605, EA009606, EA009607, EA009608, EA009609, EA009610, EA009611, EA009612, EA009613, EA009614
    Molecule: Fluroxypyr 2-(4-amino-3,5-dichloro-6-fluoropyridin-2-yl)oxyacetic acid
    Entries: EA013401, EA013402, EA013403, EA013404, EA013405, EA013406, EA013407, EA013408, EA013409, EA013410, EA013411, EA013412, EA013413, EA013414, EA013451, EA013452, EA013453, EA013454, EA013455, EA013456, EA013457, EA013458, EA013459, EA013460, EA013461, EA013462, EA013463, EA013464
    Molecule: N4-Acetylsulfathiazole N-[4-(1,3-thiazol-2-ylsulfamoyl)phenyl]acetamide
    Entries: EA024901, EA024902, EA024903, EA024904, EA024905, EA024906, EA024907, EA024908, EA024909, EA024910, EA024911, EA024912, EA024913, EA024914, EA024951, EA024952, EA024953, EA024959, EA024960
    Molecule: Dimethachlor 2-chloranyl-N-(2,6-dimethylphenyl)-N-(2-methoxyethyl)ethanamide
    Entries: EA070701, EA070702, EA070703, EA070704, EA070705, EA070706, EA070707, EA070708, EA070709, EA070710, EA070711, EA070712, EA070713, EA070714
    Molecule: N,N-Dimethylsulfamide [methyl(sulfamoyl)amino]methane MAS
    Entries: EA034101, EA034102, EA034103, EA034104, EA034105, EA034106, EA034107, EA034108, EA034109, EA034110, EA034111, EA034112, EA034113, EA034114
    Molecule: Propamocarb N-[3-(dimethylamino)propyl]carbamic acid propyl ester
    Entries: EA294501, EA294502, EA294503, EA294504, EA294505, EA294506, EA294507, EA294508, EA294509, EA294510, EA294511, EA294512, EA294513, EA294514
    Molecule: Acesulfame 2,2-diketo-6-methyl-oxathiazin-4-one
    Entries: EA275651, EA275652, EA275653, EA275654, EA275655, EA275656, EA275657, EA275658, EA275659, EA275660, EA275661, EA275662, EA275663, EA275664
    Molecule: N,N-Dimethyl-N'-p-tolylsulphamide 1-(dimethylsulfamoylamino)-4-methyl-benzene
    Entries: EA034001, EA034002, EA034003, EA034004, EA034005, EA034006, EA034007, EA034008, EA034009, EA034010, EA034011, EA034012, EA034013, EA034014, EA034051, EA034052, EA034053, EA034058, EA034059, EA034064
    Molecule: Rosuvastatin (E,3R,5R)-7-[4-(4-fluorophenyl)-2-(methyl-methylsulfonylamino)-6-propan-2-ylpyrimidin-5-yl]-3,5-dihydroxyhept-6-enoic acid
    Entries: EA280201, EA280202, EA280203, EA280204, EA280205, EA280206, EA280207, EA280208, EA280209, EA280210, EA280211, EA280212, EA280213, EA280214, EA280251, EA280252, EA280253, EA280254, EA280255, EA280256, EA280257, EA280258, EA280259, EA280260, EA280261, EA280262, EA280263, EA280264
    Molecule: Metolachlor 2-chloranyl-N-(2-ethyl-6-methyl-phenyl)-N-(1-methoxypropan-2-yl)ethanamide
    Entries: EA026801, EA026802, EA026803, EA026804, EA026805, EA026806, EA026807, EA026808, EA026809, EA026810, EA026811, EA026812, EA026813, EA026814
    Molecule: Fipronil 5-amino-1-[2,6-dichloro-4-(trifluoromethyl)phenyl]-4-(trifluoromethylsulfinyl)-3-pyrazolecarbonitrile
    Entries: EA266601, EA266602, EA266603, EA266604, EA266605, EA266606, EA266607, EA266608, EA266609, EA266610, EA266611, EA266612, EA266613, EA266614, EA266651, EA266652, EA266653, EA266654, EA266655, EA266656, EA266657, EA266658, EA266659, EA266660, EA266661, EA266662, EA266663, EA266664
    Molecule: Pravastatin (3R,5R)-7-[(1S,2S,6S,8S,8aR)-2-methyl-8-[(2S)-2-methylbutanoyl]oxy-6-oxidanyl-1,2,6,7,8,8a-hexahydronaphthalen-1-yl]-3,5-bis(oxidanyl)heptanoic acid
    Entries: EA285951, EA285952, EA285953, EA285954, EA285955, EA285956, EA285957, EA285958, EA285959, EA285960, EA285961, EA285962, EA285963, EA285964
    Molecule: Sebuthylazine 2-N-butan-2-yl-6-chloro-4-N-ethyl-1,3,5-triazine-2,4-diamine
    Entries: EA067001, EA067002, EA067003, EA067004, EA067005, EA067006, EA067007, EA067008, EA067009, EA067010, EA067011, EA067012, EA067013, EA067014
    Molecule: Dimethachlor OXA 2-[(2,6-dimethylphenyl)(2-methoxyethyl)amino]-2-oxo-acetic acid
    Entries: EA253601, EA253602, EA253603, EA253604, EA253605, EA253606, EA253607, EA253608, EA253609, EA253610, EA253611, EA253612, EA253613, EA253614, EA253651, EA253652, EA253653, EA253654, EA253655, EA253656, EA253657, EA253658, EA253659, EA253660, EA253661, EA253662, EA253663, EA253664
    Molecule: Linuron 3-(3,4-dichlorophenyl)-1-methoxy-1-methyl-urea
    Entries: EA016001, EA016002, EA016003, EA016004, EA016005, EA016006, EA016007, EA016008, EA016009, EA016010, EA016011, EA016012, EA016013, EA016014

