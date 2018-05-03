# a test set of MGF files for ChemDistiller
Files with extension name  '.mgf' can be converted to chemdistiller input files like txt files under folder 
<chemdistiller folder>/data/MassBankTestSpectra/ .

Open the terminal and type in the following command:
`python MGF2ChemDistiller.py <input_folder_fullname or input_file_fullname> (<output_folder_fullname>)`

MGF2ChemDistiller.py module will create a single file for each spectra.

The default output folder is a created dirctory named 'default' under the same directory as the input folder or file, and 
new created files have the same name as its parent mgf files.

Professor Fiehn’s laboratory owns datasets under this directory.