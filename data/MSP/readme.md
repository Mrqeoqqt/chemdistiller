# a test set of MSP files for ChemDistiller
Files with extension name  '.txt' and MSP format can be converted to chemdistiller input files like txt files under folder 
<chemdistiller folder>/data/MassBankTestSpectra/ .

Open the terminal and type in the following command:
`python MSP2ChemDistiller.py <input_folder_fullname or input_file_fullname> (<output_folder_fullname>)`

MSP2ChemDistiller.py module will create a single file for each spectra.

The default output folder is a created dirctory named 'default' under the same directory as the input folder or file, and 
new created files are named after its molecular name or its filename if molecular name is null.

Professor Fiehn’s laboratory owns datasets under this directory.