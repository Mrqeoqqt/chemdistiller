# README #

##ChemDistiller v0.1.

__ChemDistiller__ is the high-throughput annotation engine for
tandem MS spectra annotation. 

##Installing:

__ChemDistiller__ supports Python 2 and 3 (64-bit version recommended) and requires _SciPy_, _NumPy_, _h5py_ libraries.

We strongly recommend installing one of the versions of 
[Anaconda](https://www.anaconda.com/download/) python from: https://www.anaconda.com/download/
as it comes with many dependencies already pre-installed.

[h5py](http://www.h5py.org/) can be installed in [Anaconda](https://www.anaconda.com/download/) environment by calling: 

```
conda install -c anaconda h5py 
```

__ChemDistiller__ benefits from installed and available [RDKit](http://www.rdkit.org/),
however it is optional ([RDKit](http://www.rdkit.org/) is currently used to generate
images of 2D chemical structures of compounds for the report).

You can learn how to install [RDKit](http://www.rdkit.org/) for [Anaconda](https://www.anaconda.com/download/) here:
http://www.rdkit.org/docs/Install.html

##Getting Started:
 1. To use __ChemDistiller__ - check it out or download into one of your folders 
    on your PC or Mac (e.g. /home/yourname/ChemDistiller/.). "cd" to your __ChemDistiller__ directory.
    Check that your python is installed and available by typing: 

    ```
    python --help
    ```

 2. Try running __ChemDistiller__ from the command line:
    
    ```
    python annotate.py --help
    ```

    This should give you the list of currently available command line arguments.

 3. The simplest way to check the engine is working is to run it in test mode:

    ```
    python annotate.py --test
    ```

    This will initiate a test annotation run using default settings and the test
    data provided in .../ChemDistiller/data/MassBankTestSpectra/*.
    The test should take a few minutes for 863 spectra provided and the results
    will be placed in ../ChemDistiller/testlogs/ in subfolder named with current date and time.
    Since these compounds also have the known chemical structures (known compounds)
    the retrieval statistics should also be printed.
    The output folder should contain:

       * annotated_spectra.json - JSON formatted annotations;

       * log.txt - log of the current run, duplicating console output;
 
       * annotated_specta - folder containing individual annotated spectra in internal hierarchial text format;
 
       * Report - folder containing HTML annotation report. In it there should be _index.html_ file which you can open in any web browser such as _FireFox_,
         _Chrome_, _Opera_ or _Internet Explorer_ and navigate the hyperlinks to view 
         different aspects of the report. 

 4. To analyse your spectra you would need to run ChemDistiller as follows:

    ```
    python annotate.py Input_Spectra_folder [output_folder] --ncpu 5 --delta_mz 20 --max_results 100
    ```

    where:

      * _Input\_Spectra\_folder_ is the full path to the folder containing your spectra to be annotated;
      * Optional _output\_folder_ specifies where to write your results to (if not supplied will be directed to a subfolder in .../ChemDistiller/testlogs/ named with current date and time); 
      * _--ncpu_ tells how many processors to use in the multiprocessor/multicore environment. Default: 1, maximum for your PC will be shown in help screen (see `python annotate.py --help`);
      * _--delta\_mz_ - tolerance for the MS1 peak _m/z_ values when searching in the databases, _i.e._ candidates should be withing +/- _delta\_mz_ from MS1 peak _m/z_ value (_delta\_mz_ is in ppm);
      * _--max\_results_ - limit the number of candidate compounds per MS1 peak to this value. Default: 10, but you may wish to increase it for large compound databases.


At the moment your spectra will need to be in __ChemDistiller__ internal format, which is very simple
hierarchical format (see [int_format.md](int_format.md)). We are working hard to implement conversions from other formats and they should be available in the next few weeks workload permitting. 

##Additonal Chemical Compound Databases and Support Vector Machines (SVMs) for FingerScorer in ChemDistiller

Current version of __ChemDistiller__ is supplied with only a small set of __chemical compound databases__ and
simple __SVM models__ due to size limitations. In order to perform larger scale analysis and use
more complex (and thus more slow and large) __SVMs__ - please download them from:

* __Compound Databases__ in HDF5 format (*.h5): https://www.mediafire.com/folder/v5l4380gqbvie/DBs .
  Databases when downloaded and unpacked where needed (use [7-Zip](http://www.7-zip.org/) or similar) should be copied to .../ChemDistiller/DBs/. subfolder. They will be picked up automatically when you run __ChemDistiller__ the
  next time.

* __SVMs__: https://www.mediafire.com/folder/v4lb8s2nns9c6/SVMs .
  You can choose _Linear\_weighted_ (already installed) or _Radial\_weighted_ __SVMs__ here. _Linear_ are smaller, faster, but prediction power is worse than for _Radial_. _Radial_ are better (and thus strongly recommended), but larger and slower. If replacing the default __SVMs__ - make sure to delete both __"1"__ and __"-1"__ subfolders in .../ChemDistiller/SVMs/. to avoid confusion prior to copying the new ones. When unpacked, __SVMs__ (subfolders __"1"__ and __"-1"__) should be copied into .../ChemDistiller/SVMs/. to be used by default. Alternatively you can unpack it somewhere on your disk
  and point __ChemDistiller__ to it by adding a command line argument _--svm\_folder_ _"Your\_folder\_with\_SVMs"_ . 
  
  Example:
  ```
  python annotate.py /home/me/my_spectra /home/me/my_results --svm_folder /home/me/my_svms
  ```

##Legality and Disclaimers
Please note that our databases are generated by our scripts using original chemical compound databases openly available online. 
Some of the original databases may have individual restrictions regarding commercial use and
thus it is not always clear if the derivative database could be used for commercial purposes either.
So, please consult _readme.txt_ in https://www.mediafire.com/folder/v5l4380gqbvie/DBs
and legal disclaimers of the corresponding source compound databases prior to attemting
to use them for commercial purposes. Also please do not forget to properly cite the original 
databases as well as ChemDistiller in your publications:

##Acknowledgments

Any papers describing data analytics using this package, or whose results were significantly aided by the use of this package (except when the use was internal to a larger program), should include an acknowledgment and citation to the following manuscript(s):

1. Ivan Laponogov, Noureddin Sadawi, Dieter Galea, Reza Mirnezami & Kirill A. Veselkov. __ChemDistiller: an engine for metabolite annotation in mass spectrometry__ (2017) _submitted_. 

If the package is used prior to its open source release and publication, please contact [kirill.veselkov04@imperial.ac.uk](mailto:kirill.veselkov04@imperial.ac.uk) to ensure that all contributors are fairly acknowledged.


Furthermore, provided engine, generated __databases__ and __SVMs__ are provided __AS IS__ and 
authors accept no responsibility for whatever they are going to be used for or any consequences.
Use them at your own risk. 


##Authors

Lead Developer: Dr. Ivan Laponogov [i.laponogov@imperial.ac.uk](mailto:i.laponogov@imperial.ac.uk) 
Chief project investigator: Dr. Kirill Veselkov [kirill.veselkov04@imperial.ac.uk](mailto:kirill.veselkov04@imperial.ac.uk)

##Issues

Please report any bugs or requests that you have using the Bitbucket issue tracker!

##Development

Please include TODO list here! 

####License:
The code is released under permissive BSD type license.

####Funding:
This project was funded by [__European Comission__](https://ec.europa.eu/commission/index_en) as part of [__Horizon2020__](https://ec.europa.eu/programmes/horizon2020/) programme
and is a part of [__METASPACE__](http://metaspace2020.eu/) project.
