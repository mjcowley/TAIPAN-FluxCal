![alt text](http://web.science.mq.edu.au/~mcowley/taipan/taipan.jpg "TAIPAN LOGO")

Python code for flux calibration in the TAIPAN Spectroscopic Survey

Original by Michael Cowley

Copyright 2016, [TAIPAN Survey Team](http://www.taipan-survey.org/)

---
### Directories and files

**SEDs/** - seven synthetic SEDs for F stars   
**filters/** - instrument filter curves  
**photometry/** - photometry of F stars to determine g-r colours (.csv)  
**extinction/aao_extinction.tab** - file with atomospheric extinction at the AAT 

### Quick start

To generate spectral response curves from flux calibration, run:

```python
python fstars.py
```

### Command-line options

The following command-line options are available:

**-c** is the configuration file for directories and files seen above (fstars.config)  
**-s** is the file that associates a star of a given colour to a spectral type (g-r.txt)  
**-p** is to produce diagnostic plots  
**-l** is a magnitude limit (i.e. stars need to be brighter than this to be accepted)  
**-f** is the smoothing filter used on the sensitivity curves (Savitzky_Golay or Spline)  
**-d** is the directory containing the two spectra directories (ccd_1 and ccd_2)  
**-o** is the name of the output file (written in -d directory)  
**-a** is to arm to process (ccd_1 or ccd_2)  
**-z** is the zero point limit of the observations  

To generate spectral response curves from flux calibration, run:

```python
python fstars.py -c fstars.config -s g-r.txt -p -l 18.0 -f Spline -d spectra/ -o src_ccd_1.fits -a ccd_1 -z -30.5
```

If selected, the following diagnostics plots will be generated: 

**figure_1** - The individual sensitivity functions and the variance  
![alt text](http://web.science.mq.edu.au/~mcowley/taipan/figure1.png "Figure 1")

**figure_2** - The scatter in the sensitivity function at three wavelengths  
![alt text](http://web.science.mq.edu.au/~mcowley/taipan/figure2.png "Figure 2")

**figure_3** - The mean sensitivity function  
![alt text](http://web.science.mq.edu.au/~mcowley/taipan/figure3.png "Figure 3")

### Help Message

The help message can be seen using ```python fstars.py --help```:

```
         ___  ____  ____  ____  ____                
_________\__\/  __\/  __\/  __\/  __\_______________
__TAIPAN_______/  /__/  /__/  /__/  /__FluxCal______
           \__/  \__/  \__/  \__/  \   <> \         
                                    \_____/--<      

Usage: fstars.py [options]

Options:
  -c CONFIG, --config=CONFIG, Configuration file
  -d DIRECTORY, --directory=DIRECTORY, Spectra Directory
  -s SENS, --sens=SENS, Sensitivity function
  -o OUTPUT, --output=OUTPUT, Sensitivity function
  -l LIMIT, --limit=LIMIT, Magnitude limit
  -f FILTER, --filter=FILTER, Filter
  -a ARM, --arm=ARM, Arm to process
  -z ZPLIMIT, --zplimit=ZPLIMIT, Zero point limit
  -n, --noextinction, Extinction Correction
  -p, --plot, Plot diagnostics

```

For each star, the code performs the following:

* Select star if bright enough and the conditions were good enough  
* From the F star catalogue, determine the g-r colour  
* Using the g-r colour, find the best spectral template  
* Warp the spectral template to match the colours of the star exactly  
* Apply an extinction correction  
* Compute the sensitivity curve for that star and smooth it  
* Compute the mean sensitivity curve and its variance  
  
### Help Message
 
 Diagrammatic representation of TAIPAN-FluxCal
![alt text](http://web.science.mq.edu.au/~mcowley/taipan/taipan_flowchart.png "Flow Chart")

This code is a modified version of the OsDES flux calibration code, written by Chris Lidman of the [AAO](https://www.aao.gov.au/) . 

