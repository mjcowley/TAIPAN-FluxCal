![alt text](http://web.science.mq.edu.au/~mcowley/taipan/taipan.jpg "TAIPAN LOGO")

Python program to compute sensitivity functions for [TAIPAN](http://www.taipan-survey.org/) . 
---
### Directories and files

**spectra/blue/** - spectra of F stars in the blue arm extracted from the files that 2dfdr produces (.fits)  
**spectra/red/** - spectra of F stars in the red arm extracted from the files that 2dfdr produces (.fits)  
**SEDs/** - seven synthetic SEDs for F stars   
**filters/** - filters curves (currently for DES and OmegaCam)  
**photometry/** - photometry of F stars to determine g-r colours (.csv)  
**extinction/aao_extinction.tab** - file with atomospheric extinction at the AAT 

### Quick start

To create the sensitivity curves, run:

```python
python fstars.py -c fstars.config -s g-r.txt -p -l 18.0 -o red.fits -f Savitzky_Golay -a red -z -31.2
python fstars.py -c fstars.config -s g-r.txt -p -l 18.0 -o blue.fits -f Savitzky_Golay -a blue -z -30.5
```

You'll get three plots when running the code: 

**figure_1.png** - The individual sensitivity functions and the variance  
![alt text](http://web.science.mq.edu.au/~mcowley/taipan/figure1.png "Figure 1")
**figure_2.png** - The scatter in the sensitivity function at three wavelengths  
![alt text](http://web.science.mq.edu.au/~mcowley/taipan/figure2.png "Figure 2")
**figure_3.png** - The mean sensitivity function  
![alt text](http://web.science.mq.edu.au/~mcowley/taipan/figure3.png "Figure 3")

**-c** is the directory configuration file (fstars.config)  
**-s** is the file that associates a star of a given colour to a spectral type (g-r.txt)  
**-p** is to produce the plots  
**-l** is a magnitude limit (i.e. stars need to be brighter than this to be accepted)  
**-o** is the name of the output file  
**-f** is the smoothing filter used on the sensitivity curves (Savitzky_Golay, Gaussian or Poly)  
**-a** is the arm to process (red or blue)  
**-z** the zero point limit of the observations (i.e. observations taken with cloud or in bad seeing will have low ZPs  
and are excluded NOTE: The ZP is not computed by 2dfdr, but is included in this DES example)  

To see the above, run:

```python
python fstars.py --help
```

For each star, the code performs the following:

* Select star if bright enough and the conditions were good enough  
* From the F star catalogue, determine the g-r colour  
* Using the g-r colour, find the best spectral template  
* Warp the spectral template to match the colours of the star exactly  
* Apply an extinction correction  
* Compute the sensitivity curve for that star and smooth it  
* Compute the mean sensitity curve and its variance  
  
This code is a modified version of the OsDES flux calibration code, written by Chris Lidman of the [AAO](https://www.aao.gov.au/) . 

