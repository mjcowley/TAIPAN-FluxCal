# coding=utf-8
from optparse import OptionParser
import library as LIB
import os
import pyfits
import numpy
import warnings
import scipy as sp
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib as mpl
import matplotlib.pyplot as plt 
from scipy.stats import norm
from scipy import ndimage
from math import factorial
mpl.rcParams['font.family'] = 'RomanD'
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams.update({'font.size': 10})

# ======================================================================================================================
# ======  ╔╦╗╔═╗╦╔╗╔  ╔═╗╔═╗╔╦╗╔═╗  ====================================================================================
# ======  ║║║╠═╣║║║║  ║  ║ ║ ║║║╣   ====================================================================================
# ======  ╩ ╩╩ ╩╩╝╚╝  ╚═╝╚═╝═╩╝╚═╝  ====================================================================================
# ======================================================================================================================

'''
Python program to compute sensitivity functions for Taipan
By Michael Cowley, michael.cowley@students.mq.edu.au
'''

# ======================================================================================================================
# ===== Supress Warnings ===============================================================================================
# ======================================================================================================================

warnings.filterwarnings('ignore', message='Overwriting existing file')
numpy.seterr(divide='ignore', invalid='ignore')

# ======================================================================================================================
# ===== Globals ========================================================================================================
# ======================================================================================================================

filterNames=['u','g','r','i','z']
filterTransCurves={'u':'DES_u.dat','g':'DES_g.dat','r':'DES_r.dat','i':'DES_i.dat','z':'DES_z.dat'}

# ======================================================================================================================
# ===== Input Parameters ===============================================================================================
# ======================================================================================================================

parser = OptionParser()
parser.add_option("-c", "--config", dest="config", default=None,
                  help="Configuration file")
parser.add_option("-s", "--sens", dest="sens", default=None,
                  help="Sensitivity function")
parser.add_option("-o", "--output", dest="output", default='test.fits',
                  help="Sensitivity function")
parser.add_option("-l", "--limit", dest="limit", default=None,
                  help="Magnitude Limit")
parser.add_option("-f", "--filter", dest="filter", default='Savitzky_Golay',
                  help="filter")
parser.add_option("-a", "--arm", dest="arm", default='spliced',
                  help="Arm to process")
parser.add_option("-z", "--zplimit", dest="zplimit", default=-100.0,
                  help="Zero point limit")
parser.add_option("-n", "--noextinction", dest="noextinction", action="store_true", default=False,
                  help="Extinction correction")
parser.add_option("-p", "--plot", dest="plot", action="store_true", default=False,
                  help="plot intermediate results")
(options, args) = parser.parse_args()

param=LIB.buildDictionary(options.config)

# ======================================================================================================================
# ===== Functions ======================================================================================================
# ======================================================================================================================

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Savitzky Golay Smoothing Algorithm"""
    try:
        window_size = numpy.abs(numpy.int(window_size))
        order = numpy.abs(numpy.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    b = numpy.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = numpy.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    firstvals = y[0] - numpy.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + numpy.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = numpy.concatenate((firstvals, y, lastvals))

    return numpy.convolve( m[::-1], y, mode='valid')


def correctExtinction(hdr,extinction):
    """Correct for Atmospheric Extinction"""
    wave=hdr['CRVAL1']+(numpy.arange(hdr['NAXIS1']) - hdr['CRPIX1']+1) * hdr['CDELT1']
    airmass=1.0/numpy.cos((hdr['ZDSTART']+hdr['ZDEND'])/2.*numpy.pi/180.)
    correction=interp1d(extinction.wave,10**(extinction.extinction * airmass/ 2.5))(wave)

    return correction


def sens(param,filterTransCurves,extinction,templates,catalogue,options):
    """Generate Sensitivity Curves"""
    fits=numpy.array([],float)
    nfit=0
    print options.arm
    directory=param['outputdir']+'/'+options.arm+'/'
    start={'red':5650,'blue':3740,'spliced':3740}
    end={'red':9000,'blue':5860,'spliced':9000}

    notInUse=0
    if options.arm=='blue':
        scale=0.9
    else:
        scale=1.0

    for star in os.listdir(directory):
        if 'fits' in star:
            data,hdr=pyfits.getdata(directory+star,0,header=True)
            var=pyfits.getdata(directory+star,1,header=False)
            mag=getmag(hdr['OBJECT'],catalogue)
            use=False
            ZPkey={'red':'REDZP','blue':'BLUZP'}
            cts = numpy.nanmedian(data[1300:1350])
            if len(mag) > 3:
                if mag[2] > 10. and mag[2] < float(options.limit) and mag[1]>0 and hdr[ZPkey[options.arm]] < float(options.zplimit):
                    use=True
            if not use:
                notInUse+=1
                print 'Not using %s with ZP %6.2f' % (star,hdr[ZPkey[options.arm]])
            else:
                # Use one of the colours to select the approrpiate F star template
                print 'Using %s with a zeropoint of %6.2f' % (star,hdr[ZPkey[options.arm]])
                template=param['seds']+'/'+findTemplate(templates,mag[1]-mag[2])+'.fits'
                fstar=LIB.spectrum()
                fstar.read(template)
                # Scale the template to the g magnitude of the star
                # Template spectra are in F_lambda
                fstar.flux=fstar.flux*3.631e-9/10**(mag[1]/2.5)
                # Warp the synthetic spectrum so that it matches the observed colours
                # and rebin to the wavelength scale of the data.
                synthetic=warp(fstar,filterTransCurves,mag,hdr)

                # Compute the extinction correction
                if not options.noextinction:
                    correction=correctExtinction(hdr,extinction)
                    data=data*correction
                    var=var*correction**2.
                # Filter the sensitivity funciotn
                if options.filter =='Savitzky_Golay':
                    # Smooth using the Savitzky Golay algorithm
                    sg=computeSens(synthetic,data,var)
                    wave=numpy.arange(start[options.arm],end[options.arm],scale)
                    fit=interp1d(synthetic.wave,sg,bounds_error=False,fill_value=numpy.nan)(wave)
                elif options.filter =='Gaussian': # 
                    # Smooth using a Gaussian filter after removing points with excessive variance
                    sg=computeSensGauss(synthetic,data,var)
                    wave=numpy.arange(start[options.arm],end[options.arm],scale)
                    fit=interp1d(synthetic.wave,sg,bounds_error=False,fill_value=numpy.nan)(wave)
                elif options.filter == 'Poly':  #
                    # Smooth using a least squares polynomial fit.
                    sg = computeSensPoly(synthetic, data, var)
                    wave = numpy.arange(start[options.arm], end[options.arm], scale)
                    fit = interp1d(synthetic.wave, sg, bounds_error=False, fill_value=numpy.nan)(wave)

                # Scale the fits by the median value
                fits=numpy.append(fits,fit/numpy.nanmedian(fit))
                nfit+=1


    fit=fits.reshape(nfit,len(fits)/nfit)
    print 'Using %d stars' % nfit
    print 'Exlcuded %d stars' % notInUse
    print
    # Compute a robust RMS, exclude Nans
    rms=numpy.zeros(len(fits)/nfit,float)
    medianfit=numpy.zeros(len(fits)/nfit,float)
    for i in range(len(fits)/nfit):
        sample=fit[:,i]
        good=~numpy.isnan(sample)
        y=numpy.sort(sample[good])
        if len(y) > 0:
            index_low=int(0.185 * len(y))
            index_high=len(y)-index_low-1
            rms[i]=(y[index_high]-y[index_low]) / 2.0
            medianfit[i]=y[int(0.5*len(y))]
        else:
            rms[i]=numpy.nan
            medianfit[i]=numpy.nan
            
    #print numpy.isnan(medianfit).any()

    if options.plot:

        fig=plt.figure() 
        ax=fig.add_subplot(212)
        #plot the average fit and the percentage deviation from the average fit (two plots)
        ax.plot(wave[10:-10],(rms/medianfit)[10:-10])
        print 'Median RMS %4.3f' % numpy.nanmedian((rms/medianfit)[10:-10])
        ax.set_xlabel('Wavelength')
        ax.set_ylim(0,0.5)
        ax.set_ylabel('scatter')
        ax.set_title('Fractional variation in the sensitivity function')
        ax=fig.add_subplot(211)
        ax.set_ylim(0,3.0)
        for i in range(nfit):
            ax.set_ylabel('Sensitivity')
            ax.plot(wave[10:-10],fit[i,10:-10])


        # Compare distribution at some location with that of a Gaussian
        fig=plt.figure()
        plots=[{'figure':311,'loc':500},{'figure':312,'loc':1000},{'figure':313,'loc':1500}]
        for plot in plots:
            ax=fig.add_subplot(plot['figure'])
            step=0.01
            loc=plot['loc']
            bin=numpy.arange(medianfit[loc]-0.5,medianfit[loc]+0.5,step)
            ax.hist(fit[:,loc],bin)
            rv=norm(medianfit[loc],rms[loc])
            ax.plot(bin,len(fit[:,loc])*rv.pdf(bin)*step,label='%d' % loc)
            ax.legend()
            ax.set_xlabel('sensitivity')
            ax.set_ylabel('frequency')
        
    

        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(wave,medianfit,label=options.filter+' filtered')
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Sensitivity')
        ax.set_ylim([0, 1.8])
        ax.legend(loc='upper left')


        plt.show()
        plt.close()


    # Write out the sensitivity curve
    hdr['CDELT1']=scale
    hdr['CRVAL1']=wave[0]
    hdr['CRPIX1']=1.0
    
    pyfits.writeto(options.output,medianfit,hdr,clobber=True)

    return


def getmag(fstar,catalogue):
    """Determine Star's Magnitude"""
    mag=[]
    # If there is no measurement, the results is set to -99
    for star in catalogue:
        if star['Candidate Name']==fstar and star['removed not Fstar']=='':
            for filter in ['DES <u> SingleEpoch','CoAdd mag_auto g','CoAdd mag_auto r','CoAdd mag_auto i','CoAdd mag_auto z']:
                if filter in star.keys():
                    if float(star[filter]) > 10.0 and float(star[filter]) < 30.0:
                        mag.append(float(star[filter]))
                    else:
                        mag.append(-99.9)
                else:
                    mag.append(-99.9)
            break # Some entries are in the catalogue twice, we take the first entry
             

    return numpy.array(mag)


def computeSens(synthetic,data,var):
    """Compute the Sensitivity Curve with Savitsky Golay"""
    # Divide by the synthetic curve
    # First use a median filter
    # good=~(numpy.isnan(data) | numpy.isnan(var))
    medfilter=21
    y=medfilt(data/synthetic.flux/1.e19,medfilter)
    return savitzky_golay(y, 31, 1)


def computeSensGauss(synthetic,data,var):
    """Compute the Sensitivity Curve with Gaussian"""
    # interpolate over pixels with excess variance 
    vmed = numpy.nanmedian(var) 
    mask = var > 3*vmed     # true for values that need to be interpolated over
    xp = numpy.arange(len(data), dtype=int)[~mask]
    fp = data[~mask]
    x = numpy.arange(len(data), dtype=int)[mask]	# indices to interpolate over
    data[mask] = numpy.interp(x, xp, fp)
    return ndimage.gaussian_filter(data/synthetic.flux/1.e19, sigma=(51), order=0)


def computeSensPoly(synthetic,data,var):
    """Compute the Sensitivity Curve with Polynomial"""
    # High-order polynomials may oscillate wildly
    x = numpy.arange(len(data), dtype=int)
    y = data / synthetic.flux / 1.e19
    p = sp.polyfit(x, y, deg=15)
    return sp.polyval(p, x)


def medfilt (x, k):
    """Apply a length-k median filter to a 1D array x. Boundaries are extended by repeating endpoints"""
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = numpy.zeros ((len (x), k), dtype=x.dtype)
    y[:,k2] = x
    for i in range (k2):
        j = k2 - i
        y[j:,i] = x[:-j]
        y[:j,i] = x[0]
        y[:-j,-(i+1)] = x[j:]
        y[-j:,-(i+1)] = x[-1]
    return numpy.nanmedian (y, axis=1)

def findTemplate(templates,colour):
    """Find Best Template"""
    selectedTemplate = None
    min=1.0
    for template in templates:
        diff=abs(float(template['Colour']) - colour) 
        if diff < min :
            min=diff
            selectedTemplate=template['SpectralType'] 

    return selectedTemplate

def warp(fstar,filters,mag,hdr):
    """Warp the spectral template to match the broad band photometry"""
    # First, compute AB magnitudes and effective wavelengths
    magAB=numpy.array([],float)
    effWavelength=numpy.array([],float)
    for filter in filters.keys():
        magAB=numpy.append(magAB,filters[filter].computeABmag(fstar))
        effWavelength=numpy.append(effWavelength,filters[filter].computeEffectiveWavelength(fstar))

    # Take care with the way filters.keys() orders the filters
    # Sort in terms of increasing wavelength, otherwise the spline interpolation does not work
    indices=numpy.argsort(effWavelength)
    magDiff=numpy.array([],float)
    for i in range(len(mag)):
        if mag[i] < 0:
            magDiff=numpy.append(magDiff,0.0)
        else:
            magDiff=numpy.append(magDiff,mag[i]-magAB[indices[i]])
    
    s = InterpolatedUnivariateSpline(effWavelength[indices], (magDiff))
    synthetic=LIB.spectrum()
    synthetic.wave=hdr['CRVAL1']+(numpy.arange(hdr['NAXIS1']) - hdr['CRPIX1']+1) * hdr['CDELT1']
    ynew=s(synthetic.wave)
    synthetic.flux=interp1d(fstar.wave,fstar.flux)(synthetic.wave) / 10.**(ynew / 2.5)

    return synthetic

# ======================================================================================================================
# ===== Main Code ======================================================================================================
# ======================================================================================================================

if options.sens != None:
    # Read in the Fstar catalogue
    catalogue,keys=LIB.readData(param['catalogue'],'csv')
    # Read in the filter curves
    filterCurves=LIB.readFilterCurves(filterNames,param['filters']+'/',filterTransCurves)
    # Read in colours and spectral types
    templates,keys=LIB.readData(options.sens,'asciicolumn')
    # Read in the extinction curve
    extinction=LIB.extinction()
    extinction.read(param['extinction'])
    # Compute sensitivity functions
    sens(param,filterCurves,extinction,templates,catalogue,options)
    
