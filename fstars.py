# coding=utf-8
#import library as lib
import os
import pyfits
import numpy
import warnings
import scipy as sp
import matplotlib.pyplot as plt
from optparse import OptionParser
from termcolor import colored
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.stats import norm
from scipy import ndimage
from math import factorial


print ""
print colored('########    ###    #### ########     ###    ##    ##', 'red')
print colored('   ##      ## ##    ##  ##     ##   ## ##   ###   ##', 'yellow')
print colored('   ##     ##   ##   ##  ##     ##  ##   ##  ####  ##', 'red')
print colored('   ##    ##     ##  ##  ########  ##     ## ## ## ##', 'yellow')
print colored('   ##    #########  ##  ##        ######### ##  ####', 'red')
print colored('   ##    ##     ##  ##  ##        ##     ## ##   ###', 'yellow')
print colored('   ##    ##     ## #### ##        ##     ## ##    ##', 'red')
print colored('         ___  ____  ____  ____  ____                ', 'green')
print colored('_________\__\/  __\/  __\/  __\/  __\_______________', 'green')
print colored('____________/  /__/  /__/  /__/  /__________________', 'green')
print colored('           \__/  \__/  \__/  \__/  \   <> \         ', 'green')
print colored('                                    \_____/--<      ', 'green')
print ""

'''
Python program to compute sensitivity functions for Taipan
By Michael Cowley, michael.cowley@students.mq.edu.au
'''

# ======================================================================================================================
# ===== Suppress Overwrite Warning =====================================================================================
# ======================================================================================================================

warnings.filterwarnings('ignore', message='Overwriting existing file')
numpy.seterr(divide='ignore', invalid='ignore')

# ======================================================================================================================
# ===== Filter Curves == Currently set to DES for testing ==============================================================
# ======================================================================================================================

filter_names = ['u', 'g', 'r', 'i', 'z']
filter_curves = {'u': 'DES_u.dat', 'g': 'DES_g.dat', 'r': 'DES_r.dat', 'i': 'DES_i.dat', 'z': 'DES_z.dat'}

# ======================================================================================================================
# ===== Input Parameters ===============================================================================================
# ======================================================================================================================

parser = OptionParser()

parser.add_option("-c", "--config", dest="config", default=None, help="Configuration file")
parser.add_option("-s", "--sens", dest="sens", default=None, help="Sensitivity function")
parser.add_option("-o", "--output", dest="output", default='test.fits', help="Sensitivity function")
parser.add_option("-l", "--limit", dest="limit", default=None, help="Magnitude limit")
parser.add_option("-f", "--filter", dest="filter", default='Savitzky_Golay', help="Filter")
parser.add_option("-a", "--arm", dest="arm", default='spliced', help="Arm to process")
parser.add_option("-z", "--zplimit", dest="zplimit", default=-100.0, help="Zero point limit")
parser.add_option("-n", "--noextinction", dest="noextinction", action="store_true", default=False,
                  help="Extinction Correction")
parser.add_option("-p", "--plot", dest="plot", action="store_true", default=False, help="plot intermediate results")

(options, args) = parser.parse_args()

# ======================================================================================================================
# ===== Classes and Functions ==========================================================================================
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
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    b = numpy.mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
    m = numpy.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
    firstvals = y[0] - numpy.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + numpy.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = numpy.concatenate((firstvals, y, lastvals))

    return numpy.convolve(m[::-1], y, mode='valid')


def correct_extinction(hdr, ext):
    """Correct for Atmospheric Extinction"""

    wave = hdr['CRVAL1'] + (numpy.arange(hdr['NAXIS1']) - hdr['CRPIX1'] + 1) * hdr['CDELT1']
    airmass = 1.0 / numpy.cos((hdr['ZDSTART'] + hdr['ZDEND']) / 2. * numpy.pi / 180.)
    correction = interp1d(ext.wave, 10 ** (ext.extinction * airmass / 2.5))(wave)

    return correction

def getmag(fstar, catalogue):
    """Determine Star's Magnitude"""
    mag = []
    # If there is no measurement, the results is set to -99
    for star in catalogue:
        if star['Candidate Name'] == fstar and star['removed not Fstar'] == '':
            for filter in ['DES <u> SingleEpoch', 'CoAdd mag_auto g', 'CoAdd mag_auto r', 'CoAdd mag_auto i',
                           'CoAdd mag_auto z']:
                if filter in star.keys():
                    if float(star[filter]) > 10.0 and float(star[filter]) < 30.0:
                        mag.append(float(star[filter]))
                    else:
                        mag.append(-99.9)
                else:
                    mag.append(-99.9)
            break  # Some entries are in the catalogue twice, we take the first entry

    return numpy.array(mag)


def computeSens(synthetic, data, var):
    """Compute the Sensitivity Curve with Savitsky Golay"""
    # Divide by the synthetic curve
    # First use a median filter
    # good=~(numpy.isnan(data) | numpy.isnan(var))
    medfilter = 21
    y = medfilt(data / synthetic.flux / 1.e19, medfilter)
    return savitzky_golay(y, 31, 1)


def computeSensGauss(synthetic, data, var):
    """Compute the Sensitivity Curve with Gaussian"""
    # interpolate over pixels with excess variance 
    vmed = numpy.nanmedian(var)
    mask = var > 3 * vmed  # true for values that need to be interpolated over
    xp = numpy.arange(len(data), dtype=int)[~mask]
    fp = data[~mask]
    x = numpy.arange(len(data), dtype=int)[mask]  # indices to interpolate over
    data[mask] = numpy.interp(x, xp, fp)
    return ndimage.gaussian_filter(data / synthetic.flux / 1.e19, sigma=(51), order=0)


def computeSensPoly(synthetic, data, var):
    """Compute the Sensitivity Curve with Polynomial"""
    # High-order polynomials may oscillate wildly
    x = numpy.arange(len(data), dtype=int)
    y = data / synthetic.flux / 1.e19
    p = sp.polyfit(x, y, deg=15)
    return sp.polyval(p, x)


def medfilt(x, k):
    """Apply a length-k median filter to a 1D array x. Boundaries are extended by repeating endpoints"""
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = numpy.zeros((len(x), k), dtype=x.dtype)
    y[:, k2] = x
    for i in range(k2):
        j = k2 - i
        y[j:, i] = x[:-j]
        y[:j, i] = x[0]
        y[:-j, -(i + 1)] = x[j:]
        y[-j:, -(i + 1)] = x[-1]
    return numpy.nanmedian(y, axis=1)


def findTemplate(templates, colour):
    """Find Best Template"""
    selectedTemplate = None
    min = 1.0
    for template in templates:
        diff = abs(float(template['Colour']) - colour)
        if diff < min:
            min = diff
            selectedTemplate = template['SpectralType']
    return selectedTemplate


def warp(fstar, filters, mag, hdr):
    """Warp the spectral template to match the broad band photometry"""
    # First, compute AB magnitudes and effective wavelengths
    magAB = numpy.array([], float)
    effWavelength = numpy.array([], float)
    for filter in filters.keys():
        magAB = numpy.append(magAB, filters[filter].computeABmag(fstar))
        effWavelength = numpy.append(effWavelength, filters[filter].computeEffectiveWavelength(fstar))
    # Take care with the way filters.keys() orders the filters - sort in terms of increasing wavelength
    indices = numpy.argsort(effWavelength)
    magDiff = numpy.array([], float)
    for i in range(len(mag)):
        if mag[i] < 0:
            magDiff = numpy.append(magDiff, 0.0)
        else:
            magDiff = numpy.append(magDiff, mag[i] - magAB[indices[i]])
    s = InterpolatedUnivariateSpline(effWavelength[indices], (magDiff))
    synthetic = spectrum()
    synthetic.wave = hdr['CRVAL1'] + (numpy.arange(hdr['NAXIS1']) - hdr['CRPIX1'] + 1) * hdr['CDELT1']
    ynew = s(synthetic.wave)
    synthetic.flux = interp1d(fstar.wave, fstar.flux)(synthetic.wave) / 10. ** (ynew / 2.5)
    return synthetic


def readData(catalogue, type):
    """Read in catalogue data"""
    global keys, data
    file = open(catalogue)
    listing = file.readlines()
    file.close()
    if type == 'tst' or type == 'tsv':
        firstline = True
        keys = []
        data = []
        data = []
        for line in listing:
            if line[0] not in '#-\n':
                if firstline:
                    firstline = False
                    for key in line.strip().split('\t'):
                        keys.append(key.strip())
                else:
                    elements = line.strip().split('\t')
                    data.append(dict([(keys[i], elements[i].strip()) for i in range(len(elements))]))
    if type == 'csv':
        firstline = True
        keys = []
        data = []
        for line in listing:
            if line[0] not in '#-\n':
                if firstline:
                    firstline = False
                    for key in line.strip().split(','):
                        keys.append(key.strip())
                else:
                    elements = line.strip().split(',')
                    data.append(dict([(keys[i], elements[i].strip()) for i in range(len(elements))]))
    if type == 'NED':
        # find the line with the key
        keyline = False
        keys = []
        data = []
        for line in listing:
            if keyline:
                elements = line.strip().split('\t')
                data.append(dict([(keys[i], elements[i].strip()) for i in range(len(elements))]))
            if len(line.strip()) > 0:
                if line.split()[0].strip() == 'No.':
                    keyline = True
                    for key in line.strip().split('\t'):
                        keys.append(key.strip())
        if len(data) > 50000:
            print 'Warning, request may have been truncated'
    if type == 'vimos':
        firstline = True
        keys = []
        data = []
        for line in listing:
            if line[0] not in '#-d\n':
                if firstline:
                    firstline = False
                    for key in line.strip().split('\t'):
                        keys.append(key.strip())
                else:
                    elements = line.strip().split('\t')
                    data.append(dict([(keys[i], elements[i].strip()) for i in range(len(elements))]))
    elif type == 'txt':
        firstline = True
        keys = []
        data = []
        for line in listing:
            if firstline:
                firstline = False
                for key in line.strip().split():
                    if key != '#':
                        keys.append(key.strip())
            else:
                elements = line.strip().split()
                data.append(dict([(keys[i], elements[i].strip()) for i in range(len(elements))]))
    elif type == 'asciicolumn':
        keys = []
        data = []
        for line in listing:
            if line[0] == '#' and line.split()[1] == "Column":
                keys.append(line.split()[2])
            elif line[0] != '#':
                elements = line.strip().split()
                data.append(dict([(keys[i], elements[i].strip()) for i in range(len(elements))]))
    return data, keys


class spectrum:
    """Generate spectrum"""

    def __init__(self):
        self.wave = numpy.array([], 'float')
        self.flux = numpy.array([], 'float')
        return

    def read(self, template):
        flux, hdr = pyfits.getdata(template, 0, header=True)
        self.flux = flux
        self.wave = (numpy.arange(hdr['NAXIS1']) + 1.0) * hdr['CDELT1'] + hdr['CRVAL1']
        return


class extinction:
    """Generate extinction curve"""

    def __init__(self):
        self.wave = numpy.array([], 'float')
        self.extinction = numpy.array([], 'float')
        return

    def read(self, curve):
        file = open(curve)
        lines = file.readlines()
        file.close()
        for line in lines:
            if line[0] != '#':
                entries = line.split()
                self.wave = numpy.append(self.wave, float(entries[0]))
                self.extinction = numpy.append(self.extinction, float(entries[1]))
        return


class filterCurve:
    """Generate filter curves"""

    def __init__(self):
        self.wave = numpy.array([], 'float')
        self.trans = numpy.array([], 'float')
        self.effectiveWavelength = 0.
        self.magAB = 0.
        return

    def read(self, file):
        # DES filter curves express the wavelengths in nms
        if 'DES' in file:
            factor = 10.
        else:
            factor = 1.
        file = open(file, 'r')
        for line in file.readlines():
            if line[0] != '#':
                entries = line.split()
                self.wave = numpy.append(self.wave, float(entries[0]))
                self.trans = numpy.append(self.trans, float(entries[1]))
        file.close()
        # We use Angstroms for the wavelength in the filter transmission file
        self.wave = self.wave * factor
        return

    def computeEffectiveWavelength(self, template):
        flux1 = interp1d(template.wave, template.flux)(self.wave) * self.trans * self.wave
        flux2 = interp1d(template.wave, template.flux)(self.wave) * self.trans
        self.effectiveWavelength = flux1.sum() / flux2.sum()

        return self.effectiveWavelength

    def computeABmag(self, template):
        # To convert this into f_lambda, divide by c/lambda^2
        c = 2.992792e18  # Angstrom/s
        flux = interp1d(template.wave, template.flux)(self.wave) * self.trans * self.wave ** 2. / 2.992792e18

        # An object with an AB mag of zero has a constant flux density of 3631 Jy
        # 3631 Jy # 1 Jy = 1e-26 W / m^2 / Hz = 3.631 e-20 erg / s /cm^2 / Hz
        # NB 48.60=2.5*log10(3.631e-20)

        const = self.trans
        self.magAB = -2.5 * numpy.log10(flux.sum() / const.sum()) - 48.60

        return self.magAB


def readFilterCurves(filter_names, filterDir, filter_curves):
    filterCurves = {}
    for filter in filter_names:
        filterCurves[filter] = filterCurve()
        filterCurves[filter].read(filterDir + filter_curves[filter])
    return filterCurves


def buildDictionary(file):
    file = open(file, 'r')
    lines = file.readlines()
    file.close()
    parameters = []

    for line in lines:
        if line[0] != '#':
            parameters.append(line.split()[0:2])

    return dict(parameters)


def sens(param, filter_curves, extinct, templates, catalogue, options):
    """Generate Sensitivity Curves"""

    fits = numpy.array([], float)
    nfit = 0
    print options.arm
    directory = param['outputdir'] + '/' + options.arm + '/'
    start = {'red': 5650, 'blue': 3740, 'spliced': 3740}
    end = {'red': 9000, 'blue': 5860, 'spliced': 9000}

    notInUse = 0
    if options.arm == 'blue':
        scale = 0.9
    else:
        scale = 1.0

    for star in os.listdir(directory):
        if 'fits' in star:
            data, hdr = pyfits.getdata(directory + star, 0, header=True)
            var = pyfits.getdata(directory + star, 1, header=False)
            mag = getmag(hdr['OBJECT'], catalogue)
            use = False
            ZPkey = {'red': 'REDZP', 'blue': 'BLUZP'}
            cts = numpy.nanmedian(data[1300:1350])
            if len(mag) > 3:
                if mag[2] > 10. and mag[2] < float(options.limit) and mag[1] > 0 and hdr[ZPkey[options.arm]] < float(
                        options.zplimit):
                    use = True
            if not use:
                notInUse += 1
                print 'Not using %s with ZP %6.2f' % (star, hdr[ZPkey[options.arm]])
            else:
                # Use one of the colours to select the approrpiate F star template
                print 'Using %s with a zeropoint of %6.2f' % (star, hdr[ZPkey[options.arm]])
                template = param['seds'] + '/' + findTemplate(templates, mag[1] - mag[2]) + '.fits'
                fstar = spectrum()
                fstar.read(template)
                # Scale the template to the g magnitude of the star
                # Template spectra are in F_lambda
                fstar.flux = fstar.flux * 3.631e-9 / 10 ** (mag[1] / 2.5)
                # Warp synthetic spectrum to match the observed colours and rebin to wavelength scale of data
                synthetic = warp(fstar, filter_curves, mag, hdr)
                # Compute the extinction correction
                if not options.noextinction:
                    correction = correct_extinction(hdr, extinct)
                    data = data * correction
                    var = var * correction ** 2.
                # Filter the sensitivity function
                if options.filter == 'Savitzky_Golay':
                    # Smooth using the Savitzky Golay algorithm
                    sg = computeSens(synthetic, data, var)
                    wave = numpy.arange(start[options.arm], end[options.arm], scale)
                    fit = interp1d(synthetic.wave, sg, bounds_error=False, fill_value=numpy.nan)(wave)
                elif options.filter == 'Gaussian':  #
                    # Smooth using a Gaussian filter after removing points with excessive variance
                    sg = computeSensGauss(synthetic, data, var)
                    wave = numpy.arange(start[options.arm], end[options.arm], scale)
                    fit = interp1d(synthetic.wave, sg, bounds_error=False, fill_value=numpy.nan)(wave)
                elif options.filter == 'Poly':  #
                    # Smooth using a least squares polynomial fit
                    sg = computeSensPoly(synthetic, data, var)
                    wave = numpy.arange(start[options.arm], end[options.arm], scale)
                    fit = interp1d(synthetic.wave, sg, bounds_error=False, fill_value=numpy.nan)(wave)
                # Scale the fits by the median value
                fits = numpy.append(fits, fit / numpy.nanmedian(fit))
                nfit += 1
    fit = fits.reshape(nfit, len(fits) / nfit)
    print 'Using %d stars' % nfit
    print 'Exlcuded %d stars' % notInUse
    print
    # Compute RMS, excluding nans
    rms = numpy.zeros(len(fits) / nfit, float)
    medianfit = numpy.zeros(len(fits) / nfit, float)
    for i in range(len(fits) / nfit):
        sample = fit[:, i]
        good = ~numpy.isnan(sample)
        y = numpy.sort(sample[good])
        if len(y) > 0:
            index_low = int(0.185 * len(y))
            index_high = len(y) - index_low - 1
            rms[i] = (y[index_high] - y[index_low]) / 2.0
            medianfit[i] = y[int(0.5 * len(y))]
        else:
            rms[i] = numpy.nan
            medianfit[i] = numpy.nan

    if options.plot:
        fig = plt.figure()
        ax = fig.add_subplot(212)
        # plot the average fit and the percentage deviation from the average fit (two plots)
        ax.plot(wave[10:-10], (rms / medianfit)[10:-10], color='navy')
        ax.set_xlabel('Wavelength')
        ax.set_ylim(0, 0.5)
        ax.set_ylabel('scatter')
        ax.set_title('Fractional variation in the sensitivity function')
        ax = fig.add_subplot(211)
        ax.set_ylim(0, 3.0)
        for i in range(nfit):
            ax.set_ylabel('Sensitivity')
            ax.plot(wave[10:-10], fit[i, 10:-10])
        # Compare distribution at some location (5650A + loc) with Gaussian
        fig = plt.figure()
        plots = [{'figure': 311, 'loc': 500}, {'figure': 312, 'loc': 1000}, {'figure': 313, 'loc': 1500}]
        for plot in plots:
            ax = fig.add_subplot(plot['figure'])
            step = 0.01
            loc = plot['loc']
            wl = plot['loc'] + start[options.arm]
            bin = numpy.arange(medianfit[loc] - 0.5, medianfit[loc] + 0.5, step)
            ax.hist(fit[:, loc], bin, range=(0, 2), color='salmon')
            rv = norm(medianfit[loc], rms[loc])
            ax.plot(bin, len(fit[:, loc]) * rv.pdf(bin) * step, color='navy', label='%d' % wl)
            ax.legend()
            ax.set_xlabel('sensitivity')
            #ax.set_xlim([0.4, 1.6])
            ax.set_ylabel('frequency')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(wave, medianfit, label=options.filter + ' filtered', color='navy')
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Sensitivity')
        ax.set_ylim([0, 1.8])
        ax.legend(loc='upper left')
        plt.show()
        plt.close()

    # Write out the sensitivity curve
    hdr['CDELT1'] = scale
    hdr['CRVAL1'] = wave[0]
    hdr['CRPIX1'] = 1.0
    pyfits.writeto(options.output, medianfit, hdr, clobber=True)

    return

# ======================================================================================================================
# ===== Main Code ======================================================================================================
# ======================================================================================================================

param = buildDictionary(options.config)

if options.sens != None:
    # Read in the Fstar catalogue
    catalogue,keys=readData(param['catalogue'],'csv')
    # Read in the filter curves
    filterCurves=readFilterCurves(filter_names,param['filters']+'/',filter_curves)
    # Read in colours and spectral types
    templates,keys=readData(options.sens,'asciicolumn')
    # Read in the extinction curve
    extinction=extinction()
    extinction.read(param['extinction'])
    # Compute sensitivity functions
    sens(param,filterCurves,extinction,templates,catalogue,options)