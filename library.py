import numpy
import pyfits
from scipy.interpolate import interp1d

def readData(catalogue,type):
    global keys, data
    file=open(catalogue)
    listing=file.readlines()
    file.close()

    if type=='tst' or type=='tsv':
        firstline=True
        keys=[]
        data=[]
        data=[]
        for line in listing:
            if line[0] not in '#-\n':
                if firstline:
                    firstline=False
                    for key in line.strip().split('\t'):
                        keys.append(key.strip())
                else:
                    elements=line.strip().split('\t')
                    data.append(dict([(keys[i],elements[i].strip()) for i in range(len(elements))]))

    if type=='csv':
        firstline=True
        keys=[]
        data=[]
        for line in listing:
            if line[0] not in '#-\n':
                if firstline:
                    firstline=False
                    for key in line.strip().split(','):
                        keys.append(key.strip())
                else:
                    elements=line.strip().split(',')
                    data.append(dict([(keys[i],elements[i].strip()) for i in range(len(elements))]))

    if type=='NED':
        # find the line with the key
        keyline=False
        keys=[]
        data=[]
        for line in listing:
            if keyline:
                elements=line.strip().split('\t')
                data.append(dict([(keys[i],elements[i].strip()) for i in range(len(elements))]))
            if len(line.strip()) > 0:
                if line.split()[0].strip() == 'No.':
                    keyline=True
                    for key in line.strip().split('\t'):
                        keys.append(key.strip())
        if len(data) > 50000:
            print 'Warning, request may have been truncated'

    if type=='vimos':
        firstline=True
        keys=[]
        data=[]
        for line in listing:
            if line[0] not in '#-d\n':
                if firstline:
                    firstline=False
                    for key in line.strip().split('\t'):
                        keys.append(key.strip())
                else:
                    elements=line.strip().split('\t')
                    data.append(dict([(keys[i],elements[i].strip()) for i in range(len(elements))]))

    elif type=='txt':
        firstline=True
        keys=[]
        data=[]
        for line in listing:
            if firstline:
                firstline=False
                for key in line.strip().split():
                    if key != '#':
                        keys.append(key.strip())
            else:
                elements=line.strip().split()
                data.append(dict([(keys[i],elements[i].strip()) for i in range(len(elements))]))

    elif type=='asciicolumn':
        keys=[]
        data=[]
        for line in listing:
            if line[0]== '#' and line.split()[1]=="Column":
                keys.append(line.split()[2])
            elif line[0]!='#':
                elements=line.strip().split()
                data.append(dict([(keys[i],elements[i].strip()) for i in range(len(elements))]))

    return data,keys

def buildDictionary(file):
    """????????????"""
    file=open(file,'r')
    lines=file.readlines()
    file.close()
    parameters=[]

    for line in lines:
        if line[0]!='#':
            parameters.append(line.split()[0:2])

    return dict(parameters)

class spectrum:
    """A spectrum"""
    def __init__(self):
        self.wave=numpy.array([],'float')
        self.flux=numpy.array([],'float')

        return

    def read(self,template):
        flux,hdr=pyfits.getdata(template,0,header=True)
        self.flux=flux
        self.wave=(numpy.arange(hdr['NAXIS1'])+1.0)*hdr['CDELT1']+hdr['CRVAL1']

        return

class extinction:
    """A extinction curve"""
    def __init__(self):
        self.wave=numpy.array([],'float')
        self.extinction=numpy.array([],'float')

        return

    def read(self,curve):
        file=open(curve)
        lines=file.readlines()
        file.close()
        for line in lines:
            if line[0]!='#':
                entries=line.split()
                self.wave=numpy.append(self.wave,float(entries[0]))
                self.extinction=numpy.append(self.extinction,float(entries[1]))
        return

class filterCurve:
    """A filter"""
    def __init__(self):
        self.wave=numpy.array([],'float')
        self.trans=numpy.array([],'float')
        self.effectiveWavelength=0.
        self.magAB=0.

        return

    def read(self,file):
        # DES filter curves express the wavelengths in nms
        if 'DES' in file:
            factor=10.
        else:
            factor=1.
        file=open(file,'r')
        for line in file.readlines():
            if line[0]!='#':
                entries=line.split()
                self.wave=numpy.append(self.wave,float(entries[0]))
                self.trans=numpy.append(self.trans,float(entries[1]))
        file.close()
        # We use Angstroms for the wavelength in the filter transmission file
        self.wave=self.wave * factor
        return

    def computeEffectiveWavelength(self,template):
        flux1=interp1d(template.wave,template.flux)(self.wave)*self.trans*self.wave
        flux2=interp1d(template.wave,template.flux)(self.wave)*self.trans
        self.effectiveWavelength=flux1.sum() / flux2.sum()

        return self.effectiveWavelength

    def computeABmag(self,template):
        # To convert this into f_lambda, divide by c/lambda^2
        c=2.992792e18 # Angstrom/s
        flux=interp1d(template.wave,template.flux)(self.wave)*self.trans*self.wave**2./2.992792e18

        # An object with an AB mag of zero has a constant flux density of 3631 Jy
        # 3631 Jy # 1 Jy = 1e-26 W / m^2 / Hz = 3.631 e-20 erg / s /cm^2 / Hz
        # NB 48.60=2.5*log10(3.631e-20)

        const=self.trans
        self.magAB=-2.5*numpy.log10(flux.sum()/const.sum())-48.60

        return self.magAB

def readFilterCurves(filterNames,filterDir,filterTransCurves):
    filterCurves={}
    for filter in filterNames:
        filterCurves[filter]=filterCurve()
        filterCurves[filter].read(filterDir+filterTransCurves[filter])
    return filterCurves

