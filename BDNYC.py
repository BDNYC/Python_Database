#!/usr/bin/env python
# encoding: utf-8

# I ++++++++++++++++++++++ GENERAL DOCUMENTATION ++++++++++++++++++++++++++++++
"""
The module **BDNYC** is a set of classes and methods to handle data from the BDNYC database, which is stored in the BDNYCData.txt file.

**The BDNYC Database**
The structure of the database consists of two types of class instances, with the actual data in the form of nested dictionaries.

**The Database Classes**
An instance of the BDNYCData class is created to hold the whole database. This object holds instances of the Target class in the form of a Python list. The BDNYCData instance also includes methods to handle data from the Target instances.

An instance of the Target class corresponds to a target in the sky. Each target has six attributes that describe it. They are 'name', 'unum', 'ra', 'dec', 'sptype', and 'standard'. The unique identifier is unum. All these are described in detail in the Target class documentation.

Each target also has three dictionaries in the form of attributes, where spectra and photometry are stored. They are 'opt', 'nir', and 'mir'. The keys for each one of these dictionaries are ['high', 'med', 'low', 'phot'], which refer to spectral resolution levels and photometries.

Each one of these keys in turn refer to a dictionary as well. For the 'high', 'med', and 'low' dictionaries, which store spectra, the keys are the instruments.

Each instrument in turn is a dictionary as well. Its keys are the dates of observation. If a date is not known, the key will be '0000xxx00'.

Each observation in turn is a dictionary as well. Its keys are ['wl', 'flux', 'uncertainty', 'snr'].

For the 'phot' dictionary, its keys are the surveys. Each survey in turn is a dictionary as well. Its keys are the bands. Each band in turn is a dictionary as well. Its keys are ['val', 'err'].

For full documentation on the structure tree of the BDNYC database, you can refer to the Google Presentation named 'Python Database Structure'.

**New Targets**
Whenever a new instance of Target is added to the database, two levels of dictionaries are initialized as empty. That is, a target will have the following, regardless of the data that is added:
*target*.opt['high':{}, 'med':{}, 'low':{}, 'phot':{}]
*target*.nir['high':{}, 'med':{}, 'low':{}, 'phot':{}]
*target*.mir['high':{}, 'med':{}, 'low':{}, 'phot':{}]
*target*.phot[{}]


:Authors:
	Dan Feldman, Alejandro N |uacute| |ntilde| ez

:Date of Last Update:
    2012/06/27, Dan

:Repository:
    https://github.com/BDNYC/Python_Database

:Contact: 
    bdnyc.labmanager@gmail.com
     
:Requirements:
    The following modules should already be installed in your computer: `matplotlib`_, `numpy`_.
"""

# II ++++++++++++++++++++++++ EXTERNAL MODULES ++++++++++++++++++++++++++++++++
# External Python modules used by functions and classes

# Basic Python modules
import pdb

# Third party Python modules
import matplotlib.pyplot as plt
import numpy as np

# III +++++++++++++++++++++++ PUBLIC CLASSES ++++++++++++++++++++++++++++++++++
# Classes meant to be used by end users. Capitalize class names.
class Target:
    """
    An object containing data for a specific target object to be put into the BDNYC database. This flavor produces one extra level of nesting.
    
    Parameters:
    
    *name*
      A string with any known name(s) for the target. If there is no common name associated with the target, please use None.
    *unum*
      A string with the U-number of the target. It acts as the unique identifier (e.g. "U10000".)
    *ra*
      A string containing the right ascension of the target (e.g. 13 25 10.3). If unknown, please use placeholders ('XX XX XX.X').
    *dec*
      A string containing the declination of the target, (e.g. +13 25 10.3). If unknown, please use placeholders ('+XX XX XX.X').
    *sptype*
      A string containing the spectral type of the target. It can include detailed descriptors (e.g. "L1.5:b"). If spectral type is unknown, use the placeholder 'XX'.
    *opt*
      A dictionary containing optical spectra and photometry. The keys for this dictionary are described in the general documentation.
    *nir*
      A dictionary containing near-infrared spectra and photometry. The keys for this dictionary are described in the general documentation.
    *mir*
      A dictionary containing mid-infrared spectra and photometry. The keys for this dictionary are described in the general documentation.
    *standard*
      A string identifying whether the target is a standard or a candidate standard. It can be 'Yes' or 'No'.
    """
    
    def __init__(self,name,unum,ra,dec,sptype,opt,nir,mir,standard):
        self.name = name
        self.unum = unum
        self.ra = ra
        self.dec = dec
        self.sptype = sptype
        self.opt = opt
        self.nir = nir
        self.mir = mir
        self.standard = standard
    


class BDNYCData:
    """
    An object that holds all of the targets (i.e. Target instances) —as a Python list— and methods for data handling.
    """
    
    def __init__(self):
        self.targets = []
    
    def addTarget(self, targetObj, init=True):
        """
        Add a new target to the database.
        
        *targetObj*
          The instance of the Target class being added to the database.
        *init*
          Boolean: If True, runs the res_initializer function after the target is added.
        """
        self.targets.append(targetObj)
        if init:
            self.res_initializer()
        return
    
    def dateList(self, obsType, res, surv_instr):
        """
        Sort through the dates and list all targets with observations on given 
        date.
        
        *obsType*
          The spectral range of observation. It can be 'opt', 'nir', or 'mir'.
        *res*
          The resolution of observation. It can be 'low', 'med', or 'high'.
        *surv_instr*
          The name of survey or instrument of observation.
        """
        
        dates = {}
        if obsType != 'nir' and obsType != 'mir' and obsType != 'opt':
            print "Invalid Observation Type. Must be opt, nir, or mir."
            return
        if res != 'low' and res != 'med' and res != 'high':
            print "Invalid Resolution Type. Must be low, med, or high."
            return
        
        # Build the dates dictionary by looping through all existing dates and making keys
        # for each date, and then making a list of all targets observed on those dates  
        for target in self.targets:
            name = target.name
            if obsType == 'nir':
                if surv_instr not in target.nir[res].keys():
                    continue
                tDates = target.nir[res][surv_instr].keys()
            if obsType == 'mir':
                if surv_instr not in target.mir[res].keys():
                    continue
                tDates = target.mir[res][surv_instr].keys()
            if obsType == 'opt':
                if surv_instr not in target.opt[res].keys():
                    continue
                tDates = target.opt[res][surv_instr].keys()
            for d in tDates:
                if d in dates.keys():
                    dates[d].append(name)
                else:
                    dates[d] = [name]
        
        # At this point, dates is a complete dictionary that needs to be sorted, then printed:
        sortedDatesByYear = np.array(sorted(dates.keys(), reverse=True))
        yearsSorted = np.array([], dtype='S4')
        for i in sortedDatesByYear:
            year = i[0:4]
            yearsSorted = np.append(year, yearsSorted)
        yearsUnique = np.array(sorted(np.unique(yearsSorted), reverse=True))
        sortedDatesFull = np.array([], dtype='S9')
        for y in yearsUnique:
            ind = np.where(yearsSorted==y)[0]
            months = np.array([], dtype='S3')
            for num in ind:
                months = np.append(months,sortedDatesByYear[num][4:7])
            monthStrings = ['jan','feb','mar','apr','may','jun','jul','aug','sep',\
            'oct','nov','dec']
            mIndices = np.array([], dtype='I32')
            for m in monthStrings:
                mInd = np.where(months==m)[0]
                mIndices = np.append(mInd, mIndices)
            sortedDatesByMonth = sortedDatesByYear[ind][mIndices]
            sortedDatesFull = np.append(sortedDatesFull,sortedDatesByMonth)
        
        # Time to output the dates and targets to a file for viewing!
        File = open('dateSortOutput.txt', 'w')
        for date in sortedDatesFull:
            File.write(date)
            File.write(':\n')
            targetNames = '\n'.join(dates[date])
            File.write(targetNames)
            File.write('\n\n')
        File.close()
        return
    
    def giveSpectrum(self, unum=None, otype=None, res=None, instr=None, date=None, snr=False, order=None, Filter=None):
        """
        Output a spectrum in a numpy array as [wl, flux, snr(optional)]. The necessary info about the spectrum must be provided as keyword arguments.
        
        *unum*
          U number associated with the target whose spectrum you want (e.g. 'U10000').
        *otype*
          Observation type. It can be 'opt', 'nir', or 'mir'.
        *res*
          Resolution of the spectrum. It can be 'low', 'med', or 'high'.
        *instr*
          Instrument or survey (e.g. NIRSPEC).
        *date*
          Observation date (e.g. 2009jan23).
        *snr*
          Boolean, whether or not to include the snr array.
        *order*
          Integer. If high resolution desired, you must specify the order (e.g. 38).
        *filter*
          If medium or low resolution desired, you must specify the filter (e.g. 'JHK').
        """
        
        specInd = self.matchUNum(unum)
        if otype=='opt':
            
            if res=='high':
                specArr = \
                    np.array([self.targets[specInd].opt['high'] \
                                           [instr][date][order]['wl'], \
                              self.targets[specInd].opt['high'] \
                                           [instr][date][order]['flux']])
                if snr:
                    specArr = np.append(specArr, \
                              [self.targets[specInd].opt['high'] \
                                           [instr][date][order]['snr']], axis=0)
                return specArr
            
            if res=='med':
                specArr = \
                    np.array([self.targets[specInd].opt['med'] \
                                           [instr][date][Filter]['wl'], \
                              self.targets[specInd].opt['med'] \
                                           [instr][date][Filter]['flux']])
                if snr:
                    specArr = np.append(specArr, \
                              [self.targets[specInd].opt['med'] \
                                          [instr][date][Filter]['snr']], axis=0)
                return specArr
            
            if res=='low':
                specArr = \
                    np.array([self.targets[specInd].opt['low'] \
                                           [instr][date]['wl'], \
                              self.targets[specInd].opt['low'] \
                                         [instr][date][Filter]['flux']], axis=0)
                if snr:
                    specArr = np.append(specArr, \
                              [self.targets[specInd].opt['low'] \
                                          [instr][date][Filter]['snr']], axis=0)
                return specArr
        
        elif otype=='nir':
            
            if res=='high':
                specArr = \
                    np.array([self.targets[specInd].nir['high'] \
                                           [instr][date][order]['wl'], \
                              self.targets[specInd].nir['high'] \
                                           [instr][date][order]['flux']])
                if snr:
                    specArr = np.append(specArr, \
                              [self.targets[specInd].nir['high'] \
                                           [instr][date][order]['snr']], axis=0)
                return specArr
            
            if res=='med':
                specArr = \
                    np.array([self.targets[specInd].nir['med'] \
                                           [instr][date][Filter]['wl'], \
                              self.targets[specInd].nir['med'] \
                                           [instr][date][Filter]['flux']])
                if snr:
                    specArr = np.append(specArr, \
                              [self.targets[specInd].nir['med'] \
                                          [instr][date][Filter]['snr']], axis=0)
                return specArr
            
            if res=='low':
                specArr = \
                    np.array([self.targets[specInd].nir['low'] \
                                           [instr][date]['wl'], \
                              self.targets[specInd].nir['low'] \
                                         [instr][date][Filter]['flux']], axis=0)
                if snr:
                    specArr = np.append(specArr, \
                              [self.targets[specInd].nir['low'] \
                                          [instr][date][Filter]['snr']], axis=0)
                return specArr
        
        elif otype=='mir':
            
            if res=='high':
                specArr = \
                    np.array([self.targets[specInd].mir['high'] \
                                           [instr][date][order]['wl'], \
                              self.targets[specInd].mir['high'] \
                                           [instr][date][order]['flux']])
                if snr:
                    specArr = np.append(specArr, \
                              [self.targets[specInd].mir['high'] \
                                           [instr][date][order]['snr']], axis=0)
                return specArr
            
            if res=='med':
                specArr = \
                    np.array([self.targets[specInd].mir['med'] \
                                           [instr][date][Filter]['wl'], \
                              self.targets[specInd].mir['med'] \
                                           [instr][date][Filter]['flux']])
                if snr:
                    specArr = np.append(specArr, \
                              [self.targets[specInd].mir['med'] \
                                          [instr][date][Filter]['snr']], axis=0)
                return specArr
            
            if res=='low':
                specArr = \
                    np.array([self.targets[specInd].mir['low'] \
                                           [instr][date]['wl'], \
                              self.targets[specInd].mir['low'] \
                                         [instr][date][Filter]['flux']], axis=0)
                if snr:
                    specArr = np.append(specArr, \
                              [self.targets[specInd].mir['low'] \
                                          [instr][date][Filter]['snr']], axis=0)
                return specArr
    
    def matchUNum(self, unum, array=False, index=True, verbose=False):
        """
        Match a U number to a target so the code knows which object you are referring to. It returns the array of U numbers, the index of the target object in question, or both, in that order.
        
        *unum*
          U number associated with the target whose spectrum you want (e.g. 'U10000').
        *array*
          Boolean: If True, returns an array containing all of the U numbers in the database.
        *index*
          Boolean: If True, returns the index associated with the given U number.
        *verbose*
          Boolean: If True, prints warning messages.
        """
        
        uNumbers = np.array([], dtype='S6')
        for target in self.targets:
            uNumbers = np.append(uNumbers, target.unum)
        try:
            uNumInd = np.where(uNumbers==unum)[0][0]
        except IndexError:
            if verbose:
                print 'matchUNum: U-number not in database.'
            return
#            raise ValueError('matchUNum: U Number not in database. Could not match')
        if (array==False and index==True):
            return uNumInd
        elif (array==True and index==False):
            return uNumbers
        elif (array==True and index==True):
            return uNumbers, uNumInd
        else:
            if verbose:
                print "matchUNum: You really want to have nothing returned?"
            return
    
    def plotOrders(self, unum, instr, date):
        """
        Plots all of the orders for a given high-resolution spectrum on a single graph.
        
        *unum*
          U number associated with the target whose spectrum you want (e.g. 'U10000').
        *instr*
          Instrument or survey (e.g. NIRSPEC).
        *date*
          Observation date (e.g. 2009jan23).
        """
        
        ind = self.matchUNum(unum)
        for i in self.targets[ind].nir['high'][instr][date].keys():
            plt.plot(self.targets[ind].nir['high'][instr][date][i]['wl'], \
            self.targets[ind].nir['high'][instr][date][i]['flux'])
        return
    
    def res_initializer(self):
        """
        Look at each target in the database and initilizes the keys at the resolution level. So, for instance, if 'med' and 'phot' are the only two resolutions present, it will initialize 'high' and 'low' as empty dictionaries.
        """
        
        for target in self.targets:
            try:
                Buffer = target.opt['low']
            except KeyError:
                target.opt['low'] = {}
            try:
                Buffer = target.opt['med']
            except KeyError:
                target.opt['med'] = {}
            try:
                Buffer = target.opt['high']
            except KeyError:
                target.opt['high'] = {}
            try:
                Buffer = target.opt['phot']
            except KeyError:
                target.opt['phot'] = {}
            
            try:
                Buffer = target.nir['low']
            except KeyError:
                target.nir['low'] = {}
            try:
                Buffer = target.nir['med']
            except KeyError:
                target.nir['med'] = {}
            try:
                Buffer = target.nir['high']
            except KeyError:
                target.nir['high'] = {}
            try:
                Buffer = target.nir['phot']
            except KeyError:
                target.nir['phot'] = {}
            
            try:
                Buffer = target.mir['low']
            except KeyError:
                target.mir['low'] = {}
            try:
                Buffer = target.mir['med']
            except KeyError:
                target.mir['med'] = {}
            try:
                Buffer = target.mir['high']
            except KeyError:
                target.mir['high'] = {}
            try:
                Buffer = target.mir['phot']
            except KeyError:
                target.mir['phot'] = {}
        
        print 'RES INITIALIZER: Database updated but unsaved. ' \
              + 'Please save changes when finished.'
        return
    
    def show_data(self, unum, dump=False):
        '''
        Show all existing data in database for a target.
        
        *unum*
           String with the U-number of target (e.g. 'U10000').
         *dump*
           Boolean, whether to return output on a list. If False, show_data will only print existing data in terminal window.
        '''
        
        # 1. Initialize variables
        resolutions = ['high', 'med', 'low']
        
        # 2. Check for correct input
        try:
            length = len(unum)
        except TypeError:
            print 'U-number must be a string. Try again.'
            return
        
        if length != 6:
            print 'U-number is invalid. Try again.'
            return
        
        # 3. Check if target with requested U-number exists
        idxU = self.matchUNum(unum, verbose=False)
        if idxU is None:
            print 'U-number not in database.'
            return
        
        else:
            curTgt = self.targets[idxU]
        
        # 4. Create structure with all data available for target
        data = []
        
        # 4.1 Add target attributes to output
        data.append(['unum',unum])
        data.append(['index', idxU])
        if curTgt.name is not None:
            data.append(['name', curTgt.name])
        data.append(['sptype', curTgt.sptype])
        data.append(['ra', curTgt.ra])
        data.append(['dec', curTgt.dec])
        data.append(['standard', curTgt.standard])
        
        # 4.2 Add spectra & photometry indicators to output
        # OPT range spectra
        rng = 'opt'
        for res in resolutions:
            specs = curTgt.opt[res]
            if specs != {}:
                for instr in specs.keys():
                    for date in specs[instr]:
                        data.append([rng, res, instr, date])
        # OPT range photometry
        photos = curTgt.opt['phot']
        if photos != {}:
            for survey in photos.keys():
                bands = []
                for band in photos[survey].keys():
                    bands.append(band)
                data.append([rng, 'phot', survey, bands])
        
        # NIR range spectra
        rng = 'nir'
        for res in resolutions:
            specs = curTgt.nir[res]
            if specs != {}:
                for instr in specs.keys():
                    for date in specs[instr]:
                        data.append([rng, res, instr, date])
        # NIR range photometry
        photos = curTgt.nir['phot']
        if photos != {}:
            for survey in photos.keys():
                bands = []
                for band in photos[survey].keys():
                    bands.append(band)
                data.append([rng, 'phot', survey, bands])
        
        # MIR range spectra
        rng = 'mir'
        for res in resolutions:
            specs = curTgt.mir[res]
            if specs != {}:
                for instr in specs.keys():
                    for date in specs[instr]:
                        data.append([rng, res, instr, date])
        # MIR range photometry
        photos = curTgt.mir['phot']
        if photos != {}:
            for survey in photos.keys():
                bands = []
                for band in photos[survey].keys():
                    bands.append(band)
                data.append([rng, 'phot', survey, bands])
        
        # 5. Give output
        if dump:
            return data
        else:
            for row in data:
                print row
            return
    
