#!/usr/bin/env python
# encoding: utf-8

# I ++++++++++++++++++++++ GENERAL DOCUMENTATION ++++++++++++++++++++++++++++++
"""
The module **BDNYC** is a set of classes and methods to handle data from the BDNYC database, which is stored in the BDNYCData.txt file.

:Authors:
	Dan Feldman, Alejandro N |uacute| |ntilde| ez

:Last Update:
    2012/07/03, Alejo
    
:Repository:
    https://github.com/BDNYC/Python_Database

:Contact:
    bdnyc.labmanager@gmail.com

:Requirements:
    The following modules should already be installed in your computer: `astrotools`_, `asciidata`_, `matplotlib`_, `numpy`_.

**The BDNYC Database**

The structure of the database consists of two types of class instances, with the actual data in the form of nested dictionaries.

**The Database Classes**

An instance of the BDNYCData class is created to hold the whole database. This object holds instances of the Target class in the form of a Python list. The BDNYCData instance also includes methods to handle data from the Target instances.

An instance of the Target class corresponds to a target in the sky. Each target has six attributes that describe it. They are 'name', 'unum', 'ra', 'dec', 'sptype', and 'standard'. The unique identifier is unum. All these are described in detail in the Target class documentation.

Each target also has three dictionaries in the form of attributes, where spectra and photometry are stored. They are 'opt', 'nir', and 'mir'. The keys for each one of these dictionaries are ['high', 'med', 'low', 'phot'], which refer to spectral resolution levels and photometries.

Each one of these keys in turn refer to a dictionary as well. For the 'high', 'med', and 'low' dictionaries, which store spectra, the keys are the instruments.

Each instrument in turn is a dictionary as well. Its keys are the dates of observation. If a date is not known, the key will be '0000xxx00'.

Each date in turn is a dictionary as well. Its keys are the order numbers (for high resolution spectra) or the filter names (for medium and low resolution spectra).

Each order/filter in turn is a dictionary as well. Its keys are ['wl', 'flux', 'uncertainty', 'snr'].

For the 'phot' dictionary, its keys are the surveys. Each survey in turn is a dictionary as well. Its keys are the bands. Each band in turn is a dictionary as well. Its keys are ['val', 'err'].

For full documentation on the structure tree of the BDNYC database, you can refer to the Google Presentation named 'Python Database Structure'.

**New Targets**

Whenever a new instance of Target is added to the database, two levels of dictionaries are initialized as empty. That is, a target will have the following, regardless of the data that is added:

*target*.opt['high':{}, 'med':{}, 'low':{}, 'phot':{}]

*target*.nir['high':{}, 'med':{}, 'low':{}, 'phot':{}]

*target*.mir['high':{}, 'med':{}, 'low':{}, 'phot':{}]
"""

# II ++++++++++++++++++++++++ EXTERNAL MODULES ++++++++++++++++++++++++++++++++
# External Python modules used by functions and classes

# Basic Python modules
import pdb
import pickle

# BDNYC modules
import astrotools as at

# Third party Python modules
import asciidata as ad
import matplotlib.pyplot as plt
import numpy as np

# III +++++++++++++++++++++++ PUBLIC CLASSES ++++++++++++++++++++++++++++++++++
# Classes meant to be used by end users. Capitalize class names.
class Target:
    """
    An object containing data for a specific target in the sky to be put into the BDNYC database. This flavor produces one extra level of nesting.
    
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
    
    def add_data(self, overwrite=False):
        """
        Add spectral and/or photometry data for a target. If target does not exist in database, it will create a Target instance automatically using the *add_target* method. The data to upload is read form the *upload.txt* file located in the same folder as the database file in your computer. Read *upload.txt* header for more information.
        
        *overwrite*
          Boolean, whether to overwrite existing data or not.
        """
        
        # 1. Initialize variables ---------------------------------------------
        DB_FILE = 'BDNYCData.txt'
        UP_FILE = 'upload.txt' # ascii file with upload data
        UP_HEADERS = ('unum','name','ra', 'dec', 'sptype', 'standard', \
                      'rng', 'res', 'instr', 'date', 'ord_filt', 'fitsname', \
                      'survey', 'band_1', 'val_1', 'err_1', \
                      'band_2', 'val_2', 'err_2', \
                      'band_3', 'val_3', 'err_3')
        NULL = ''   # Null character in ascii file
        DEL  = '\t' # Delimiter character in ascii file
        COMM = '#'  # Comment character in ascii file
        NUM_BANDS = 3
        
        colNmUnum = UP_HEADERS[0]
        colNmName = UP_HEADERS[1]
        colNmRa   = UP_HEADERS[2]
        colNmDec  = UP_HEADERS[3]
        colNmSptype = UP_HEADERS[4]
        colNmStd    = UP_HEADERS[5]
        colNmRng    = UP_HEADERS[6]
        colNmRes    = UP_HEADERS[7]
        colNmInstr  = UP_HEADERS[8]
        colNmDate   = UP_HEADERS[9]
        colNmOrdfilt = UP_HEADERS[10]
        colNmFits    = UP_HEADERS[11]
        colNmSurvey  = UP_HEADERS[12]
        colNmBand1 = UP_HEADERS[13]
        colNmVal1  = UP_HEADERS[14]
        colNmErr1  = UP_HEADERS[15]
        colNmBand2 = UP_HEADERS[16]
        colNmVal2  = UP_HEADERS[17]
        colNmErr2  = UP_HEADERS[18]
        colNmBand3 = UP_HEADERS[19]
        colNmVal3  = UP_HEADERS[20]
        colNmErr3  = UP_HEADERS[21]
        
        # 2. Load ascii file --------------------------------------------------
        dataRaw = ad.open(UP_FILE, null=NULL, delimiter=DEL, comment_char=COMM)
        # Store ascii data in a dictionary-type object
        data = {}.fromkeys(UP_HEADERS)
        for colIdx, colData in enumerate(dataRaw):
            data[UP_HEADERS[colIdx]] = colData.tonumpy()
        
        if data[colNmUnum] is None:
            print 'Upload file empty.'
            return
        
        # 3. Upload data to database ------------------------------------------
        somethingAdded = False
        for row in range(len(data[colNmUnum])):
            
            # 3.1 Check if target already exists in database
            unum = data[colNmUnum][row].upper()
            try:
                unum + 0
                print 'U-number invalid.'
                continue
            except TypeError:
                if len(unum) != 6:
                    print 'U-number invalid.'
                    continue
                if unum[0] != 'U':
                    print 'U-number invalid.'
                    continue
            dbIdx = self.match_unum(unum)
            if dbIdx is None:
                newTgt = True
            else:
                newTgt = False
            
            # 3.2 Get target attributes
            if newTgt:
                name = data[colNmName][row]
                if name == '':
                    name = None
                sptype = data[colNmSptype][row].capitalize()
                if sptype == '':
                    sptype = 'XX'
                ra = data[colNmRa][row]
                if len(ra) < 8 or ra == '':
                    print unum + ' Right ascension invalid.'
                    continue
                dec = data[colNmDec][row]
                if len(dec) < 9:
                    print unum + ' Declination invalid.'
                    continue
                std = data[colNmStd][row].capitalize()
                if not (std == 'Yes' or std == 'No'):
                    print unum + ' Standard column must be Yes or No.'
                    continue
            
            # 3.3 Get range of data
            rng = data[colNmRng][row].lower()
            if not (rng == 'opt' or rng == 'nir' or rng == 'mir'):
                print unum + ' Range invalid.'
                continue
            
            # 3.4 Get spectrum data
            # 3.4.1 Open & read fits file
            specAdd = False
            fitsName = data[colNmFits][row]
            if fitsName != '':
                fitsRaw = at.read_spec(fitsName, errors=True)
                if fitsRaw[0] is not None:
                    wl   = fitsRaw[0][0]
                    flux = fitsRaw[0][1]                                                                                                                                                                                                                                                                                  
                    
                    # 3.4.2 Determine if 3rd dimension is uncertainty or snr
                    errNm = None
                    if len(fitsRaw[0]) == 3:
                        errVals = None
                        med = np.median(flux/fitsRaw[0][2])
                        if med < 10**3 and med > 10**-3:
                            errVals = fitsRaw[0][2]
                            errNm = 'uncertainty'
                        else:
                            errVals = fitsRaw[0][2]
                            errNm = 'snr'
                    
                    # 3.4.3 Get spectrum attributes
                    attsOK = True
                    res = data[colNmRes][row].lower()
                    if not (res == 'high' or res == 'med' or res == 'low'):
                        print unum + ' Resolution invalid.'
                        attsOK = False
                    instr = data[colNmInstr][row]
                    if instr == '':
                        print unum + ' Must provide Instrument for spectrum.'
                        attsOK = False
                    date  = data[colNmDate][row].lower()
                    if len(date) != 9:
                        print unum + ' Date invalid.'
                        attsOK = False
                    if date[-1] == '\r':
                        date = date[:-1]
                    ord_filt = data[colNmOrdfilt][row]
                    if ord_filt == '':
                        print unum + ' Must provide Order/Filter for spectrum.'
                        attsOK = False
                    
                    # 3.4.4 Create dictionary structure with spectrum data
                    if attsOK:
                        specAdd = True
                        specDict = {instr:{date:{ord_filt:{'wl':wl, \
                                                'flux':flux, errNm:errVals}}}}
            
            # 3.5 Get target photometry attributes
            # 3.5.1 Check if photometry data was provided
            photAdd = False
            try:
                survey = data[colNmSurvey][row]
            except TypeError:
                survey = ''
            # 3.5.2 Get data for bands
            if survey != '':
                bands = []
                vals  = []
                errs  = []
                photOK = True
                
                for bndNum in range(1,NUM_BANDS + 1):
                    bndNm = 'band_' + str(bndNum)
                    valNm = 'val_' + str(bndNum)
                    errNm = 'err_' + str(bndNum)
                    try:
                        # Get band name
                        bands.append(data[bndNm][row])
                        if bands[bndNum - 1] == '':
                            photOK = False
                        else:
                            try:
                                bands[bndNum - 1][0]
                            except IndexError:
                                photOK = False
                        # Get band photometry value
                        try:
                            vals.append(data[valNm][row])
                            try:
                                vals[bndNum - 1] + 0
                            except TypeError:
                                if vals[bndNum - 1][-1] == '\r':
                                    vals[bndNum - 1] = vals[bndNum - 1][:-1]
                                    try:
                                        vals[bndNum - 1] + 0
                                    except TypeError:
                                        photOK = False
                                else:
                                    photOK = False
                        except TypeError:
                            photOK = False
                        # Get band photometry value error
                        try:
                            errs.append(data[errNm][row])
                            if errs[bndNum - 1] == '':
                                errs[bndNum - 1] = None
                            else:
                                try:
                                    errs[bndNum - 1] + 0
                                except TypeError:
                                    if errs[bndNum - 1][0] == '\r':
                                        errs[bndNum - 1] = None
                                    else:
                                        photOK = False
                        except TypeError:
                            errs.append(None)
                    except TypeError:
                        if bndNum == 1:
                            photOK = False
                
                # 3.5.3 Create dictionary structure with photometry data
                if photOK:
                    photDict = {survey:{}}
                    for bdIdx, band in enumerate(bands):
                        if band != '':
                            photAdd = True
                            photDict[survey][band] = {'val': vals[bdIdx], \
                                                      'err': errs[bdIdx]}                    
                else:
                    print unum + ' Photometry data invalid.'
            
           # 3.6 Create range-level dictionary with all data
            if photAdd and specAdd:
                rngDict = {res:specDict, 'phot':photDict}
            elif photAdd and not specAdd:
                rngDict = {'phot':photDict}
            elif not photAdd and specAdd:
                rngDict = {res:specDict}
            else:
                print 'No data to add for ' + unum
                continue
            
            # 3.7 Add new target to database if necessary
            if newTgt:
                # 3.7.1 Create Target instance
                if rng == 'opt':
                    target = Target(name, unum, ra, dec, sptype, \
                                         rngDict, {}, {}, std)
                elif rng == 'nir':
                    target = Target(name, unum, ra, dec, sptype, \
                                         {}, rngDict, {}, std)
                elif rng == 'mir':
                    target = Target(name, unum, ra, dec, sptype, \
                                         {}, {}, rngDict, std)
                
                # 3.7.2 Add to database
                self.add_target(target, verbose=False)
                somethingAdded = True
                print unum + ' new target added to database: Index # ' + \
                      str(len(self.targets) - 1)
            
            # 3.8 Add new data to database if target already exists
            else:
                addedPhot = False
                addedSpec = False
                if rng == 'opt':
                    currentTgt = self.targets[dbIdx].opt
                elif rng == 'nir':
                    currentTgt = self.targets[dbIdx].nir
                elif rng == 'mir':
                    currentTgt = self.targets[dbIdx].mir
                
                # 3.8.1 Check dictionary level where to add data
                if currentTgt != {}:
                    # For photometry data
                    if photAdd:
                        try:
                            currentTgt['phot'][survey]
                            for bnd in rngDict['phot'][survey].keys():
                                try:
                                    currentTgt['phot'][survey][bnd]
                                    # Overwrite existing data if requested
                                    if overwrite:
                                        addedPhot = True
                                        currentTgt['phot'][survey][bnd] = \
                                                   rngDict['phot'][survey][bnd]
                                    else:
                                        print unum  + ' ' + bnd + \
                                              ' photometry already exists' + \
                                              ' in database.'
                                except KeyError:
                                    addedPhot = True
                                    currentTgt['phot'][survey][bnd] = \
                                                   rngDict['phot'][survey][bnd]
                        except KeyError:
                            addedPhot = True
                            currentTgt['phot'][survey] = \
                                                    rngDict['phot'][survey]
                    # For spectrum data
                    if specAdd:
                        try:
                            currentTgt[res][instr]
                            try:
                                currentTgt[res][instr][date]
                                try:
                                    currentTgt[res][instr][date][ord_filt]
                                    # Overwrite existing data if requested
                                    if overwrite:
                                        addedSpec = True
                                        currentTgt[res][instr][date][ord_filt] \
                                           = rngDict[res][instr][date][ord_filt]
                                    else:
                                        print unum + ', ' + \
                                        data[colNmFits][row] + \
                                        ' spectrum already exists in database.'
                                except KeyError:
                                    addedSpec = True
                                    currentTgt[res][instr][date][ord_filt] = \
                                            rngDict[res][instr][date][ord_filt]
                            except KeyError:
                                addedSpec = True
                                currentTgt[res][instr][date] = \
                                                    rngDict[res][instr][date]
                        except KeyError:
                            addedSpec = True
                            currentTgt[res][instr] = rngDict[res][instr]
                else:
                    addedSpec = True
                    currentTgt = rngDict
                
                if addedSpec or addedPhot:
                    somethingAdded = True
                    print unum + ' new data added to target in database.'
        
        # 4. Commit additions to database file --------------------------------
        if somethingAdded:
            # Check that database txt file exists in current folder
            try:
                f = open(DB_FILE,'rb')
            except IOError:
                print DB_FILE + ' could not be loaded. Check that it is ' + \
                      'in the current folder. Process stopped.'
                return
            f.close()
            f = open(DB_FILE,'wb')
            print 'Updating ' + DB_FILE + '...'
            pickle.dump(self, f)
            f.close()
            
            print 'Remember to push updated ' + DB_FILE + ' to github.'
        
        return
    
    def add_target(self, targetObj, init=True, verbose=True):
        """
        Add a new Target instance to the database. The Target instance must have the correct structure, as illustrated in the Python Database Structure' Google Presentation.
        
        *targetObj*
          The Target instance to add to the database.
        *init*
          Boolean: If True, runs the *res_initializer* method after the target is added.
        *verbose*
          Boolean, whether to print messages.
        """
        self.targets.append(targetObj)
        if init:
            self.res_initializer(verbose=verbose)
        return
    
    def date_list(self, obsType, res, surv_instr):
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
    
    def find_unum(self, ra=None, dec=None, name=None, dump=False):
        '''
        Find potential target matches based on ra, dec, or name. Only one of these three inputs is needed. If more than one is provided, *find_unum* will return targets that match either one. It returns the attributes of all potentially matched targets (unum, index, name, sptype, ra, dec, and standard).
        
        *ra*
          String with the right ascension of the target to be found. It can be whole or part (e.g. 13 25 10.3, or 13 25).
        *dec*
          String with the declination of the target to be found. It can be whole or part (e.g. +13 25 10.3, or +13 25).
        *name*
          String with a known name of the target. it can be whole or part.
        *dump*
          Boolean, whether to return output as a Python list. If False, *find_unum* will only print potential matches in terminal window.
        '''
        if ra is None and dec is None and name is None:
            print 'No input provided.'
            return
        
        matches = []
        for tgt in self.targets:
            # Match by name
            if name is not None and tgt.name is not None and tgt.name != '':
                nmMatch = name.upper()
                loc = tgt.name.upper().find(nmMatch)
                if loc != -1:
                    matches.append(tgt.unum)
                else:
                    loc = nmMatch.upper().find(tgt.name.upper())
                    if loc != -1:
                        matches.append(tgt.unum)
            
            # Match by ra
            if ra is not None:
                raMatch = ra.strip()
                if tgt.ra.startswith(raMatch) or raMatch.startswith(tgt.ra):
                    matches.append(tgt.unum)
            
            # Match by dec
            if dec is not None:
                decMatch = dec.strip()
                if tgt.dec.startswith(decMatch) or decMatch.startswith(tgt.dec):
                    matches.append(tgt.unum)
            
        # Clean up match results
        if matches != []:
            matches = list(set(matches)) # Eliminates duplicates
        else:
            print 'No matches found.'
            return
        
        # Fetch matches attributes
        matched = []
        counter = 1
        for umatch in matches:
            matchAtts = self.show_data(umatch, dump=True)
            matchAtts = matchAtts[:7]
            if dump:
                matched.append(matchAtts)
            else:
                print '-----Potential match ' + str(counter) + ': -----'
                counter = counter + 1
                for row in matchAtts:
                    print row
                print ' '
        
        return matched
    
    def get_data(self, unum, ids, errors=True, header=False):
        """
        Return target data from database, specified using the data ids displayed in the output of the *show_data* method. This method acts as a centralized station to get either spectra or photometry of a target, as opposed to using *give_spectrum* or *give_photometry* methods directly.
        
        *unum*
          U-number of target whose data you want (e.g. 'U10000').
        *ids*
          The data id number as integer or list of integers. These ids are displayed in the output of the *show_data* as the first item in each data row. If more than one id, get_data will return the corresponding data as a list in ascending id order, regardless of the order in the input.
        *errors*
          Boolean, whether to include spectral error values in output.
        *header*
          Boolean, whether to include the data identifiers in the output, that is, its range, resolution, instrument, date, and order/filter.
        """
        
        # Initialize variables
        INIT_ROWS = 7
        data = []
        
        # Convert to list if input is only one id
        try:
            ids[0]
        except TypeError:
            ids = [ids]
        ids = sorted(ids)
        
        # Get data headers in database for target
        dataIdent = self.show_data(unum, dump=True)
        if dataIdent is None:
            return
        
        # Match id requests with available identifiers
        for idnum in ids:
            matched = False
            for rowIdent in dataIdent[INIT_ROWS:]:
                if idnum == rowIdent[0]:
                    matched = True
                    rng = rowIdent[1]
                    res = rowIdent[2]
                    if res == 'phot':
                        survey = rowIdent[3]
                    else:
                        instr = rowIdent[3]
                        date  = rowIdent[4]
                        if res == 'high':
                            order = rowIdent[5]
                            filt = None
                        elif res == 'med' or res == 'low':
                            order = None
                            filt  = rowIdent[5]
                    
                    # Get requested data
                    if res == 'phot':
                        photom = self.give_photometry(unum, rng, survey)
                        if header:
                            data.append([photom, rowIdent])
                        else:
                            data.append(photom)
                    else:
                        spec = self.give_spectrum(unum, rng, res, instr, \
                                                 date, errors, order, filt)
                        if header:
                            data.append([spec, rowIdent])
                        else:
                            data.append(spec)
            if not matched:
                print 'Id ' + str(idnum) + ' not valid.'
                data.append([])
        
        if len(ids) == 1:
            return data[0]
        else:
            return data
    
    def give_photometry(self, unum, rng, survey):
        """
        Return photometry of target as a dictionary for the range and survey requested.
        
        *unum*
          U-number of target whose photometry you want (e.g. 'U10000').
        *rng*
          Observation range. It can be 'opt', 'nir', or 'mir'.
        *survey*
          Survey that measured the photometry. (e.g. '2MASS').
        """
        
        # 1. Get the target's index in the database
        specIdx = self.match_unum(unum)
        if specIdx is None:
            print str(unum) + ' not in database.'
            return
        
        # 2. Determine data range
        if rng == 'opt':
            specDict = self.targets[specIdx].opt
        elif rng == 'nir':
            specDict = self.targets[specIdx].nir
        elif rng == 'mir':
            specDict = self.targets[specIdx].mir
        else:
            print '"' + str(rng) + '" Range invalid.'
            return
        
        # 3. Get requested photometry
        try:
            photom = specDict['phot'][survey]
        except KeyError:
            print 'Photometry requested not found.'
            return
        
        return photom
    
    def give_spectrum(self, unum, rng, res, instr, date, errors=True, order=None, Filter=None):
        """
        Return spectrum of target in a numpy array as [wl, flux, errors(optional)].
        
        *unum*
          U-number of target whose spectrum you want (e.g. 'U10000').
        *rng*
          Observation range. It can be 'opt', 'nir', or 'mir'.
        *res*
          Resolution of the spectrum. It can be 'low', 'med', or 'high'.
        *instr*
          Instrument or survey (e.g. 'NIRSPEC').
        *date*
          Observation date (e.g. '2009jan23').
        *errors*
          Boolean, whether to include flux errors (snr or uncertainties).
        *order*
          Integer. If high resolution requested, you must specify the order (e.g. 38).
        *Filter*
          If medium or low resolution requested, you must specify the filter (e.g. 'JHK').
        """
        
        # 1. Get the target's index in the database
        specIdx = self.match_unum(unum)
        if specIdx is None:
            print str(unum) + ' not in database.'
            return
        
        # 2. Do basic input checks
        if len(date) != 9:
            print '"' + str(date) + '" date invalid.'
            return
        if res == 'high':
            try:
                order + 0
            except TypeError:
                print 'Order must be an integer.'
                return
        else:
            if Filter is None:
                print 'Filter needed for low and medium resolution spectra.'
                return
        
        # 3. Use relevant key (order or filter) according to resolution
        if res == 'high':
            ord_filt = order
        elif res == 'med' or res == 'low':
            ord_filt = Filter
        else:
            print '"' + str(res) + '" Resolution invalid.'
            return
        
        # 4. Determine spectrum range
        if rng == 'opt':
            specDict = self.targets[specIdx].opt
        elif rng == 'nir':
            specDict = self.targets[specIdx].nir
        elif rng == 'mir':
            specDict = self.targets[specIdx].mir
        else:
            print '"' + str(rng) + '" Range invalid.'
            return
        
        # 5. Get spectrum
        try:
            wl = specDict[res][instr][date][ord_filt]['wl']
        except KeyError:
            print 'Data requested not found.'
            return
        flux = specDict[res][instr][date][ord_filt]['flux']
        specArr = np.array([wl, flux])
        
        # 6. Include flux errors if requested and available
        if errors:
            try:
                errArr = specDict[res][instr][date][ord_filt]['uncertainty']
            except KeyError:
                try:
                    errArr = specDict[res][instr][date][ord_filt]['snr']
                except KeyError:
                    errArr = None
                    print 'No flux errors array available for' + rng + ' ' \
                          + res + ' ' + instr + ' ' + date + ' ' + str(ord_filt)
            
            if errArr is not None:
                specArr = np.append(specArr, [errArr], axis=0)
        
        return specArr
    
    def match_unum(self, unum, array=False, index=True, verbose=False):
        """
        Find the index of a target in the database using its U-number. It has the option of returning all U-numbers in the database as a numpy array. If target U-number is not found in database, it returns *None* value.
        
        *unum*
          U-number associated with the target whose spectrum you want (e.g. 'U10000').
        *array*
          Boolean, whether to return an array containing all of the U-numbers in the database.
        *index*
          Boolean, whether to return the index associated with the given U-number.
        *verbose*
          Boolean, whether messages.
        """
        
        uNumbers = np.array([], dtype='S6')
        for target in self.targets:
            uNumbers = np.append(uNumbers, target.unum)
        try:
            uNumInd = np.where(uNumbers==unum)[0][0]
        except IndexError:
            if verbose:
                print 'match_num: U-number not in database.'
            return
#            raise ValueError('match_unum: U Number not in database. Could not match')
        if (array==False and index==True):
            return uNumInd
        elif (array==True and index==False):
            return uNumbers
        elif (array==True and index==True):
            return uNumbers, uNumInd
        else:
            if verbose:
                print "match_unum: You really want to have nothing returned?"
            return
    
    def plot_orders(self, unum, instr, date):
        """
        Plots all of the orders for a given high-resolution spectrum on a single graph.
        
        *unum*
          U-number associated with the target whose spectrum you want (e.g. 'U10000').
        *instr*
          Instrument or survey (e.g. NIRSPEC).
        *date*
          Observation date (e.g. 2009jan23).
        """
        
        ind = self.match_unum(unum)
        for i in self.targets[ind].nir['high'][instr][date].keys():
            plt.plot(self.targets[ind].nir['high'][instr][date][i]['wl'], \
            self.targets[ind].nir['high'][instr][date][i]['flux'])
        return
    
    def res_initializer(self, verbose=True):
        """
        Scans all Target instances in the database looking for missing *range* or *resolution* level dictionaries. When missing, it will create them as empty dictionaries. So, for instance, if 'med' and 'phot' are the only two resolutions present on a Target instance, it will initialize 'high' and 'low' as empty dictionaries.
        
        *verbose*
          Boolean, whether to print messages.
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
        
        if verbose:
            print 'RES INITIALIZER: Database updated but unsaved. ' \
                + 'Please save changes when finished.'
        return
    
    def show_data(self, unum, dump=False):
        '''
        Show all existing data identifiers in database for a target. Each piece of data gets assigned an id number, which is displayed as the first item in each data row. These ids can be used in the *get_data* method to get the actual data from the database.
        
        *unum*
           U-number of target (e.g. 'U10000').
        *dump*
          Boolean, whether to return output as a Python list. If False, *show_data* will only print existing data identifiers in terminal window.
        '''
        
        # 1. Initialize variables
        resolutions = ['high', 'med', 'low']
        
        # 2. Check for correct input
        try:
            length = len(unum)
        except TypeError:
            print 'U-number must be a string.'
            return
        
        if length != 6:
            print 'U-number invalid.'
            return
        
        # 3. Check if target with requested U-number exists
        idxU = self.match_unum(unum, verbose=False)
        if idxU is None:
            print 'U-number not in database.'
            return
        
        else:
            curTgt = self.targets[idxU]
        
        # 4. Create structure with all data identifiers available for target
        data = []
        
        # 4.1 Add target attributes to output
        data.append(['unum',unum])
        data.append(['index', idxU])
        if curTgt.name is not None:
            curNm = curTgt.name
        else:
            curNm = ''
        data.append(['name', curNm])
        data.append(['sptype', curTgt.sptype])
        data.append(['ra', curTgt.ra])
        data.append(['dec', curTgt.dec])
        data.append(['standard', curTgt.standard])
        
        # 4.2 Add spectra & photometry identifiers to output
        count = 0
        # OPT range spectra
        rng = 'opt'
        for res in resolutions:
            specs = curTgt.opt[res]
            if specs != {}:
                for instr in specs.keys():
                    for date in specs[instr]:
                        for order in specs[instr][date]:
                            data.append([count, rng, res, instr, date, order])
                            count = count + 1
        # OPT range photometry
        photos = curTgt.opt['phot']
        if photos != {}:
            for survey in photos.keys():
                bands = []
                for band in photos[survey].keys():
                    bands.append(band)
                data.append([count, rng, 'phot', survey, bands])
                count = count + 1
        
        # NIR range spectra
        rng = 'nir'
        for res in resolutions:
            specs = curTgt.nir[res]
            if specs != {}:
                for instr in specs.keys():
                    for date in specs[instr]:
                        for order in specs[instr][date]:
                            data.append([count, rng, res, instr, date, order])
                            count = count + 1
        # NIR range photometry
        photos = curTgt.nir['phot']
        if photos != {}:
            for survey in photos.keys():
                bands = []
                for band in photos[survey].keys():
                    bands.append(band)
                data.append([count, rng, 'phot', survey, bands])
                count = count + 1
        
        # MIR range spectra
        rng = 'mir'
        for res in resolutions:
            specs = curTgt.mir[res]
            if specs != {}:
                for instr in specs.keys():
                    for date in specs[instr]:
                        for order in specs[instr][date]:
                            data.append([count, rng, res, instr, date, order])
                            count = count + 1
        # MIR range photometry
        photos = curTgt.mir['phot']
        if photos != {}:
            for survey in photos.keys():
                bands = []
                for band in photos[survey].keys():
                    bands.append(band)
                data.append([count, rng, 'phot', survey, bands])
                count = count + 1
        
        # 5. Determine output display
        if dump:
            return data
        else:
            for row in data:
                print row
            return
    
