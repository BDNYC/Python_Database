#!/usr/bin/env python
# encoding: utf-8
"""
BDNYC.py

Last updated by: Dan 6/14/2012

The code underlying our BDNYC Python Database, based at AMNH.
"""

import sys
import os
import pylab
import numpy as np
import pdb

class Target:
    """
    An object containing data for a specific target object to be put into
    the brown dwarf database. This flavor produces one extra level of nesting.
    """
    
    def __init__(self,name,unum,ra,dec,sptype,opt,nir,mir,standard=None):
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
    Class to hold all of the target objects and methods for Brown Dwarf analysis.
    """
    
    def __init__(self):
        self.targets = []
    
    def addTarget(self, targetObj):
        """
        Adds a new target object to the database.
        """
        self.targets.append(targetObj)
    
    def matchUNum(self, unum, array=False, index=True):
        """
        Matches a U# to a target so the code knows which object you are 
        referring to. This function returns the array of U#s, the index of 
        the target object in question, or both, in that order.
        """
        
        uNumbers = np.array([], dtype='S6')
        for target in self.targets:
            uNumbers = np.append(uNumbers, target.unum)
        try:
            uNumInd = np.where(uNumbers==unum)[0][0]
        except IndexError:
            raise ValueError('matchUNum: U Number not in database. Could not match')
        if (array==False and index==True):
            return uNumInd
        elif (array==True and index==False):
            return uNumbers
        elif (array==True and index==True):
            return uNumbers, uNumInd
        else:
            print "matchUNum: You really want to have nothing returned?"
            return
    
    def plotOrders(self, unum, instr, date):
        ind = self.matchUNum(unum)
        for i in self.targets[ind].nir['high'][instr][date].keys():
            pylab.plot(self.targets[ind].nir['high'][instr][date][i]['wl'], \
            self.targets[ind].nir['high'][instr][date][i]['flux'])
    
    def overPlot(self, unum, med_instr, h_instr, dateM, dateH):
        ind = self.matchUNum(unum)
        for i in self.targets[ind].nir['high'][h_instr][dateH].keys():
            pylab.plot(self.targets[ind].nir['high'][h_instr][dateH][i]['wl'], \
            self.targets[ind].nir['high'][h_instr][dateH][i]['flux'], 'b')
        for j in self.targets[ind].nir['med'][med_instr][dateM].keys():
            pylab.plot(self.targets[ind].nir['med'][med_instr][dateM][j]['wl'], \
            self.targets[ind].nir['med'][med_instr][dateM][j]['flux'], 'r')
    
    def dateList(self, obsType, res, surv_instr):
        """
        Sort through the dates and list all targets with observations on said 
        date. obsType is opt, nir, or mir. res is what type of data (i.e. low, med,
        high). surv_instr is the survey or instrument observed with.
        """
        
        dates = {}
        if obsType != 'nir' and obsType != 'mir' and obsType != 'opt':
            print "Invalid Observation Type. Must be opt, nir, or mir!"
            return
        if res != 'low' and res != 'med' and res != 'high':
            print "Invalid Resolution Type. Must be low, med, or high!"
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
        Function to output a spectrum in a numpy array as [wl, flux, snr(optional)]
        The necessary info about the spectrum must be provided as keyword arguments.
        
        The kwargs are:
        unum   = U number associated with the target whose spectrum you want
        otype  = observation type, i.e. opt, nir, or mir
        res    = the resolution of the spectrum, i.e. low, med, or high
        instr  = the instrument or survey, e.g. NIRSPEC, Spex, Phoenix
        date   = the observation date, given in format YYYYmonDD, e.g. 2009jan23
        snr    = whether or not to include the snr array. This is True or False
        order  = if high resolution desired, then you must specify the order
        filter = if medium resolution is desired, then you must specify the filter.
        """
        
        specInd = self.matchUNum(unum)
        if otype=='opt':
            
            if res=='high':
                specArr = np.array([self.targets[specInd].opt['high'][instr][date]\
                [order]['wl'], self.targets[specInd].opt['high'][instr][date][order]['flux']])
                if snr:
                    specArr = np.append(specArr, [self.targets[specInd].opt['high'][instr]\
                    [date][order]['snr'], axis=0)
                return specArr
            
            if res=='med':
                specArr = np.array([self.targets[specInd].opt['med'][instr][date][Filter]\
                ['wl'], self.targets[specInd].opt['med'][instr][date][Filter]['flux']])
                if snr:
                    specArr = np.append(specArr, [self.targets[specInd].opt['med'][instr]\
                    [date][Filter]['snr']], axis=0)
                return specArr
            
            if res=='low':
                specArr = np.array([self.targets[specInd].opt['low'][instr][date]\
                ['wl'], self.targets[specInd].opt['low'][instr][date]['flux']])
                if snr:
                    specArr = np.append(specArr, [self.targets[specInd].opt['low'][instr]\
                    [date]['snr']], axis=0)
                return specArr
        
        if otype=='nir':
            
            if res=='high':
                specArr = np.array([self.targets[specInd].nir['high'][instr][date]\
                [order]['wl'], self.targets[specInd].nir['high'][instr][date][order]['flux']])
                if snr:
                    specArr = np.append(specArr, [self.targets[specInd].nir['high'][instr]\
                    [date][order]['snr']], axis=0)
                return specArr
            
            if res=='med':
                specArr = np.array([self.targets[specInd].nir['med'][instr][date][Filter]\
                ['wl'], self.targets[specInd].nir['med'][instr][date][Filter]['flux']])
                if snr:
                    specArr = np.append(specArr, [self.targets[specInd].nir['med'][instr]\
                    [date][Filter]['snr']], axis=0)
                return specArr
            
            if res=='low':
                specArr = np.array([self.targets[specInd].nir['low'][instr][date]\
                ['wl'], self.targets[specInd].nir['low'][instr][date]['flux']])
                if snr:
                    specArr = np.append(specArr, [self.targets[specInd].nir['low'][instr]\
                    [date]['snr']], axis=0)
                return specArr
        
        if otype=='mir':
            
            if res=='high':
                specArr = np.array([self.targets[specInd].mir['high'][instr][date]\
                [order]['wl'], self.targets[specInd].mir['high'][instr][date][order]['flux']])
                if snr:
                    specArr = np.append(specArr, [self.targets[specInd].mir['high'][instr]\
                    [date][order]['snr']], axis=0)
                return specArr
            
            if res=='med':
                specArr = np.array([self.targets[specInd].mir['med'][instr][date][Filter]\
                ['wl'], self.targets[specInd].mir['med'][instr][date][Filter]['flux']])
                if snr:
                    specArr = np.append(specArr, [self.targets[specInd].mir['med'][instr]\
                    [date][Filter]['snr']], axis=0)
                return specArr
            
            if res=='low':
                specArr = np.array([self.targets[specInd].mir['low'][instr][date]\
                ['wl'], self.targets[specInd].mir['low'][instr][date]['flux']])
                if snr:
                    specArr = np.append(specArr, [self.targets[specInd].mir['low'][instr]\
                    [date]['snr']], axis=0)
                return specArr
    
