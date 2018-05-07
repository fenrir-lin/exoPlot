#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 15:21:24 2017

@author: zifan
"""

import matplotlib.patches as mpatches
from matplotlib import pyplot as plt


##############################################################################
###### CONSTANTS #############################################################
##############################################################################
DETECTION_METHODS = ["astrometry", "microlensing", "pulsar", "transit", "imaging",
                     "other", "radialVelocity", "TTV"]
RA = 69  # index of right ascension
DEC = 70  # index of declination
MASS = 2  # index of planet mass
RADIUS = 8  # index of planet radius
DISTANCE = 76  # index of star distance
##############################################################################
MJ = 1.898E27  # mass of Jupiter in kg
ME = 5.972E24  # mass of Earth in kg
RJ = 69911.0  # radius of Jupiter in km
RE = 6371.0  # radius of Earth in km
MASS_CORRECTION_FACTOR = MJ/ME
RADIUS_CORRECTION_FACTOR = RJ/RE
##############################################################################
RED = "#FF0000"
YELLOW = "#FAD80A"
LIME = "#00FF00"
AQUA = "#00FFFF"
BLUE = "#0000FF"
MAROON = "#800000"
PURPLE = "#800080"
OLIVE = "#808000"
COLOR_LIST = [RED, YELLOW, LIME, AQUA, BLUE, MAROON, PURPLE, OLIVE]
##############################################################################
abnormal_star_list = ['55 Cnc','EPIC 211391664','HD 3167','HD 3167','EPIC 211089792',
                    'K2-97','KELT-16','HD 185603','Kepler-1651 A','Kepler-1319 A',
                    'Kepler-1649','KIC 9663113','Kepler-459','KOI-3791','KOI-3791',
                    'Kepler-64 (AB)','Kepler-86','WASP-103','WASP-47','WASP-98',
                    'XO-1','XO-2N','2M 2236+4751','51 Eri','51 Peg',' HD 1237 b',      
                    'HD 179949','HD 59686 A','YZ Cet']
##############################################################################
##############################################################################
##############################################################################


##############################################################################
###### FUNCTIONS #############################################################
##############################################################################
def readIn():
    inputDict = {}
    for method in DETECTION_METHODS:
        with open ("DetectionMethod/" + method + ".csv", "r") as infile:
            inputDict[method] = infile.readlines()
            for i in range(len(inputDict[method])):
                inputDict[method][i] = inputDict[method][i].split(",")
    return inputDict


def findIndex():
    inputDict = readIn()
    titleList = inputDict["pulsar"][0]
    print titleList.index("ra")
    print titleList.index("dec")
    print titleList.index("mass")
    print titleList.index("radius")
    print titleList.index("star_distance")
    print titleList


def exoPlot(dataDict,title):
    """ Parameters:
            dataDict: a dictionary with 8 keys, each corresbond to one detection
        method. The value of each key is a list consisted of two lists, one
        containing star right ascension, one containing declination.
            title: title of the output plot
    """
    color_dict = {'astrometry': PURPLE,
                  'microlensing': LIME,
                  'pulsar': OLIVE,
                  'transit': RED,
                  'imaging':YELLOW,
                  'other': AQUA,
                  'radialVelocity': BLUE,
                  'TTV': MAROON
                  }
    count = 0
    plt.figure(figsize=(65,30))  # figure size to be determined
    plt.title(title, fontsize = 50)
    plt.xlabel("Right ascension", fontsize = 45)
    plt.ylabel("Declination", fontsize = 45)
    plt.axis([0,360,-90,90])
    plt.xticks([0,30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360], fontsize = 40)
    plt.yticks([-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90], fontsize = 40)
    ### define the color patches indicating detection methods
    purple_patch = mpatches.Patch(color=PURPLE, label='Astrometry')
    lime_patch = mpatches.Patch(color=LIME, label='Microlensing')
    olive_patch = mpatches.Patch(color=OLIVE, label='Pulsar')
    red_patch = mpatches.Patch(color=RED, label='Transit')
    yellow_patch = mpatches.Patch(color=YELLOW, label='Imaging')
    aqua_patch = mpatches.Patch(color=AQUA, label='Other')
    blue_patch = mpatches.Patch(color=BLUE, label='Radial Velocity')
    maroon_patch = mpatches.Patch(color=MAROON, label='TTV')
    plt.legend(handles=[purple_patch,lime_patch,olive_patch,red_patch,
                           yellow_patch,aqua_patch,blue_patch,maroon_patch],
               loc='center left', bbox_to_anchor=(1, 0.5),fontsize=40)
    #leng = 0
    for method in dataDict:
        clr = color_dict[method]
        coordinateList  = dataDict[method]
        #leng += len(coordinateList[0])
        plt.plot(coordinateList[0], coordinateList[1], 'o',
                 markersize = 20, color = clr)
        count += 1
    filename = "plot/" + title + '.eps'
    plt.savefig(filename)
    #print leng


def getRaDec(inputDict):
    """ Reads in a dictionary in the style of output dictionary of readIn,
    and extracts the right ascension and declination information from it. """
    leng = 0
    outputDict = {"astrometry":[[],[]],
                  "microlensing":[[],[]],
                  "pulsar":[[],[]],
                  "transit":[[],[]],
                  "imaging":[[],[]],
                  "other":[[],[]],
                  "radialVelocity":[[],[]],
                  "TTV":[[],[]]
                  }
    for method in DETECTION_METHODS:
        infoList = inputDict[method]
        for planetInfo in infoList:
            if planetInfo[0] != "name":
                try:
                    ra = planetInfo[RA]
                    dec = planetInfo[DEC]
                    if ra in abnormal_star_list:
                        ra = planetInfo[RA+1]
                        dec = planetInfo[DEC+1]
                    ra = float(ra)
                    ra = 360.0 - ra  # the data given by exoplanet.eu have reversed right ascension
                    dec = float(dec)
                    assert ra < 360 and ra > 0
                    assert dec > -90 and dec < 90
                    outputDict[method][0].append(ra)
                    outputDict[method][1].append(dec)
                    leng += 1
                except:
                    pass
#                    print planetInfo[RA]
    print "all " + str(leng)
    return outputDict


def getEarthLike(inputDict):
    """ Returns a dictionary containing the information of all earth-like exoplanets.
    The definition of earth-like planet in this case is:
        (1) R < 2Re, and
        (2) M < 10Me
    """
    leng = 0
    outputDict = {"astrometry":[[],[]],
              "microlensing":[[],[]],
              "pulsar":[[],[]],
              "transit":[[],[]],
              "imaging":[[],[]],
              "other":[[],[]],
              "radialVelocity":[[],[]],
              "TTV":[[],[]]
              }
    for method in DETECTION_METHODS:
        infoList = inputDict[method]
        for planetInfo in infoList:
            if planetInfo[0] != "name":
                mass = planetInfo[MASS]
                radius = planetInfo[RADIUS]
                ra = planetInfo[RA]
                dec = planetInfo[DEC]
                try:
                    if mass != '':
                        mass = float(mass)
                        mass = mass * MASS_CORRECTION_FACTOR
                    if radius != '':
                        radius = float(radius)
                        radius = radius * RADIUS_CORRECTION_FACTOR
                    if ra in abnormal_star_list:
                        ra = planetInfo[RA+1]
                        dec = planetInfo[DEC+1]
                    if ra != '' and dec != '':
                        ra = float(ra)
                        ra = 360.0 - ra  # the data given by exoplanet.eu have reversed right ascension
                        dec = float(dec)
                    assert ra < 360 and ra > 0
                    assert dec > -90 and dec < 90
                except:
                    pass
#                    print "conversion error!"
#                    print mass
#                    print radius
#                    print ra
#                    print dec
                if (radius != ''  and radius < 2.0) or (mass != '' and mass < 10.0) \
                and ra != '' and dec != '':
                    outputDict[method][0].append(ra)
                    outputDict[method][1].append(dec)
                    leng += 1
    print "earth size " + str(leng)
    return outputDict


def getWithin10pc(inputDict):
    leng = 0
    outputDict = {"astrometry":[[],[]],
                  "microlensing":[[],[]],
                  "pulsar":[[],[]],
                  "transit":[[],[]],
                  "imaging":[[],[]],
                  "other":[[],[]],
                  "radialVelocity":[[],[]],
                  "TTV":[[],[]]
                  }
    for method in DETECTION_METHODS:
        infoList = inputDict[method]
        for planetInfo in infoList:
            if planetInfo[0] != "name":
                try:
                    ra = planetInfo[RA]
                    dec = planetInfo[DEC]
                    distance = planetInfo[DISTANCE]
                    if ra in abnormal_star_list:
                        ra = planetInfo[RA+1]
                        dec = planetInfo[DEC+1]
                    ra = float(ra)
                    ra = 360.0 - ra  # the data given by exoplanet.eu have reversed right ascension
                    dec = float(dec)
                    if distance != '':
                        distance = float(distance)
                    assert ra < 360 and ra > 0
                    assert dec > -90 and dec < 90
                    if distance != '' and distance < 10.0:
                        outputDict[method][0].append(ra)
                        outputDict[method][1].append(dec)
                        leng += 1
                except:
                    pass
#                    print planetInfo[RA]
#                    print planetInfo[DEC]
#                    print planetInfo[DISTANCE]
    print "<10pc " + str(leng)
    return outputDict


def sumLength(dict):
    leng = 0
    for key in dict:
        leng += len(dict[key][0])
    return leng
##############################################################################
##############################################################################
##############################################################################
        
        
def test():
    #testDict = {"astrometry": [[19,200,312],[5,-6,76]], "microlensing": [
    #        [122,311,248],[-82,70,33]]}
    #exoPlot(testDict, "testTitle")
    inputDict = readIn()
    fullDict = getRaDec(inputDict)
    earthLikeDict = getEarthLike(inputDict)
    within10pcDict = getWithin10pc(inputDict)
    #print within10pcDict
    exoPlot(fullDict,"All Stars with Confirmed Exoplanets")
    exoPlot(earthLikeDict,"All Confirmed Earth-Size Exoplanets")
    exoPlot(within10pcDict,"All Confirmed Exoplanets within 10pc")
#    for li in inputDict['TTV']:
#        print li
#        print
#        print li[65:75]
#        print
#        print
    #print inputDict['TTV']
    #print outputDict

test()