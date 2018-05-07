#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 15:21:24 2017

@author: zifan
"""

import math
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
STAR_T = 92  # effective temperature of star
STAR_R = 85  # radius of star (unit in solar radius)
STAR_M = 82  # mass of star (unit in solar mass)
PLANET_A = 14  # planet's orbit semi-major axis in AU
ORBIT_P = 11  # planet's orbital period in days
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
########################################################################
###### READ IN FUNCTIONS ###############################################
########################################################################
def readIn():
    inputDict = {}
    for method in DETECTION_METHODS:
        with open ("DetectionMethod/" + method + ".csv", "r") as infile:
            inputDict[method] = infile.readlines()
            for i in range(len(inputDict[method])):
                inputDict[method][i] = inputDict[method][i].split(",")
    return inputDict


def readInHabitable():
    input_list_RV = []
    input_list_transit = []
    with open("earth-like habitable RV.csv") as infile_RV:
        input_list_RV = infile_RV.readlines()
        #print len(input_list_RV)
        for i in range(len(input_list_RV)):
            input_list_RV[i] = input_list_RV[i].split(',')
        #print len(input_list_RV)
    with open("earth-like habitable transit.csv") as infile_transit:
        input_list_transit = infile_transit.readlines()
        #print len(input_list_transit)
        for i in range(len(input_list_transit)):
            input_list_transit[i] = input_list_transit[i].split(',')
        #print len(input_list_transit)
    return input_list_RV, input_list_transit


def findIndex():
    inputDict = readIn()
    titleList = inputDict["pulsar"][0]
    print titleList.index("ra")
    print titleList.index("dec")
    print titleList.index("mass")
    print titleList.index("radius")
    print titleList.index("star_distance")
    print titleList


########################################################################
###### PLOT FUNCTIONS ##################################################
########################################################################
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
    plt.figure(figsize=(65, 30))  # figure size to be determined
    plt.title(title, fontsize = 50)
    plt.xlabel("Right ascension", fontsize = 50)
    plt.ylabel("Declination", fontsize = 50)
    plt.axis([0,360,-90,90])
    plt.xticks([0,30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360], fontsize = 40)
    plt.yticks([-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90], fontsize = 40)
    ### define the color patches indicating detection methods
    #purple_patch = mpatches.Patch(color=PURPLE, label='Astrometry')
    #lime_patch = mpatches.Patch(color=LIME, label='Microlensing')
    #olive_patch = mpatches.Patch(color=OLIVE, label='Pulsar')
    red_patch = mpatches.Patch(color=RED, label='Transit')
    #yellow_patch = mpatches.Patch(color=YELLOW, label='Imaging')
    #aqua_patch = mpatches.Patch(color=AQUA, label='Other')
    blue_patch = mpatches.Patch(color=BLUE, label='Radial Velocity')
    #maroon_patch = mpatches.Patch(color=MAROON, label='TTV')
    #plt.legend(handles=[purple_patch,lime_patch,olive_patch,red_patch,
    #                       yellow_patch,aqua_patch,blue_patch,maroon_patch],
    #           loc='center left', bbox_to_anchor=(1, 0.5),fontsize=40)
    plt.legend(handles=[red_patch,blue_patch],
               loc='center left', bbox_to_anchor=(1, 0.5),fontsize=40)
    #leng = 0
    for method in dataDict:
        clr = color_dict[method]
        coordinateList  = dataDict[method]
        #leng += len(coordinateList[0])
        plt.plot(coordinateList[0], coordinateList[1], 'o',
                 markersize = 30, color = clr)
        count += 1
    filename = "plot/" + title + '.eps'
    plt.savefig(filename)
    #print leng


def effectiveFluxPlot(planet_list, bline_list, title):
    """ Parameters:
            1) planet_list is a list of length-4 lists. In those lists,
            the first element is name of the planet. The second element is the
            effective flux incident on the planet. The third element is Teff of
            host star. The fourth element is mass (if RV) or radius (if transit)
            of the planet in Earth mass or Earth radius.
            2) title is the title of the output plot.
    """
    bline_y = bline_list[0]
    inline_x = bline_list[1]
    outline_x = bline_list[2]
    
    plt.figure(figsize=(65, 30))  # figure size to be determined
    plt.title(title, fontsize = 50)
    plt.xlabel("Effective Flux Incident on the Planet", fontsize = 50)
    plt.ylabel("Stellar Effective Temperature (K)", fontsize = 50)
    plt.axis([2.5,0,2400,6200])
    plt.xticks([2, 1.5, 1, 0.5, 0], fontsize = 40)
    plt.yticks([2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000], fontsize = 40)
    plt.plot(inline_x, bline_y, 'r', linewidth=7)
    plt.plot(outline_x, bline_y, 'b', linewidth=7)
    for planet in planet_list:
        #print planet[1], "   ", planet[2]
        #print planet[1], " its mass/radius is: ", planet[3]
        plt.plot(planet[1], planet[2], 'o', markersize = planet[3]*15, color='#808080')
        plt.annotate(planet[0], (planet[1]+0.08, planet[2]-120), fontsize = 40)
    filename = "plot/" + title + "_run2" + ".eps"
    plt.savefig(filename)


########################################################################
###### PLANET FILTER FUNCTIONS #########################################
########################################################################
def getRaDec(inputDict):
    """ Reads in a dictionary in the style of output dictionary of readIn,
    and extracts the right ascension and declination information from it. """
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
                except:
                    pass
#                    print planetInfo[RA]
    return outputDict


def getEarthLike(inputDict):
    """ Returns a dictionary containing the information of all earth-like exoplanets.
    The definition of earth-like planet in this case is:
        (1) R < 2Re, and
        (2) M < 10Me
    """
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
                if radius != '' and mass != '' and radius < 2.0 and mass < 10.0 \
                and ra != '' and dec != '':
                    outputDict[method][0].append(ra)
                    outputDict[method][1].append(dec)
    return outputDict


def getWithin10pc(inputDict):
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
                except:
                    pass
#                    print planetInfo[RA]
#                    print planetInfo[DEC]
#                    print planetInfo[DISTANCE]
    return outputDict


def getHabitableZonePlanets(inputDict):
    """ Reads in a dictionary in the style of output dictionary of readIn,
    and extracts the right ascension and declination information from it. """
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
                planet_name = planetInfo[0]
                ra = planetInfo[RA]
                dec = planetInfo[DEC]
                star_radius = planetInfo[STAR_R]
                star_Teff = planetInfo[STAR_T]
                planet_a = planetInfo[PLANET_A]
                mass = planetInfo[MASS]
                radius = planetInfo[RADIUS]
                star_mass = planetInfo[STAR_M]
                orbital_period = planetInfo[ORBIT_P]
#                print method
#                print planetInfo[0]
#                print "Rstar: " + star_radius
#                print "Teff: " + star_Teff
                try:
                    ra = float(ra)
                    ra = 360.0 - ra  # the data given by exoplanet.eu have reversed right ascension
                    dec = float(dec)
                except:
                    ra = -1
                    dec = 0
                assert (ra < 360 and ra > 0) or ra == -1
                assert dec > -90 and dec < 90
                if inHabitableZone(star_radius, star_Teff, planet_a, star_mass, orbital_period) and ra != -1:
                    if isEarthLike(mass, radius):
                        print planet_name
                        outputDict[method][0].append(ra)
                        outputDict[method][1].append(dec)
    return outputDict


########################################################################
###### CALCULATION FUNCTIONS ###########################################
########################################################################
def inHabitableZone(R, T, a, M, P):
    """ return a boolean, true if the planet is in the habitable zone of its
        host star, false otherwise.
        r - star radius
        T - star effective temperature
        a - planet orbit semi-major axis
        M - star mass
        P - orbital period"""
#    print "R: " + R
#    print "T: " + T
#    print "a: " + a
#    print
    if a=='' and P!='' and M!='':
        P = float(P)
        M = float(M)
        a = getSemiMajorAxis(P, M)
    try:
        R = float(R)
        T = float(T)
        a = float(a)
    except:
        return False
    luminosity = calculateLuminosity(T, R)
    #print luminosity
    # if luminosity < 0.2 or luminosity > 1.7:
    #     return False
    Din, Dout = calculateDistance(T, R)
#    print "distance in: " + str(Din)
#    print "distance out: " + str(Dout)
#    print "semi major axis: " + str(a)
    if a >= Din and a <= Dout:
        #print name
        return True
    return False


def isEarthLike(m, r):
    """ return true if the planet's mass < 10 earth mass and radius < 2 earth
        radius. """
    try:
        m = float(m)
        m = m * MASS_CORRECTION_FACTOR
    except:
        pass
    try:
        r = float(r)
        r = r * RADIUS_CORRECTION_FACTOR
    except:
        pass
    if r<2.0 or m<10.0:
        #print name
        return True
    else:
        return False


def getSemiMajorAxis(P, M):
    """ calculate the planet's semi-major axis using Kepler's third law.
        P is orbital period in days, M is star mass in solar mass. """
    SOLAR_MASS = 1.98855E30  # solar mass in kg
    AU_CORRECTION_FACTOR = 1495978.70700  # astronomical unit in km
    EARTH_YEAR = 365.242189  # 1 sidereal year in days
    P = P/EARTH_YEAR  # convert the unit of P to earth year, so that the result will be in AU
    #SECONDS_IN_YEAR = 60 * 60 * 24 * 365.242189
    #P = P * SECONDS_IN_YEAR  # convert orbital period into seconds
    M = M * SOLAR_MASS  # convert star mass to unit kg
    G = 6.67408E-11  # gravitational constant
    a = ((G*M*(P**2.0))/(4.0*(math.pi**2.0)))**(1.0/3.0)
    return a/AU_CORRECTION_FACTOR


def calculateDistance(Teff, R):
    """ return the inner and outer edge distance of the star's habitable zone
        (unit in AU), given effective temperature of the star. """
    if Teff < 2000 or Teff > 12000:
        return 0, 0
    SUN_IN = 1.7665
    A_IN = 1.3351E-4
    B_IN = 3.1515E-9
    C_IN = -3.3488E-12
    SUN_OUT = 0.324
    A_OUT = 5.3221E-5
    B_OUT = 1.4288E-9
    C_OUT = -1.1049E-12
    
    T = Teff - 5780
    Seff_in = SUN_IN + A_IN*T + B_IN*(T**2) + C_IN*(T**3)
    Seff_out = SUN_OUT + A_OUT*T + B_OUT*(T**2) + C_OUT*(T**3)
    luminosity = calculateLuminosity(Teff, R)
    distance_in = math.sqrt(luminosity/Seff_in)
    distance_out = math.sqrt(luminosity/Seff_out)
    return distance_in, distance_out


def calculateEffectiveFluxBoundary(Teff):
    """ Use to get the two boundary lines of HZ. """
    SUN_IN = 1.7665
    A_IN = 1.3351E-4
    B_IN = 3.1515E-9
    C_IN = -3.3488E-12
    SUN_OUT = 0.324
    A_OUT = 5.3221E-5
    B_OUT = 1.4288E-9
    C_OUT = -1.1049E-12
    
    T = Teff - 5780
    Seff_in = SUN_IN + A_IN*T + B_IN*(T**2) + C_IN*(T**3)
    Seff_out = SUN_OUT + A_OUT*T + B_OUT*(T**2) + C_OUT*(T**3)
    return Seff_in, Seff_out


def calculateEffectiveFluxOnPlanet(inlist):
    """ Return the effective flux on planet. Inlist is a line of exoplanet.eu
        data of a certain planet. """
    Teff = float(inlist[STAR_T])
    R = float(inlist[STAR_R])
    try:
        a = float(inlist[PLANET_A])
    except:
        a = getSemiMajorAxis(float(inlist[ORBIT_P]), float(inlist[STAR_M]))
    #print inlist[0]
    #print Teff
    #print R
    #print a
    
    lum = calculateLuminosity(Teff, R)
    Seff = lum / a**2
    return Seff


def calculateLuminosity(Teff, R):
    """ return luminosity of star, unit in solar luminosity. """
    R_SUN = 695.7E6  # radius of sun in meters
    SB_CONSTANT = 5.670367E-8  # the Stefan–Boltzmann constant
    L_SUN = 3.828E26  # solar luminosity
    R = R_SUN * R  # unit conversion from solar radius to meter
    area = 4.0*math.pi*(R**2.0)
    lumniosity = SB_CONSTANT * area * (Teff**4.0)
    #print lumniosity
    lumniosity = lumniosity/L_SUN
    #print lumniosity
    return lumniosity
    

def sumNumber(di):
    """ return the number of elements in a dictionary where all elements are
        lists. """
    count = 0
    for key in di:
        count += len(di[key][0])
    return count
##############################################################################
##############################################################################
##############################################################################

def testCalculateDistance():
    print calculateDistance(5780, 1)
    print calculateDistance(3000, 0.6)
    print calculateDistance(3050, 0.141)
    print calculateDistance(2550, 0.117)
    print calculateDistance(5710, 0.17)
    print calculateDistance(4723.0, 0.706)
    print calculateDistance(4249.0, 0.56)


def testGetSemiMajorAxis():
    print getSemiMajorAxis(365*0.24, 1)
    print getSemiMajorAxis(365*0.62, 1)
    print getSemiMajorAxis(365, 1)
    print getSemiMajorAxis(365*1.88, 1)
    print getSemiMajorAxis(4332, 1)
    print getSemiMajorAxis(365*11.86, 1)
    print getSemiMajorAxis(365*29.46, 1)
    print getSemiMajorAxis(365*84.01, 1)
    print getSemiMajorAxis(4.1876244, 1.272)  # expected 0.05511
    print getSemiMajorAxis(34.94, 0.912)  # expected 0.2055


def testInHabitableZone():
    print inHabitableZone(1, 5780, 1, 1, 1)
    print inHabitableZone(1, 5780, 2, 1, 1)
    print inHabitableZone(1, 5780, 0.5, 1, 1)
    print inHabitableZone(1, 5780, 0.8, 1, 1)
    print inHabitableZone(1, 5780, 1.5, 1, 1)
    print inHabitableZone(0.141, 3050, 0.04, 1, 1)


def testIndex():
    inputDict = readIn()
    print inputDict["transit"][0][STAR_M]
    print inputDict["transit"][0][ORBIT_P]
    
        
def test():
    inputDict = readIn()
    habitableDict = getHabitableZonePlanets(inputDict)
    print sumNumber(habitableDict)
    #print habitableDict
    exoPlot(habitableDict, "Earth-sized Planets in Habitable Zone")


def runPlotHabitble():
    RV_list, transit_list = readInHabitable()
    RV_outlist = []
    transit_outlist = []
    #print RV_list
    #print transit_list
    for planet in RV_list:
        #print calculateEffectiveFluxOnPlanet(planet)
        RV_outlist.append([planet[0],
                           calculateEffectiveFluxOnPlanet(planet),
                           float(planet[STAR_T]),
                           float(planet[MASS])*MASS_CORRECTION_FACTOR])
    for planet in transit_list:
        #print calculateEffectiveFluxOnPlanet(planet)
        transit_outlist.append([planet[0][1:],
                                calculateEffectiveFluxOnPlanet(planet),
                                float(planet[STAR_T]),
                                float(planet[RADIUS])*RADIUS_CORRECTION_FACTOR])
    #print RV_outlist
    #print transit_outlist
    
    boundary_line = [[], [], []]
    for i in range(2400, 6500):
        Sin, Sout = calculateEffectiveFluxBoundary(i)
        boundary_line[0].append(i)
        boundary_line[1].append(Sin)
        boundary_line[2].append(Sout)
    effectiveFluxPlot(RV_outlist, boundary_line, "Earth-sized Planets in Habitable Zone (RV)")
    effectiveFluxPlot(transit_outlist, boundary_line, "Earth-sized Planets in Habitable Zone (transit)")


def run():
    inputDict = readIn()
    fullDict = getRaDec(inputDict)
    earthLikeDict = getEarthLike(inputDict)
    within10pcDict = getWithin10pc(inputDict)
    exoPlot(fullDict,"All Stars with Confirmed Exoplanets")
    exoPlot(earthLikeDict,"All Confirmed Earth-Size Exoplanets")
    exoPlot(within10pcDict,"All Confirmed Exoplanets within 10pc") 


#test()
#testCalculateDistance()
#testGetSemiMajorAxis()
#testIndex()
#testInHabitableZone()
runPlotHabitble()

#run()
