# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 15:26:39 2019

@author: Bennett

This program is made to plot and analyze the data gathered from flowing argon through a small pressure cell and measuring the phase shift
of an interferometer. 
METHODS:
This code is made for the following data collection process:
    1) Collect an initial data by rolling through a full 2pi phase shift at an initial guage pressure.
    2) Collect data at some initial pressure value..
    3) Collect data at a second pressure value.
    4) Collect a final 2pi phase shift at the second pressure.
FILE FORMATTING:
To ensure the program works as it is currently written, there should be three file types from data collection:
    1) 'initialfit'= name for the 2pi ellipse data taken at the initial pressure value.  
    2) 'C10.00'= name for the data taken at the first pressure value. Should be the input channel (C1/C3 for double photodiode setup) and a numerical value. 
    3) 'finalfit'=name for the 2pi ellipse data taken at the final pressure value. 
NOTE: THIS CODE WILL NOT WORK FOR A PHASE SHIFT BEYOND 2PI RADS!!! The arctan relationship cannot account for more than one full oscillation. 
"""
# =============================================================================
# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from matplotlib.patches import Ellipse
from copy import copy
import ellipsefit as ef #import ellipsefit code
# =============================================================================
# Variables for plotting #
#create a number of arrays. Two x,y arrays for the experimental data. Two x,y arrays for full elipse fitting
pressures=[] #Will end up as a list of the guage pressures measured
xlist1=[] #x values of first pressure
ylist1=[] #y values of first pressure
xlist2=[] #x values of second pressure
ylist2=[] #y values of second pressure
xfiti=[] #first pressure x fit data
yfiti=[] #first pressure y fit data
xfitf=[] #second pressure x fit data
yfitf=[] #second pressure y fit data
# =============================================================================
# Path names for data #
d1=r'C:\Users\Bennett\Desktop\datascan4' #change path as needed
# =============================================================================
# Messing with the data a bit to format into arrays #
#put words "fit" and "initial"/"final" in file name for data for fitting. 
#Code assumes one of the files will be C10.00 or C30.00
for filename in os.listdir(d1):
    if 'C1' in str(filename) and 'fit' not in str(filename):
        name=str(filename)
        pressures.append(float(name[2:-4]))
        full_path = os.path.join(d1, filename)
        data=pd.read_csv(full_path,sep=" ",header=5,names=['v'])
        if '0.00' in str(filename):
            for item in data.v:
                xlist1.append(item)
        else:
            for item in data.v:
                xlist2.append(item)
    if 'C3' in str(filename) and 'fit' not in str(filename):
        full_path = os.path.join(d1, filename)
        data=pd.read_csv(full_path,sep=" ",header=5,names=['v'])
        if '0.00' in str(filename):
            for item in data.v:
                ylist1.append(item)
        else:
            for item in data.v:
                ylist2.append(item)
    if 'C1' in str(filename) and 'initial' in str(filename):
        full_path = os.path.join(d1, filename)
        data=pd.read_csv(full_path,sep=" ",header=5,names=['v'])
        for item in data.v:
            xfiti.append(item)
    if 'C3' in str(filename) and 'initial' in str(filename):
        full_path = os.path.join(d1, filename)
        data=pd.read_csv(full_path,sep=" ",header=5,names=['v'])
        for item in data.v:
            yfiti.append(item)
    if 'C1' in str(filename) and 'final' in str(filename):
        full_path = os.path.join(d1, filename)
        data=pd.read_csv(full_path,sep=" ",header=5,names=['v'])
        for item in data.v:
            xfitf.append(item)
    if 'C3' in str(filename) and 'final' in str(filename):
        full_path = os.path.join(d1, filename)
        data=pd.read_csv(full_path,sep=" ",header=5,names=['v'])
        for item in data.v:
            yfitf.append(item)
# =============================================================================
#Curve Fitting#
# uses algorithm from paper Fitzgibbon 1996. http://cseweb.ucsd.edu/~mdailey/Face-Coord/ellipse-specific-fitting.pdf
dataf=np.column_stack((xfitf,yfitf))
datai=np.column_stack((xfiti,yfiti))
#outputs parameters and Ellipse patches. See ellipsefit.py code.
paramsf,ellf,ell1=ef.fitEllipse(dataf,'2.00psig Elliptical Fit')
paramsi,elli,ell2=ef.fitEllipse(datai,'0.00psig Elliptical Fit')
#create two copy variables for plotting
new1=copy(ell1)
new2=copy(ell2)
tabdata=[np.round(paramsi,3),np.round(paramsf,3)]
xsf=[]
ysf=[]
xsi=[]
ysi=[]
for item in ellf:
    xsf.append(item[0])
    ysf.append(item[1])
for item in elli:
    xsi.append(item[0])
    ysi.append(item[1])
# =============================================================================    
#Converting to x,y points to angles
#find max x,y values for fit ellipse to account for angular shift/non-circular geometry
xmax=np.sqrt((paramsi[2]**2)*(np.cos(paramsi[4])**2)+(paramsi[3]**2)*(np.sin(paramsi[4])**2))
ymax=np.sqrt((paramsi[2]**2)*(np.sin(paramsi[4])**2)+(paramsi[3]**2)*(np.cos(paramsi[4])**2))
def getThetai(values,xmax,ymax):
    thetas=[]
    for item in values:
        thetas.append(np.arctan(((item[0]-paramsi[0])/(item[1]-paramsi[1]))*(ymax/xmax)))
    meantheta=np.mean(thetas)
    stdtheta=np.std(thetas)
    return meantheta,stdtheta
def getThetaf(values,xmax,ymax):
    thetas=[]
    for item in values:
        thetas.append(np.arctan(((item[0]-paramsf[0])/(item[1]-paramsf[1]))*(ymax/xmax)))
    meantheta=np.mean(thetas)
    stdtheta=np.std(thetas)
    return meantheta,stdtheta
#zip x,y points for use in getTheta functions
points1=zip(xlist1,ylist1)
points2=zip(xlist2,ylist2)
meantheta1i, stdtheta1i=getThetai(points1,xmax,ymax)
meantheta2i, stdtheta2i=getThetai(points2,xmax,ymax)
points1=zip(xlist1,ylist1)
points2=zip(xlist2,ylist2)
meantheta1f, stdtheta1f=getThetaf(points1,xmax,ymax)
meantheta2f, stdtheta2f=getThetaf(points2,xmax,ymax)
phaseshifti=np.absolute(meantheta1i-meantheta2i)
stdtoti=stdtheta1i+stdtheta2i
phaseshiftf=np.absolute(meantheta1f-meantheta2f)
stdtotf=stdtheta1f+stdtheta2f
#quick print statement to give general idea of data collected. Checks to ensure that the outputs using both initial and final fits line up. 
print("Phaseshift initial: ", str(np.degrees(phaseshifti))+" pm "+str(np.degrees(stdtoti)), "Phaseshift final: ",  str(np.degrees(phaseshiftf))+" pm "+str(np.degrees(stdtotf)))
deltapi=(14.6959/.000281)*((635*10**(-9)*phaseshifti)/(2*np.pi*2*10**(-3))+1.000293-1)
deltapf=(14.6959/.000281)*((635*10**(-9)*phaseshiftf)/(2*np.pi*2*10**(-3))+1.000293-1)
deltapuncert=(14.6959/.000281)*((635*10**(-9)*stdtoti)/(2*np.pi*2*10**(-3))+1.000293-1)
print("P initial: ", deltapi,"P final: " , deltapf)
# =============================================================================    
# Plotting #
plt.figure(1) #Figure displays elliptical fitting data and calculated function
plt.rcParams.update({'font.size': 36})
plt.scatter(xfiti,yfiti, c='b', alpha=.1, label='0.00psig Fit Data')
plt.scatter(xfitf,yfitf, c='r', alpha=.1, label='2.00psig Fit Data')
plt.gca().add_patch(ell1)
plt.gca().add_patch(ell2)
plt.xlabel("Photodiode1 Reading (V)")
plt.ylabel("Photodiode2 Reading (V)")
plt.legend(loc='upper right')

plt.figure(2) #Table of parameters calculated from algorithm
plt.title('Ellipse Fitting Parameters')
plt.axis('off')
columns=('Major axis','Minor axis','cx','cy','Rotation (degrees)')
rows=(str(pressures[0])+'psig',str(pressures[1])+'psig')
table=plt.table(cellText=tabdata,rowLabels=rows,colLabels=columns, loc='center')
table.auto_set_font_size(False)
table.set_fontsize(36)
table.scale(1,2)

plt.figure(3) #experimental data points and calculated values
plt.rcParams.update({'font.size': 36})
plt.gca().add_patch(new1)
plt.gca().add_patch(new2)
plt.scatter(xlist1,ylist1, alpha=.1, label=str(pressures[0])+ 'psig backing pressure')
plt.scatter(xlist2,ylist2, alpha=.1, label=str(pressures[1])+r'$\pm$'+'.01 psig backing pressure')
plt.xlabel("Photodiode1 Reading (V)")
plt.ylabel("Photodiode2 Reading (V)")
plt.text(.67,5.1,'Relative '+r'$\Delta\phi$'+'(rads): '+ str(round(phaseshifti,3)) + r'$\pm$' + str(round(stdtoti,3)))
plt.text(.67,4.9,'P(psi): '+str(round(deltapi,3))+ r'$\pm$' + str(round(deltapuncert,3)))
plt.legend()

plt.show