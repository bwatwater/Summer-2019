# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:26:25 2019

@author: Bennett
Based off of Fitzgibbon, Pilu, Fisher (1996): http://cseweb.ucsd.edu/~mdailey/Face-Coord/ellipse-specific-fitting.pdf
Produces an elliptical fit of an array of data points. 
"""

import numpy
import pandas as pd
import matplotlib
from matplotlib.patches import Ellipse
#cont=paired x,y data in an array: [[x1,y1],[x2,y2]...[xn,yn]], name=string name associated with the elliptical fit. 
def fitEllipse(cont,name):

    x=cont[:,0]
    y=cont[:,1]

    x=x[:,None]
    y=y[:,None]
#create an array, D, of elliptical polynomial values from points
#a1x**2+a2xy+a3y**2+a4x+a5y+a6=0
    D=numpy.hstack([x*x,x*y,y*y,x,y,numpy.ones(x.shape)]) 
#perform a series of linear algebra steps to find lagrangian multipliers.
    S=numpy.dot(D.T,D)
    C=numpy.zeros([6,6])
    C[0,2]=C[2,0]=2
    C[1,1]=-1
    E,V=numpy.linalg.eig(numpy.dot(numpy.linalg.inv(S),C))
    n=numpy.argmax(E)
#Output various polynomial parameters for best fit stored in array, 'a'.
    a=V[:,n]

    #-------------------Fit ellipse-------------------
    b,c,d,f,g,a=a[1]/2., a[2], a[3]/2., a[4]/2., a[5], a[0]
#Use a number of geometric identities to convert the polynomial coefficients to major and minor axis, center point, and tilt angle values.
    num=b*b-a*c
    cx=(c*d-b*f)/num
    cy=(a*f-b*d)/num
    angle=0.5*numpy.arctan(2*b/(a-c))*180/numpy.pi
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*numpy.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*numpy.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    a=numpy.sqrt(abs(up/down1))
    b=numpy.sqrt(abs(up/down2))

    #---------------------Get path---------------------
#returns a matplotlib patch associated with an ellipse.
    ell=Ellipse((cx,cy),a*2.,b*2.,angle, lw=4, fill=False, color=(numpy.random.randint(100)/100,numpy.random.randint(100)/100,numpy.random.randint(100)/100), label=name)
    ell_coord=ell.get_verts()
#params values: cx=center point x value, cy=center point y value, a=major axis length, b=minor axis length, angle=tilt angle
    params=[cx,cy,a,b,angle]
    return params,ell_coord,ell