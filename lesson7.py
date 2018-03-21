# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 15:30:31 2017

@author: bwinters
"""

def get_planet_data(filename):
    import csv
    # open the file where the planet's data is stored
    # and create a reader
    mFile=open(filename,'rt')
    reader=csv.reader(mFile)
    # read the line in the file
    for line in reader:
        name,e,period,a=line
    # convert the numerical data from string to floating point
    e=float(e)
    period=float(period)
    a=float(a)
    # return the data
    return (name,e,period,a)
    
def planet_orb(planetfile):
    import numpy as np
    # period=the number of days it takes for the planet
    # to go around the sun.
    # filename is where the output will be stored
    name,e,period,a=get_planet_data(planetfile)
    outputfile=name+'_orbit.txt'
    n=2*np.pi/period
    orbit=[]
    T=int(np.round(period,0)+1)
    for i in range(T):
        M=i*n
        nu=get_nu(M,e)
        r=a*(1-e*e)/(1+e*np.cos(nu))
        x=r*np.cos(nu)
        y=r*np.sin(nu)
        orbit.append((i,x,y))
    save_orb(orbit,outputfile)
    draw_orb(outputfile)
    return 'Done'

def get_nu(M,e):
    import numpy as np
    E=M
    tol=1e-14
    testFunction=E-e*np.sin(E)-M
    while abs(testFunction)>tol:
        der=1-e*np.cos(E)
        E=E-testFunction/der
        testFunction=E-e*np.sin(E)-M
    tnu2=np.sqrt((1+e)/(1-e))*np.tan(E/2)
    nu=2*np.arctan(tnu2)
    if nu<0:
        nu=nu+2*np.pi
    return nu 
    
def save_orb(plist,filename):
    import csv
    mFile=open(filename,'wt')
    writer=csv.writer(mFile)
    for row in plist:
        writer.writerow(row)
    mFile.close()
    return 'Orbit Saved'
    
def draw_orb(filename):
    import csv
    import pylab as plt
    X=[]
    Y=[]
    mFile=open(filename,'rt')
    reader=csv.reader(mFile)
    for row in reader:
        if row!=[]:
            i,x,y=row
            X.append(float(x))
            Y.append(float(y))
    mFile.close()
    plt.axis('equal')
    plt.plot(X,Y)
    name,ext=filename.split('.')
    pictfile=name+'.jpg'
    plt.savefig(pictfile,format='jpg')
    return 'Orbit Drawn'