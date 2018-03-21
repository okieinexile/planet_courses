# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 15:15:22 2017

@author: bwinters
"""

def cir_orb(period,filename):
    import numpy as np
    # period=the number of days it takes for the planet
    # to go around the sun.
    # filename is where the output will be stored
    M=2*np.pi/period
    orbit=[]
    for i in range(period):
        angle=i*M
        x=np.cos(angle)
        y=np.sin(angle)
        orbit.append((i,x,y))
    save_orb(orbit,filename)
    draw_orb(filename)
    return 'Done'
    
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
    