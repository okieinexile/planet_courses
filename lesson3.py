# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 13:55:38 2017

@author: bwinters
"""

def ellipse(a,b,c):
    import numpy as np
    interval=np.linspace(0,2*np.pi,100)
    points=[]
    for t in interval:
        x=a*np.cos(t)-c
        y=b*np.sin(t)
        points.append((x,y))
        
    return points 
    
def write_list(pList,filename):
    import csv
    out='Error'
    mFile=open(filename,'wt')
    writer=csv.writer(mFile)
    for row in pList:
        writer.writerow(row)
        out='Done'
    mFile.close()
           
    return out
    
def read_list(filename)  :
    import csv
    mFile=open(filename,'rt')
    S=[]
    T=[]
    reader=csv.reader(mFile)
    for row in reader:
        if len(row)==2:
            X,Y=row
            x=float(X)
            y=float(Y)
            S.append(x)
            T.append(y)
            
    mFile.close()
    return [S,T]
    
def plot_list(pair):
    import pylab as plt
    x,y=pair
    plt.axis('equal')
    plt.plot(x,y)
    return 'Done'