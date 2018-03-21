# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 10:44:09 2017

@author: bwinters
"""

def circle(radius):
    import numpy as np
    interval=np.linspace(0,2*np.pi,1000)
    C=[]
    for t in interval:
        x=radius*np.cos(t)
        y=radius*np.sin(t)
        C.append([x,y])
    return C
        
    

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
