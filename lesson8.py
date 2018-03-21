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
        name,e,period,a,perihelion=line
    # convert the numerical data from string to floating point
    e=float(e)
    period=float(period)
    a=float(a)
    date=ord_day(perihelion)
    # return the data
    return (name,e,period,a,date)
    
def planet_orb(planetfile):
    import numpy as np
    # period=the number of days it takes for the planet
    # to go around the sun.
    # filename is where the output will be stored
    name,e,period,a,perihelion=get_planet_data(planetfile)
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
        orbit.append((i+perihelion,x,y))
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
    
#####
####
    
def days_after(mdy1,mdy2):
    # enter dates mdy1 and mdy2 in the form 'mm/dd/yyyy'
    # mdy1 is assumed to be earlier than mdy2
    d1=ord_day(mdy1)
    d2=ord_day(mdy2)
    dif=d2-d1
    return dif

def ord_day(mdy):
    # mdy is a string of the form 'mm/dd/yyyy' of the day in question
    # basemdy is a string of the form 'mm/dd/yyyy' of a base day in the past
    #   which is after the change to Gregorian Time. We set it as '1/1/1800'
    basemdy='1/1/1800'
    ds=0         
    bm,bd,by=basemdy.split('/')
    bm=int(bm)
    bd=int(bd)
    by=int(by)
    m,d,y=mdy.split('/')
    m=int(m)
    d=int(d)
    y=int(y)
    if y< by:
        print('Error: Date too early')
        return None
    for i in range(by,y):
        if is_leap(i):
            ds=ds+366
        else:
            ds=ds+365
    for i in range(1,m):
        ds=ds+month_length(i,y)
    ds=ds+d
    
    return ds
    
def is_leap(year):
    if ((year%4)!=0):
        return False
    if (year%100)!=0:
        return True
    if (year%400)!=0:
        return False
    return True
    
def month_length(m,y):
    months=[31,28,31,30,31,30,31,31,30,31,30,31] 
    L=months[m-1]
    if is_leap(y) and m==2:
        L=29
    return L

def standard_date(order):
    year=1800
    yr_ln=365
    while order>yr_ln:
        order=order-yr_ln
        year=year+1
        if is_leap(year):
            yr_ln=366
        else:
            yr_ln=365
    m=1
    mn_ln=month_length(m,year)
    while order>mn_ln:
        order=order-mn_ln
        m=m+1
        mn_ln=month_length(m,year)
    d=order
    mdy=str(m)+'/'+str(d)+'/'+str(year)
    return(mdy)
     