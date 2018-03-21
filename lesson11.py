# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:42:38 2017

@author: bwinters
"""
def celestial_coords(planet1,planet2):
    import csv
    import numpy as np
    mFrame=read_relative_position(planet1,planet2)
    N=len(mFrame.values)
    savefile=planet2+'in'+planet1+'skyTable.txt'
    sfile=open(savefile,'wt')
    writer=csv.writer(sfile)    
    for i in range(N):
        day,x,y,z=mFrame.values[i]
        rho=np.sqrt(x*x+y*y+z*z)
        x,y,z=x/rho,y/rho,z/rho
        r=np.sqrt(x*x+y*y)
        delta=np.arcsin(z)
        delta=rad_to_dms(delta)
        alpha=2*np.arctan((1-x/r)/(y/r))
        if alpha<0:
            alpha=alpha+2*np.pi
        
        alpha=rad_to_hms(alpha)
        #delta=180*delta/np.pi
        date=standard_date(day)
        writer.writerow([date, alpha, delta])

    
    return 'Done'
    

def read_relative_position(planet1,planet2):
    import pandas
    import csv
    datafile='coords_'+planet2+'_relativeto_'+planet1+'.txt'
    file=open(datafile,'rt')
    reader=csv.reader(file)
    x=[]
    y=[]
    z=[]
    day=[]
    for line in reader:
        if line!=[]:
            d0,x0,y0,z0=line
            x.append(float(x0))
            y.append(float(y0))
            z.append(float(z0))
            day.append(int(d0))
    data={
    'day':day,
    'x':x,
    'y':y,
    'z':z}
    frame=pandas.DataFrame(data)
    return frame 
def rad_to_hms(rad):
    import numpy as np
    hours=12*rad/np.pi
    h=int(hours-hours%1)
    hours=(hours-h)*60
    m=int(hours-hours%1)
    hours=(hours-m)*60
    s=int(hours-hours%1)
    out=str(h)+'h'+str(m)+'m'+str(s)+'s'
    return out

def rad_to_dms(rad):
    import numpy as np
    sign=np.sign(rad)
    rad=np.abs(rad)
    degrees=180*rad/np.pi
    d=int(degrees-degrees%1)
    degrees=(degrees-d)*60
    m=int(degrees-degrees%1)
    degrees=(degrees-m)*60
    s=int(degrees-degrees%1)
    out=str(d)+'d'+str(m)+'m'+str(s)+'s'
    if sign>0:
        out=out+' North'
    if sign<0:
        out=out+' South'
    return out    
    
    
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
    d=int(order)
    mdy=str(m)+'/'+str(d)+'/'+str(year)
    return(mdy)
     