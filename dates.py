# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 15:02:36 2017

@author: bwinters
"""
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
    