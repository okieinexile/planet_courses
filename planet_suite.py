# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 16:44:17 2017

@author: bwinters
"""
import numpy as np
import pandas as pd

def planet_table(planet):
    planet1='earth'
    planet2=planet
    planet_orb(planet1)
    planet_orb(planet2)
    put_in_position(planet1)
    put_in_position(planet2)
    relative_position(planet1,planet2)
    celestial_coords(planet)
    
    return 'Done'

#==============================================================================
# This functions will take a planet that is in ecliptic coordinates
# with respect to the earth and will transform them to
# equatorial coordinates with respect to the earth. They
# will take date from <planet>_relativeto_earth_ecliptic.txt
# and will save it at <planet>_table.txt.  An intermediate file
# <planet>_relativeto_earth_equatorial.txt will be
# created.
#==============================================================================
def celestial_coords(planet):
    import csv
    import numpy as np
    M=rotation(23.44,[1,0,0])
    mFrame=read_relative_position(planet)
    N=len(mFrame.values)
    savefile=planet+'_table.txt'
    intmdfile=planet+'_reltativeto_earth_equatorial.txt'
    sfile=open(savefile,'wt')
    writer=csv.writer(sfile)    
    ifile=open(intmdfile,'wt')
    writeri=csv.writer(ifile)
    for i in range(N):
        day,x,y,z=mFrame.values[i]
        w=np.matrix([x,y,z])
        v=M*w.T
        x,y,z=v[0,0],v[1,0],v[2,0]
        writeri.writerow([day,x,y,z])
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
    sfile.close()
    ifile.close()

    
    return 'Done'
    

def read_relative_position(planet):
    import pandas
    import csv
    datafile=planet+'_relativeto_earth_ecliptic.txt'
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
    


#==============================================================================
# These functions will put one planet into ecliptic coordinates with respect 
# to another.  The input will be files 
# <planet1>_ecliptic.txt  and <planet2>_ecliptic.txt 
# that contain the ecliptic coordinate 
# data of the respective planets.  These will be output to a file
# <planet2>_relativeto_<planet1>_ecliptic.txt. 
#==============================================================================

def relative_position(planet1, planet2):
    # planet1 is to be the base planet; in most cases this will be earth
    # planet2 will be the planet whose position we are plotting
    import csv
    frame1=read_planet(planet1)
    frame2=read_planet(planet2)
    days1=frame1.day.values
    days2=frame2.day.values
    begin=max([min(days1),min(days2)])
    end=min([max(days1),max(days2)])
    day=begin
    Flt1=((frame1.day>=begin)&(frame1.day<=end))
    F1=frame1[Flt1]
    Flt2=((frame2.day>=begin)&(frame2.day<=end))
    F2=frame2[Flt2]
    savefile=planet2+'_relativeto_'+planet1+'_ecliptic.txt'
    sfile=open(savefile,'wt')
    writer=csv.writer(sfile)    
    for i in range(end-begin):
        x1,y1,z1=F1.values[i][1],F1.values[i][2],F1.values[i][3]
        x2,y2,z2=F2.values[i][1],F2.values[i][2],F2.values[i][3]
        v1,v2,v3=x2-x1,y2-y1,z2-z1
        writer.writerow([day,v1,v2,v3])
        day=day+1
    sfile.close()
    return None

def read_planet(planet):
    file=open(planet+'_ecliptic.txt')
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
  

#==============================================================================
# These functions will transform planetary orbits from their planar form
# which are called from <planet>_orbit.txt to ecliptic coordinates.  
# They will save those coordinates to a file
# <planet>_ecliptic.txt They will use a file <planet>_data2.txt
#  that contains orbital elements of the given planet.
# These orbital elements will be put in a rotation matrix and saved in a file
# <planet>_matrix.txt.   
#==============================================================================
    
def put_in_position(planet):
    get_matrix(planet)
    mFmt='{}_matrix.csv'
    matfile=mFmt.format(planet)
    mat=pd.read_csv(matfile,usecols=[1,2,3])
    #print(mat)
    R=np.matrix(mat)
    o_fmt='{}_orbit.csv'
    orbfile=o_fmt.format(planet)
    o_df=pd.read_csv(orbfile)
    day=o_df['day']
    X=list(o_df['X'])
    Y=list(o_df['Y'])
    Z=len(X)*[0.0]
    V=zip(X,Y,Z)
    V=list(map(lambda v:np.matrix(v).transpose(), V))
    W=list(map(lambda v:R*v,V))
    Wx=list()
    Wy=list()
    Wz=list()
    for w in W:
        Wx.append(w[0,0])
        Wy.append(w[1,0])
        Wz.append(w[2,0])
    data={
        'day':day,
        'x':Wx,
        'y':Wy,
        'z':Wz
    }
    pos_frame=pd.DataFrame(data)
    p_fmt='{}_ecliptic.csv'
    posfile=p_fmt.format(planet)
    pos_frame.to_csv(posfile)
#
#    for line in reader:
#        if line!=[]:
#            day,x,y=line
#            x=float(x)
#            y=float(y)
#            v=np.matrix([x,y,0])
#            w=R*v.T
#            writer.writerow([day,w[0,0],w[1,0],w[2,0]])
#    ofile.close()
#    pfile.close()
    return None

def get_matrix(planet):
    fmt='{}_data2.txt'
    filename=fmt.format(planet)
    with open(filename,'rt') as myFile:
        name,Omega,i,omega=myFile.read().split(',')
        Omega=float(Omega)
        i=float(i)
        omega=float(omega)
    R1=rotation(Omega,[0,0,1])
    R2=rotation(i,[np.cos(Omega),np.sin(Omega),0])
    R3=rotation(omega,[R2[0,2],R2[1,2],R2[2,2]])
    R=R3*R2*R1
    mFmt='{}_matrix.csv'
    outfile=mFmt.format(planet)
    mat_frame=pd.DataFrame(R)
    mat_frame.to_csv(outfile)
    return 'Matrix Done'

def rotation(theta, axisList):
    import numpy as np
    axis=np.matrix(axisList)
    phi=theta*np.pi/180
    u=(1/(axis*axis.T)[0,0])*axis
    v=np.sin(phi/2)*u
    r=np.cos(phi/2)
    v1=v[0,0]
    v2=v[0,1]
    v3=v[0,2]
    p=r**2-(v*v.T)[0,0]
    M=np.matrix( ( (p+2*v1**2,  2*(v1*v2-r*v3), 2*(v1*v3+r*v2)),
                 (2*(v1*v2+r*v3), p-2*v2**2, 2*(v2*v3-r*v1)),
                 (2*(v1*v3-r*v2),  2*(v2*v3+r*v1), p+2*v3**2)))
    
    return(M)    
    
#==============================================================================
# Load planetary data from a file of filename <planet>_data1.txt. 
# Calculate the planar orbit.  
# Save that data to a file <planet>_orbit.txt.    
#==============================================================================
    
def get_planet_data(planet):
    """Given the name planet in lowercase, this will open the file where
    the data is stored and return the planet's
    name
    eccentricity
    period
    semi-major axis
    and date of perihelion
    """
    # open the file where the planet's data is stored
    # and create a reader
    fmt='{}_data1.txt'
    filename=fmt.format(planet)
    mFile=open(filename,'rt')
    name,e,period,a,perihelion=mFile.read().split(',')
    # convert the numerical data from string to floating point
    e=float(e)
    period=float(period)
    a=float(a)
    date=ord_day(perihelion)
    # return the data
    return (name,e,period,a,date)
    
def planet_orb(planet):
    """planet is the name of a planet in lowercase.  This will output to a file
    <planet>_orbit.txt. It uses get_planet_data and get_nu"""
    # period=the number of days it takes for the planet
    # to go around the sun.
    # filename is where the output will be stored
    name,e,period,a,perihelion=get_planet_data(planet)
    out_ft='{}_orbit.csv'
    outputfile=out_ft.format(planet)    
    data=dict()    
    day=list()
    X=list()
    Y=list()
    n=2*np.pi/period

    T=int(np.round(period,0)+1)
    for i in range(T):
        M=i*n
        nu=get_nu(M,e)
        r=a*(1-e*e)/(1+e*np.cos(nu))
        x=r*np.cos(nu)
        y=r*np.sin(nu)
        day.append(i+perihelion)
        X.append(x)
        Y.append(y)
    data={
        'day':day,
        'X':X,
        'Y':Y
        }
    df=pd.DataFrame(data)
    df.to_csv(outputfile)
    #draw_orb(outputfile)
    return 'Done'

def get_nu(M,e):
    """This uses 
    M=mean anomaly 
    e=eccentricity
    to calculate
    n=true anomaly
    """
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
    
#==============================================================================
# Calendar functions to calculate dates as integers. And to return those 
# integral values to standard American forme.    
#==============================================================================
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
     