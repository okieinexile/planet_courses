# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 16:44:17 2017

@author: bwinters
"""
import numpy as np
import pandas as pd
import pylab as plt
import matplotlib.pyplot as plt3
from mpl_toolkits.mplot3d import Axes3D

DATAPATH='data//'
REPORTPATH='report//'
INNERPLANETS=['mercury','venus','earth','mars']
OUTERPLANETS=['mars','jupiter','saturn','uranus','neptune']

def planet_3d(planetlist):
    fig = plt3.figure()
    ax = fig.add_subplot(111, projection='3d')

    fmt=DATAPATH+'{}_ecliptic.csv'
    for planet in planetlist:
        planet_orb(planet)
        put_in_position(planet)
        df=pd.read_csv(fmt.format(planet))
        X=list(df['x'])
        Y=list(df['y'])
        Z=list(df['z'])
        ax.plot(X,Y,Z)
    return None

def planet_graph(planetlist,listname):
    for p in planetlist:
        print(p)
        planet_orb(p)
        put_in_position(p)
        fmt=DATAPATH+'{}_ecliptic.csv'
        df=pd.read_csv(fmt.format(p))
        X=list(df['x'])
        Y=list(df['y'])
        plt.axis('equal')
        plt.plot(X,Y)
    pict_fmt=REPORTPATH+'{}.jpg'
    pictfile=pict_fmt.format(listname)
    plt.savefig(pictfile,format='jpg',dpi=1200)

def planet_table(planet):
    planet1='earth'
    planet2=planet
    planet_orb(planet1)
    planet_orb(planet2)
    put_in_position(planet1)
    put_in_position(planet2)
    relative_position(planet2,planet1)
    celestial_coords(planet)
    draw_orb(planet)
    
    return 'Done'

#==============================================================================
# This functions will take a planet that is in ecliptic coordinates
# with respect to the earth and will transform them to
# equatorial coordinates with respect to the earth. They
# will take date from <planet>_relativeto_earth_ecliptic.txt
# and will save it at <planet>_table.txt.  An intermediate file
# <planet>_relativeto_earth_equatorial.txt will be
# created.
#===============================================================================
def celestial_coords(planet):
    """This returns the file <planet>_table.csv that contains the
    celestial coordinates of the planet's orbit."""
    M=rotation(23.44,[1,0,0])
    mFrame=read_relative_position(planet)
    save_fmt=REPORTPATH+'{}_table.csv'
    savefile=save_fmt.format(planet)
    #int_fmt='{}_reltativeto_earth_equatorial.csv'
    #intmdfile=int_fmt.format(planet)
    X=np.array(mFrame['x'])
    Y=np.array(mFrame['y'])
    Z=np.array(mFrame['z'])
    V=zip(X,Y,Z)
    V=list(map(lambda x:np.matrix(x),V))
    W=np.array(map(lambda x:M*x.transpose(),V))
    date=list(map(lambda d:standard_date(d),mFrame['day']))
    ra=list()
    dc=list()
    for w in W:
        x,y,z=w[0,0],w[1,0],w[2,0]    
        rho=np.sqrt(x*x+y*y+z*z)
        x,y,z=x/rho,y/rho,z/rho
        r=np.sqrt(x*x+y*y)
        delta=np.arcsin(z)
        delta=rad_to_dms(delta)
        alpha=2*np.arctan((1-x/r)/(y/r))
        if alpha<0:
            alpha=alpha+2*np.pi        
        alpha=rad_to_hms(alpha)
        ra.append(alpha)
        dc.append(delta)
    data={   
        'date':date,
        'right_ascension':ra,
        'declination':dc  
    }
    df=pd.DataFrame(data)
    df.to_csv(savefile)    
    return 'Done'
    

def read_relative_position(planet):
    """This returns a file that containes the coordinates of planet
    with respect to earth."""
    d_fmt=DATAPATH+'{}_relativeto_earth_ecliptic.csv'
    datafile=d_fmt.format(planet)
    frame=pd.read_csv(datafile,usecols=['day','x','y','z'])
    return frame 
def rad_to_hms(rad):
    """This takes in an angle rad in radians and renders it 
    into hours, minutes, and seconds"""
    hours=12*rad/np.pi
    h=int(hours-hours%1)
    hours=(hours-h)*60
    m=int(hours-hours%1)
    hours=(hours-m)*60
    s=int(hours-hours%1)
    o_fmt='{}h {}m {}s'
    out=o_fmt.format(h,m,s)
    return out

def rad_to_dms(rad):
    """This takes in an angle rad in radians and renders it 
    into degrees, minutes, and seconds"""
    sign=np.sign(rad)
    rad=np.abs(rad)
    degrees=180*rad/np.pi
    d=int(degrees-degrees%1)
    degrees=(degrees-d)*60
    m=int(degrees-degrees%1)
    degrees=(degrees-m)*60
    s=int(degrees-degrees%1)
    o_fmt='{}d {}m {}s {}'
    if sign>0:
        direction='North'
    if sign<0:
        direction='South'
    out=o_fmt.format(d,m,s,direction)
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
    """This will give the coordinates of planet1 
    with respect to planet2. The output is in a CSV file."""
    # planet1 is to be the base planet; in most cases this will be earth
    # planet2 will be the planet whose position we are plotting
    frame1=read_planet(planet1)
    frame2=read_planet(planet2)
    days1=frame1.day.values
    days2=frame2.day.values
    begin=max([min(days1),min(days2)])
    end=min([max(days1),max(days2)])
    Flt1=((frame1.day>=begin)&(frame1.day<=end))
    F1=frame1[Flt1]
    Flt2=((frame2.day>=begin)&(frame2.day<=end))
    F2=frame2[Flt2]
    s_fmt=DATAPATH+'{}_relativeto_{}_ecliptic.csv'
    savefile=s_fmt.format(planet1,planet2)
    u=np.array(F1['x'])-np.array(F2['x'])
    v=np.array(F1['y'])-np.array(F2['y'])
    w=np.array(F1['z'])-np.array(F2['z'])
    days=[i for i in range(begin,len(u)+begin)]
    data={
        'day':days,
        'x':u,
        'y':v,
        'z':w    
    }
    pos_frame=pd.DataFrame(data)
    pos_frame.to_csv(savefile)

    return None
def read_planet(planet):
    """This returns a file that contains the coordinates of the orbit of
    planet in ecliptic coordinates."""
    f_fmt=DATAPATH+'{}_ecliptic.csv'
    frame=pd.read_csv(f_fmt.format(planet),usecols=['day','x','y','z'])

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
    """This puts the planet's orbit in its actual physical position."""
    get_matrix(planet)
    mFmt=DATAPATH+'{}_matrix.csv'
    matfile=mFmt.format(planet)
    mat=pd.read_csv(matfile,usecols=[1,2,3])
    R=np.matrix(mat)
    o_fmt=DATAPATH+'{}_orbit.csv'
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
    p_fmt=DATAPATH+'{}_ecliptic.csv'
    posfile=p_fmt.format(planet)
    pos_frame.to_csv(posfile)
    return None

def get_matrix(planet):
    """This creates the matrix that puts the planet into its
    actual physical position and saves it in a file <planet>_matrix.csv."""
    p_df=pd.read_csv(DATAPATH+'planetarydata.csv')
    filt=(p_df['name']==planet)
    df=p_df[filt]
    Omega=list(df['Omega'])[0]
    omega=list(df['omega'])[0]
    i=list(df['inclination'])[0]
    Om=(Omega/180)*np.pi
    R1=rotation(Omega,[0,0,1])
    R2=rotation(i,[np.cos(Om),np.sin(Om),0])
    R3=rotation(omega,[R2[0,2],R2[1,2],R2[2,2]])
    R=R3*R2*R1
    mFmt=DATAPATH+'{}_matrix.csv'
    outfile=mFmt.format(planet)
    mat_frame=pd.DataFrame(R)
    mat_frame.to_csv(outfile)
    return None

def rotation(theta, axisList):
    """This returns a rotation matrix of angle theta radians about
    the axis with direction numbers in axisList"""
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
    """This retrieves eccentricity, period, semimajor axis and date
    of perihelion from the file planetarydata.csv."""
    planet_data=pd.read_csv(DATAPATH+'planetarydata.csv')
    filt=(planet_data['name']==planet)
    df=planet_data[filt]
    name=list(df['name'])[0]
    e=list(df['e'])[0]
    period=list(df['period'])[0]
    a=list(df['a'])[0]
    date=ord_day(list(df['perihelion'])[0])
    return (name,e,period,a,date)
    
    
def planet_orb(planet):
    """planet is the name of a planet in lowercase.  This will output to a file
    <planet>_orbit.txt. It uses get_planet_data and get_nu"""
    # period=the number of days it takes for the planet
    # to go around the sun.
    # filename is where the output will be stored
    name,e,period,a,perihelion=get_planet_data(planet)
    out_ft=DATAPATH+'{}_orbit.csv'
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
    nu=true anomaly
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
    
    
def draw_orb(planet):
    """This draw of orbit of planet projected on the the ecliptic plane."""
    fmt=DATAPATH+'{}_ecliptic.csv'
    df=pd.read_csv(fmt.format(planet))
    X=list(df['x'])
    Y=list(df['y'])
    plt.axis('equal')
    plt.plot(X,Y)
    pict_fmt=REPORTPATH+'{}.jpg'
    pictfile=pict_fmt.format(planet)
    plt.savefig(pictfile,format='jpg')
    fmt=DATAPATH+'{}_ecliptic.csv'
    df=pd.read_csv(fmt.format(planet))
    X=list(df['x'])
    Y=list(df['y'])
    plt.axis('equal')
    plt.plot(X,Y)
    pict_fmt=REPORTPATH+'{}.jpg'
    pictfile=pict_fmt.format(planet)
    plt.savefig(pictfile,format='jpg')
    return None
    
#==============================================================================
# Calendar functions to calculate dates as integers. And to return those 
# integral values to standard American forme.    
#==============================================================================

def days_after(mdy1,mdy2):
    """Enter the days in the mm/dd/yy format and it returns the number of days 
    the second is after the first.
    """
    # enter dates mdy1 and mdy2 in the form 'mm/dd/yyyy'
    # mdy1 is assumed to be earlier than mdy2
    d1=ord_day(mdy1)
    d2=ord_day(mdy2)
    dif=d2-d1
    return dif

def ord_day(mdy):
    """ mdy is a string of the form 'mm/dd/yyyy' of the day in question
    basemdy is a string of the form 'mm/dd/yyyy' of a base day in the past
     which is after the change to Gregorian Time. We set it as '1/1/1800'
    """
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
    
    
def latex_table(planet):
    fmt=REPORTPATH+'{}_table.csv'
    df=pd.read_csv(fmt.format(planet),usecols=['date','declination','right_ascension'])
    Columns=list(df.columns)
    c_fmt='{{{}}}'
    item='{} & '
    table='\\begin{tabular}'
    line='' 
    for c in Columns:
        line=line+c_fmt.format('c')
    table=table+line+'\n '
    line=''
    for c in Columns:
        line=line+item.format(c)
    line=line[:-2]+'\\\\ \n '
    N=len(list(df['date']))
    table=table+line
    for i in range(N):
        line=''
        for j, c in enumerate(Columns):
            line=line+item.format(df.iat[i,j])
        line=line[:-2] + '\\\\ \n '
        table=table+line
    table=table[:-4]+'\n '
    table=table + '\end{tabular}'
    out_fmt='{}_table.tex'
    ofile=REPORTPATH+out_fmt.format(planet)
    with open(ofile,'w') as f:
        f.write(table)
    return  None 


    
    