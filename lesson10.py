# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 16:35:37 2017

@author: bwinters
"""

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
    savefile='coords_'+planet2+'_relativeto_'+planet1+'.txt'
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
    import pandas
    import csv
    file=open(planet+'_position.txt')
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
            
       