# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 13:50:02 2017

@author: bwinters
"""
def put_in_position(planet):
    import numpy as np
    import csv
    matfile=planet+'_matrix.txt'
    mfile=open(matfile,'rt')
    reader=csv.reader(mfile)
    mat=[]
    for line in reader:
        if line!=[]:  
            x,y,z=line  
            mat.append([float(x),float(y),float(z)])
    mfile.close()
    R=np.matrix(mat)

    orbfile=planet+'_orbit.txt'
    ofile=open(orbfile,'rt')
    reader=csv.reader(ofile)
    posfile=planet+'_position.txt'
    pfile=open(posfile,'wt')
    writer=csv.writer(pfile)
    for line in reader:
        if line!=[]:
            day,x,y=line
            x=float(x)
            y=float(y)
            v=np.matrix([x,y,0])
            w=R*v.T
            writer.writerow([day,w[0,0],w[1,0],w[2,0]])
    ofile.close()
    pfile.close()

def get_matrix(filename):
    import csv
    import numpy as np
    myFile=open(filename,'rt')
    reader=csv.reader(myFile)
    for line in reader:
        name,Omega,i,omega=line    
        Omega=float(Omega)
        i=float(i)
        omega=float(omega)
    myFile.close()
    R1=rotation(Omega,[0,0,1])
    R2=rotation(i,[np.cos(Omega),np.sin(Omega),0])
    R3=rotation(omega,[R2[0,2],R2[1,2],R2[2,2]])
    R=R3*R2*R1
    matFile=open(name+'_matrix.txt','wt')
    writer=csv.writer(matFile)
    for i in range(3):
        line=[]
        for j in range(3):
            line.append(R[i,j])
        writer.writerow(line)
    matFile.close()
    return(R)

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