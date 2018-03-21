# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 13:14:33 2017

@author: bwinters
"""

def prime_list(n):
    primes=[2,3]
    count=3
    while count<n:
        flag=True
        count=count+2
        for p in primes:
            if ((count%p)==0):
                flag=False
                break 
        if flag==True:
            primes.append(count) 
        if ((count%999)==0):
            print('Primes less than',count, 'calculated')
    return primes
    
def table(pList,filename):
    mFile=open(filename,'wt')
    mFile.writelines(filename+'\n')
    mFile.close()
    mFile=open(filename,'a')
    N=len(pList)
    for i in range(N):
        mFile.write(str(pList[i])+',')
        if i%10==0:
            mFile.write('\n')
        if i%10000==0:
            print('Primes less than',i, 'written')
    mFile.close()
    return 'Done'
    
def prime_table(n) :
    filename='tableOfPrimesLessThan'+str(n)+'.txt'
    p=prime_list(n)
    table(p,filename)
    return 'Done'
def load_primes():
    primes=[]
    filename='tableOfPrimesLessThan1000000.txt'   
    mFile=open(filename,'r')
    for line in mFile:
        entries=line.split(',')
        entries.pop()
        ent=[int(e) for e in entries]
        primes=primes+ent
    mFile.close()
    return primes
    
def get_pi():
    import numpy as np
    primes=load_primes()
    prod=1
    for p in primes:
        prod=prod*(1-1/(p*p))
    out=np.sqrt(6/prod)
    return out