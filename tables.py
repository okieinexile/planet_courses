# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 13:50:28 2017

@author: bwinters
"""
def make_list(n, function):
    p=[]
    for i in range(n):
        p.append(function(i))
    return p

def make_table(pList,filename):
    mfile=open(filename,'wt')
    for p in pList:
        x=str(p)
        line=x+'\n'
        mfile.write(line)
    mfile.close()
    return 'Done'
    
def xx(i):
    return i*i