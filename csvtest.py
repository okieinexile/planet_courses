# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 07:32:25 2017

@author: bwinters
"""

def csv_test():
    import csv
    mList=[]
    with open('primes1000.txt') as mFile:
        mReader=csv.reader(mFile)
        for row in mReader:
            row.pop()
            nrow=[int(e)/100 for e in row]
            
            nrow=[str(e) for e in nrow]
            print(nrow)
            mList.append(nrow)
    mFile.close()        
    return 'Done'
    
def panda_test():
    import pandas as pd
    Series=pd.Series
    obj=Series([.3,.4,.5,.7])
    obj.index=['a','b','c','d']
    print(obj['b'])
    return obj
    