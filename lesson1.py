# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 08:39:50 2017

@author: bwinters
"""

def save_list(pList,filename):
    # This will save a list of number in a file
    
    # Create a File to write to.  This will have the name contained
    # in the variable filename
    mFile=open(filename,'wt')
    
    
    # The variable count will keep track of the numerical order of this list
    count=0
    
    # We will do through pList one item at a time
    for item in pList:
        # Turn the item into a string and write it to the file
        mFile.write(str(item))
        
        # Put in a common as a separater
        mFile.write(',')
        
        # increment our counter
        count=count+1
        
        # Every tenth item create a new line
        if count%10==0:
            mFile.write('\n')
    # When we are out of the loop close the file.
    mFile.close()
    
    # Send the message to the user that we are done.
    return 'Saved '+filename

def read_list(filename):
    # This will read a list of numbers that have been store in the textfile
    # whose name is store in filename
    
    # Create list to store our values
    out=[]
    
    # Open the file
    mFile=open(filename,'r')
    
    # Read the file line by line
    for line in mFile:
        
        # split the line along the commas
        data=line.split(',')
        
        # Remove the last entry of the line which would be a new line character
        data.pop()
        
        # change each entry in data to an integer
        data=[float(item) for item in data]
        
        # Add our data to out
        out=out+data
    # Close the file
    mFile.close()
    
    # Send the data to the user
    return out