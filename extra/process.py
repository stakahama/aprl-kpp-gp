###_* import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from operator import mod
import re
import os
import sys

###_* --- read in the file ---
def Readfile(filename):
    def isfloatp(x):
        try:
            float(x)
            return 1
        except:
            return 0
    with open(filename) as f:                           # open a connection to the file
    ###_ . --- read the header ---
        header = []                                     # create a container (list) for the header
        i = 1                                           # initialize loop index
        for line in f:                                  # loop through (read) line by line
            fields = line.strip().split()               # strip '\n' char and split by whitespace
            isdataline = all(map(isfloatp,fields))      # test if all values are floating pt
            if isdataline:                              # if true, we have read in a line of data
                numlines = i                            # save the number of header lines
                break                                   # break out of the loop
            header += fields                            # concatenate list of field names to "header"
            i += 1                                      # increment index by one
    ###_ . --- read the data ---
        i = 2                                           # we have already read one line of data above so start at 2
        data = []                                       # create container for the data
        dataline = map(float,fields)                    # the first data line was already read in above
        for line in f:                                  # continue reading file
            dataline += map(float,line.strip().split()) # read each line, convert to float; append to list of data fields
            i += 1                                      # increment index by 1 for having read another line of the file
            if mod(i, numlines)==0:                     # test if the remainder of i divided by number of header lines is == 0
                data.append(dataline)                   # append another row from text to data
                dataline = []                           # re-initialize dataline
                i = 1                                   # re-initialize counter/index (this is actually not necessary in this case)
    return pd.DataFrame(data,columns=header) # create data frame object


###_* --- define input file name ---

filename = sys.argv[-1]
results = Readfile(filename).drop_duplicates()
results.to_csv(os.path.basename(filename).replace('.dat','_formatted.txt'),index=False)
