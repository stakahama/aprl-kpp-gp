# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 15:01:22 2014

@author: shipley
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def matmul(filename1,filename2):
    def isfloat(x):
        try:
            float(x)
            return 1
        except:
            return 0 
  
    def convertmat(filename):
        data = []
        header = []
        with open(filename) as f:
            i = 1
            for line in f:
                nums = line.strip().split(',')
                isdataline = isfloat(nums[1])
                if isdataline:
                    if filename == filename2:
                        header.append(nums[0])
                    nums.pop(0)
                    nums = map(float, nums)
                    data.append(nums)
                elif filename == filename1:
                    header = nums
                    header.pop(0)
                i += 1

        return header, data
    
    results1 = convertmat(filename1)
    header1, data1 = results1
    
    results2 = convertmat(filename2)
    header2, data2 = results2
    
    #print header2

    if header2 != header1:    
    #if [x for x in header2 if x not in header1] != []:
        print 'DANGER!!! YOU MESSED UP SOMEWHERE'
        
    return np.dot(data1,data2)


if __name__=='__main__':    
    results1 = matmul('toluene_edited.txt','r2_Groups.txt')
    #results2 = matmul('test_aer.txt','COCount.txt')
    total = results1 #+ results2
    
    with open('xxx.txt','w') as r:
        for line in total:
            r.write(str(line)+'\n')
    
    final = pd.DataFrame(total)
    plt.ion()
    plt.gcf().clear()
    # plt.plot(final.index,final['n_C'])
    plt.plot(final.index, final)
    plt.show()   
