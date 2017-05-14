# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 14:26:28 2017

@author: lu
"""

import pandas as pd
import numpy as np
import os
from os.path import join
 
def tmain() :
    dest = "G:\project2\NPM201507\data\\clusterNumber10_step10_ljy"
    resultPd = pd.DataFrame()
    i = 0
    for root, dirs, files in os.walk( dest ):
        for OneFileName in files :
            if OneFileName.find( '.txt' ) == -1 :
                continue
            OneFullFileName = join( root, OneFileName )
            pdData = pd.read_csv(OneFullFileName, sep = '\t') # filename
            resultPd[i] = pdData['1']
            i += 1
    label = pdData.columns[0]
    resultPd = pd.concat([pdData[label],resultPd],axis=1)
    print resultPd.shape
    #resultPd.to_csv('result10.txt', sep = '\t', columns = None, index = False, header = None)
    
    
    
if __name__ == "__main__" :
    tmain()
    
