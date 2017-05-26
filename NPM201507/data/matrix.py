# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 14:26:28 2017

@author: lu
"""

import pandas as pd
import numpy as np
import os
from os.path import join
import time
 
def cluster_NPM(dest, fwname) :
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
    resultPd.to_csv(fwname, sep = '\t', columns = None, index = False, header = None)
    
    
def cluster_NPM100(path1) :
    j=49  
    x = pd.DataFrame()
    for root, dirs, files in os.walk (path1) :
        resultPd = pd.DataFrame()
        i=0  
        if root==path1 :
            continue
        for filepath in files:
            if filepath.find( '.txt' ) == -1 :
                continue            
            fns = os.path.join(root,filepath)
            pdData = pd.read_csv(fns, sep = '\t')   
            resultPd[i] = pdData['1']   
            i += 1
            label = pdData.columns[0]
            x = pdData[label]
        resultPd = pd.concat([x,resultPd],axis=1)
        print resultPd.shape
        print root, os.path.split(root)[1]
        fwname ='G:\project2\\NPM201507\\data\windowMatrix\matrix'+os.path.split(root)[1]+'.txt' 
        resultPd.to_csv(fwname, sep='\t', columns = None, index = False, header = None)
        j=j+1    
       
        
def cal_weight(dest,fwname):
    
    resultPd = pd.DataFrame()
    i = 0
    window = 49
    for root, dirs, files in os.walk( dest ):
        for OneFileName in files :
            dic = {}
            if OneFileName.find( '.txt' ) == -1 :
                continue
            filename = join( root, OneFileName )
            data = pd.read_csv(filename, sep = '\t',index_col=0, header=None) # #load the matrix49, matrix50...
            for i in range(0, data.shape[0]-1) :
                for j in range(i+1, data.shape[0]):
                    weight = np.count_nonzero(data.iloc[i]==data.iloc[j])/float(data.shape[1])
                    name = (data.index[i], data.index[j])
                    dic[name] = weight
            print 'len(dic):', len(dic)
            df = pd.DataFrame(dic.items(),columns = ['gene',window])
            resultPd[window] = df[window]
            window = window+1
            print data.shape,df.shape,resultPd.shape
    label = df['gene']
    resultPd = pd.concat([label,resultPd],axis=1)
    print resultPd.shape
    resultPd.to_csv(fwname, sep = '\t',  index = False)
    '''
    data = pd.read_table(dest, sep = '\t',index_col=0, header=None)      
    dic = {}

    w=49
    m=0
    n=0
    fw = open('shao.txt','w')
    for i in range(0, data.shape[0]-1) :
        for j in range(i+1, data.shape[0]):
            n=n+1
            fw.write(str(i)+'\t'+str(j)+'\n')
            weight = np.count_nonzero(data.iloc[i]==data.iloc[j])/float(data.shape[1])
            name = (data.index[i], data.index[j])
            if dic.has_key(name) :
                m=m+1 
                print data.index[i], data.index[j]
            else :
                dic[name] = weight
    print m,n
    df = pd.DataFrame(dic.items(),columns = ['gene',w])
    print dic[('AT2G04030','AT5G62720')]
    print df.shape
    #df.to_csv('df.txt', sep='\t', index = False)
    '''
            
    
if __name__ == "__main__" :
    start = time.clock()
    #cluster_NPM("G:\project2\NPM201507\data\\clusterNumber5_step10_ljy", 'result_c5_s10.txt')
    #cluster_NPM100('G:\project2\\NPM201507\\clusterResult')
    cal_weight('G:\project2\NPM201507\data\windowMatrix','weight.txt')
    end = time.clock()
    print 'The function run time is : %.03f seconds' % (end-start)
    
    
    
    
    
    '''
    print data.tail()
    weight =  np.count_nonzero(data.iloc[0]==data.iloc[data.shape[0]-1])/float(data.shape[1])
    print data.index[0]
    name = (data.index[0], data.index[data.shape[0]-1])
    dic = {}
    print type(name)
    print data.iloc[0].head(),data.index[0]
    print data.iloc[181].head(),data.index[181]
    dic[name] = weight 
    if 'AT2G04030' in name :
        print '2'
    '''
    #w = pd.DataFrame({'x':x,'y':y})
    #print w[x==y]
    #print w[x==y].shape[0]