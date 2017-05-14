# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 10:05:46 2017

@author: lu
"""

import pandas as pd
import numpy as np
import networkx as nx
import time


def classifyByWindow(startWindow, windowSize, minGeneSize, filename):    
    dic = {}
    fw = open(filename,'w') # rem
    fw.write('wid' + '\t' + 'startWindow' + '\t' + 'stopWindow' + '\t' + 'windowSize' + '\t' + 'geneSize' + '\t' + 'groupGenes' + '\n' )
    for startWindow in range(startWindow,48):
        for stopWindow in range (startWindow + windowSize, 51):
        #e.g.startWindow = 1 stopWindow = 3, Then the window is [1,2]
            group1 = data.groupby(list(np.arange(startWindow,stopWindow)))
            for name, group in group1 :
                if group.index.size >  minGeneSize-1 :   
                    if not dic.has_key(frozenset(group.index)) :
                        dic[frozenset(group.index)] = (startWindow, stopWindow)
                    else:
                        if dic[frozenset(group.index)][0] >= startWindow and dic[frozenset(group.index)][1] <= stopWindow :
                            dic[frozenset(group.index)] = (startWindow, stopWindow)
    #Encording wid    
    dicNum = {}
    dicFlag = {}
    for i in dic :
        if dic[i][1] - dic[i][0]  <10 :
            wid = dic[i][0] * 100 + dic[i][1] - dic[i][0]
        else:
            wid = dic[i][0] * 1000 + dic[i][1] - dic[i][0]
        if not wid in dicNum.values():
            dicNum[i] = wid
            dicFlag[wid] = 'a'
        else :
            dicNum[i] = str(wid) + dicFlag[wid]
            dicFlag[wid] = chr(ord(dicFlag[wid]) + 1)            
        fw.write(str(dicNum[i])+'\t'+str(dic[i][0]) + '\t' + str(dic[i][1]) + '\t' + str(dic[i][1] - dic[i][0]) +'\t' + str(len(list(i))) + '\t' + str(list(i)) + '\n' )
    fw.close()
    # Finding edges
    dicEdge = {}
    for geneList1 in dic :       
        for geneList2 in dic : 
            if dicNum[geneList1]!= dicNum[geneList2]:           
                if geneList1.issubset(geneList2):
                    if dic[geneList1][0]<= dic[geneList2][0] and dic[geneList1][1]>= dic[geneList2][1]: #1) geneList1 is a sub of geneList2
                        dicEdge[(dicNum[geneList2],dicNum[geneList1])] = 1 
                        #fw2.write(str(dicNum[geneList2])+ '\t' + str(dicNum[geneList1]) +'\n')
            else:
                continue    
    return dic, dicNum, dicEdge
    
# remove duplicate edges    
def graph_edges(dicEdge,edgefile):
    G = nx.DiGraph()
    for key, value in dicEdge.items():
        G.add_edge(key[0],key[1])
    print 'Before removing the duplicate edges, the G.size is :',G.size(),G.number_of_nodes()
    remEdges = []
    for i in G.nodes():        
        predec = G.predecessors(i)       
        if len(predec) >1:
            for  k in range(0, len(predec)-1): 
                for j in range (k+1,len(predec)):
                    if nx.has_path(G,predec[k],predec[j]):
                        remEdges.append((predec[k],i))
                    elif nx.has_path(G,predec[j],predec[k]):
                        remEdges.append((predec[j],i))                       
    G.remove_edges_from(remEdges)
    print 'Atrer removing the duplicate edges, the G.size is:',G.size(),G.number_of_nodes()
    remoDupEdges = edgefile.split('.txt')[0]+'_Edges.txt'
    fw3 = open(remoDupEdges,'w')
    fw3.write('parent'+'\t'+'child'+'\n')
    for i in G.edges():
        fw3.write(str(i[0])+'\t'+str(i[1])+'\n')
    fw3.close() 

def chose_colums():
    filename = 'classResult1_3_4_v0.txt'
    data2 = pd.read_table(filename, sep = '\t')
    print data2.head(5)
    data3 = data2.loc[:,['wid','windowSize']]
    print data3.head(5)
    data3.to_csv('windowSize.txt',sep = '\t',index = False)
    data4 = data2.loc[:,['wid','geneSize']]
    print data4.head(5)
    data4.to_csv('geneSize.txt',sep = '\t',index = False)
      
if __name__ == '__main__':
    
    start = time.clock()
    f_result = 'G:\project2\NPM201507\data\\result10.txt'    #The results of cluster
    data = pd.read_table(f_result, header=None,index_col = 0)   #The first column is the index of each row
    startWindow = 1
    windowSize = 3
    minGeneSize = 4
    fw_filename = 'classResult_10_1_3_4_v0.txt'    
    dic, dicNum, dicEdge = classifyByWindow(startWindow, windowSize, minGeneSize, fw_filename)
    graph_edges(dicEdge,fw_filename)
    end = time.clock()
    print 'The function run time is : %.03f seconds' % (end-start)
    
    


