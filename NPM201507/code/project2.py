# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 10:05:46 2017

@author: lu
"""

import pandas as pd
import numpy as np
import networkx as nx
import time
import os
from os.path import join

 
def cluster_NPM(dest, fwname) :
    resultPd = pd.DataFrame()
    i = 1
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
    resultPd = resultPd.set_index([label])
    
    resultPd.to_csv(fwname, sep = '\t', columns = None, index = True, header = None)
    return resultPd
    
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
        #fwname ='G:\project2\\NPM201507\\data\windowMatrix\matrix'+os.path.split(root)[1]+'.txt' 
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
                    #name = (data.index[i], data.index[j])
                    name = data.index[i]+'_'+data.index[j]
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
    print len(dic), len(dicEdge)
    return dic, dicNum, dicEdge
    
# remove duplicate edges    
def graph_edges(dicEdge,edgefile):
    G = nx.DiGraph()
    print G.size(), G.number_of_nodes()
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

# for cytoscape
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

def clique(filename, threshold):
    #自顶向下
    data = pd.read_csv(filename, index_col = 0, sep = '\t' )
    print data.shape[1] #54 columns
    print data.columns[0]
    #print data[data.columns[data.shape[1]-1]].head(5)
    windowGraph = {}
    cliqueGraph = {}
    for i in range(0, data.shape[1]):
        # for each column(window)
        windowGraph[i] = nx.Graph()
        df = data[data[data.columns[i]]>=threshold][data.columns[i]]
        for edge in range(0, df.shape[0]):
            # for each edges up the threshold
            node_1, node_2 = df.index[edge].split('_')
            windowGraph[i].add_edge(node_1, node_2)
        print 'window', i, windowGraph[i].size(), windowGraph[i].number_of_nodes()   
        
        cliques = list(nx.find_cliques(windowGraph[i]))
        N = max(len(c) for c in cliques)
        print 'size of max clique:',N
        for c in sorted(cliques):
            if len(c) == N :
                #print c # nodes list
                cliqueGraph[i] = c

def clique2(filename, threshold):
    #自底向上
    '''
    data = pd.read_csv(filename, index_col = 0, sep = '\t' )
    windowGraph = {}    # generate window graph
    cliqueGraph = {}    # generate term
    
    for i in range(0, data.shape[1]):
        # Declare windowGraph and cliqueGraph for each column(window)
        windowGraph[i] = nx.Graph()
        cliqueGraph[i] = nx.Graph()
        term = 1
        lastIndex = []
        while threshold >=0:
            df = data[data[data.columns[i]]>=threshold][data.columns[i]]
            newIndex = df.index.difference(lastIndex)
            
            #Update windowGraph
            for edge in range(0, len(newIndex)):
                node_1, node_2 = newIndex[edge].split('_')
                windowGraph[i].add_edge(node_1, node_2)
            #Update cliqueGraph
            for cliques in sorted(list(nx.find_cliques(windowGraph[i]))) :
                for node in cliques :
                    cliqueGraph[i].add_edge(term, node)
                term = term +1
            lastIndex = df.index
            threshold = threshold-0.2
        print 'window', i, windowGraph[i].size(), windowGraph[i].number_of_nodes()   
     '''   
    windowGraph = {}
    cliqueGraph = {}
    windowGraph[49] = nx.Graph()
    cliqueGraph[49] = nx.DiGraph()
    data = pd.read_csv(filename, index_col = 0, sep = '\t' )
    df = data[data.columns[0]].sort_values(ascending = False)      # 默认是最小在前 若要降序 ascending = False  
    t=1
    term = 183
    dic_term = {}
    #dic_cliques  = {}

    for i in range(0, df.shape[0]):

        if df[i] == t :
            node_1, node_2 = df.index[i].split('_')
            windowGraph[49].add_edge(node_1, node_2)
        else :
            # find cliques when threshold = t 
            print i
            print 'number_of_cliques(windowGraph):',nx.graph_number_of_cliques(windowGraph[49])
            for cliques in sorted(list(nx.find_cliques(windowGraph[49]))) : 
                print 'cliques length' , len(cliques)
                gene_set = set()
                term_set = set()
                if sorted(cliques) not in dic_term.values() : 
                    #this clique is new                     
                    cliqueGraph[49].add_node(term, annotation = cliques, windowsize = [49])# generate a term                    
                    # find child
                    
                    for key,value in sorted(dic_term.items(), key=lambda d:d[0], reverse = True):
                        if set(value).issubset(cliques) :  
                            old_size = len(gene_set)    #old size
                            gene_set |= set(value)  #add term genes
                            if len(gene_set) > old_size :    #new size > old size
                                term_set.add(key)   # add useful term
                            if len(set(cliques).intersection(gene_set)) ==len(cliques) :   #gene_set == cliques
                                print term, 'all link to terms',gene_set.difference(cliques)
                                for child in term_set :
                                    cliqueGraph[49].add_edge(term, child)
                                    print term, child, len(gene_set), len(cliques)
                                break
                        else:
                            continue
                    if gene_set.issubset(cliques) and len(gene_set)<len(cliques):
                        #print len(gene_set), len(cliques)
                        #link to term
                        for child_term in term_set :
                            print 'some', term, child_term
                            cliqueGraph[49].add_edge(term, child_term)
                        # link to gene     
                        #print term,'some link to genes'
                        for child_gene in set(cliques)-gene_set:
                            #print term, child_gene
                            cliqueGraph[49].add_edge(term, child_gene)
                    dic_term[term] = sorted(cliques)
                    term = term +1                                       
                else :
                    continue
            t = df[i]
            if not t==-1:
                node_1, node_2 = df.index[i].split('_')
                windowGraph[49].add_edge(node_1, node_2)
            
            print 'dic_term',len(dic_term)
            print 'windowGraph[49].size()',windowGraph[49].size(), windowGraph[49].number_of_nodes()
            
            print 'cliqueGraph[49].size()',cliqueGraph[49].size(), cliqueGraph[49].number_of_nodes()            

    #fw = open ('49ontology_edges.txt', 'w')
    #for i in sorted(cliqueGraph[49].edges(), key=lambda d:d[0]) :
        #fw.write(str(i[0]) + '\t' + str(i[1]) +'\n')
        
    #fw2 = open('49ontology_term_annotation.txt', 'w')
    #for i in dic_term :
        #fw2.write(str(i) +'\t'+ str(len(dic_term[i])) + '\t' + ','.join(dic_term[i])+ '\n')
    
    print 'windowGraph[49].size()',windowGraph[49].size(), windowGraph[49].number_of_nodes()
    print 'cliqueGraph[49].size()',cliqueGraph[49].size(), cliqueGraph[49].number_of_nodes()      
    print 'nx.graph_clique_number(windowGraph[49]):',nx.graph_clique_number(windowGraph[49])
    print 'number_of_cliques(windowGraph):',nx.graph_number_of_cliques(windowGraph[49])
               
    '''
    data = pd.read_csv(filename, index_col = 0, sep = '\t' )
    df = data[data[data.columns[5]]>=0.7][data.columns[5]]
    df2 =data[data[data.columns[5]]>=0.8][data.columns[5]]
    windowGraph = {}
    cliqueGraph = {}
    G = nx.Graph()
    windowGraph[5] = nx.Graph()
    lastIndex = []
    
       
    term = 0
    newIndex = df2.index.difference(lastIndex)
    print len(lastIndex),len(df2.index),len(newIndex)
    for edge in range(0, len(newIndex)):
            # for each edges up the threshold
        node_1, node_2 = newIndex[edge].split('_')
        windowGraph[5].add_edge(node_1, node_2)
        

    for i in sorted(list(nx.find_cliques(windowGraph[5]))) :
            if len(i)>=2 : 
                for node in i :
                    G.add_edge(term, node)
                term = term +1
                print i
                
                
    print 'window', 5, windowGraph[5].size(), windowGraph[5].number_of_nodes()   
    cliques = list(nx.find_cliques(windowGraph[5]))
    N = max(len(c) for c in cliques)
    print 'size of max clique:',N
    for c in sorted(cliques):
        if len(c) == N :
                #print c # nodes list
            cliqueGraph[5] = c     
    lastIndex = df2.index
    
    
    newIndex = df.index.difference(lastIndex)
    print len(lastIndex),len(df.index),len(newIndex)
    for edge in range(0, len(newIndex)):
            # for each edges up the threshold
        node_1, node_2 = newIndex[edge].split('_')
        windowGraph[5].add_edge(node_1, node_2)
        G.add_edge(term, node_1)
        G.add_edge(term, node_2)
        term+=1
    print 'window', 5, windowGraph[5].size(), windowGraph[5].number_of_nodes()   
    cliques = list(nx.find_cliques(windowGraph[5]))
    N = max(len(c) for c in cliques)
    print 'size of max clique:',N
    for c in sorted(cliques):
        if len(c) == N :
                #print c # nodes list
            cliqueGraph[5] = c
    '''             
    
if __name__ == '__main__':   
    start = time.clock()
    
    # Integrate cluster files into a matrix. result_clustNum_step.txt The index is gene name
    fw_clust = 'result_c5_s10_v2'
    '''
    data = cluster_NPM("G:\project2\NPM201507\data\\clusterNumber10_step15_ljy", fw_clust+'.txt')
    #Group by window
    startWindow = 1
    windowSize = 3
    minGeneSize = 4
    fw_group = fw_clust+'_group'   
    dic, dicNum, dicEdge = classifyByWindow(startWindow, windowSize, minGeneSize, fw_group+'.txt')
    #Remove the duplicate edges
    graph_edges(dicEdge,fw_group+'.txt')
    '''
    #Integrate cluter files, each NPM window is calculate for 100 times.
    #cluster_NPM100('G:\project2\\NPM201507\\clusterResult')
    #calculate weight for each window of each gene pattern    
    fw_weight = fw_clust +'_weight'
    #cal_weight('G:\project2\\NPM201507\\data\\100_c5_s10_windowMatrix',fw_weight +'.txt')
    
    #find cliques for each window graph using weight.txt
    #clique('G:\project2\\NPM201507\\code\\result_c5_s10_v2_weight.txt',0.95)
    clique2('G:\project2\\NPM201507\\code\\result_c5_s10_v2_weight.txt',1)
    end = time.clock()
    print 'The function run time is : %.03f seconds' % (end-start)
    
    #f_result = 'G:\project2\NPM201507\data\\result_c5_s10_2.txt'    #The results of cluster
    #data = pd.read_table(f_result, header=None,index_col = 0) 
    #print data.head()
