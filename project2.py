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

def clique_v1(filename, threshold):
    #自底向上  _root_ok
    windowGraph = {}
    cliqueGraph = {}
    windowGraph[49] = nx.Graph()
    cliqueGraph[49] = nx.DiGraph()
    data = pd.read_csv(filename, index_col = 0, sep = '\t' )
    #df = data[data.columns[0]].sort_values(ascending = False)      # 默认是最小在前 若要降序 ascending = False  
    df = data[data.columns[0]].sort_values(ascending = False)
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
            #print i
            #print 'number_of_cliques(windowGraph):',nx.graph_number_of_cliques(windowGraph[49])
            for cliques in sorted(list(nx.find_cliques(windowGraph[49])), key=lambda x:len(x)) : #cliques sorted by size of each clique
                
                gene_set = set()
                term_set = set()
                if sorted(cliques) not in dic_term.values() : #[1,2,3]=[1,2,3]  [2,1,3]!=[1,2,3]
                    #this clique is new then a term is generated for this clique                     
                    cliqueGraph[49].add_node(term, annotation = cliques, windowsize = [49])# generate a term                    
                    # find child
                    #print 'cliques size:' , len(cliques)
                    for key,value in sorted(dic_term.items(), key=lambda d:d[0], reverse = True): # sorted by term id 854,853,852,851...
                        if set(value).issubset(cliques) :  
                            old_size = len(gene_set)    #old size
                            gene_set |= set(value)  #add term genes
                            if len(gene_set) > old_size :    #new size > old size
                                term_set.add(key)   # add useful term
                            if len(set(cliques).intersection(gene_set)) ==len(cliques) :   #gene_set == cliques
                                #print term, 'all link to terms',gene_set.difference(cliques)
                                for child in term_set :
                                    cliqueGraph[49].add_edge(term, child)
                                    #print term, child
                                break
                        else:
                            continue
                    if gene_set.issubset(cliques) and len(gene_set)<len(cliques):
                        #print len(gene_set), len(cliques)
                        #link to term
                        for child_term in term_set :
                            #print 'some', term, child_term
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
    print 'term:',term
    print 'Before '
    print 'windowGraph[49].size()',windowGraph[49].size(), windowGraph[49].number_of_nodes()
    print 'cliqueGraph[49].size()',cliqueGraph[49].size(), cliqueGraph[49].number_of_nodes() 
    #find terms not root but do not have a parent
    print 'find terms which is not root but do not have a parent:'
    degree_0_terms = []
    for i in cliqueGraph[49].nodes():
        if cliqueGraph[49].in_degree(i)==0:
            degree_0_terms.append(i)
    print 'degree 0 terms: ',degree_0_terms
    # link these terms to a sutible parent
    for term_id in sorted(degree_0_terms) : #从小到大 从当前向后查找
        for fp  in range(term_id+1, term ):
            if set(dic_term[term_id]).issubset(dic_term[fp]):
                cliqueGraph[49].add_edge(fp, term_id)
                break         
    '''
    #output files
    fw = open ('49ontology_edges_v3.txt', 'w')
    fw2 = open('49ontology_term_annotation_v3.txt', 'w')
    fw3 = open('49ontology_term_annotation_v3a.txt', 'w')
    fw.write('parent' +'\t' +'child' +'\n')
    fw2.write('TermId' +'\t' + ' GeneSize'+ '\t' +'Gene'+'\n')
    fw3.write('TermId' +'\t' + ' GeneSize'+ '\t' +'Gene'+'\n')
    for i in sorted(cliqueGraph[49].edges(), key=lambda d:d[0]) :
        fw.write(str(i[0]) + '\t' + str(i[1]) +'\n')
    for i in dic_term :
        fw2.write(str(i) +'\t'+ str(len(dic_term[i])) + '\t' + ','.join(dic_term[i])+ '\n')
    # input the annotation of each term by nodes attribute---annotation
    for i in sorted(cliqueGraph[49].nodes()):
        try :
            fw3.write(str(i) + '\t'+ str(len(cliqueGraph[49].node[i]['annotation'])) + '\t'+ ' '.join(cliqueGraph[49].node[i]['annotation']) +'\n')
        except :
            continue
    fw.close()
    fw2.close()
    fw3.close()
    '''
    print 'After:'
    print 'windowGraph[49].size()',windowGraph[49].size(), windowGraph[49].number_of_nodes()
    print 'cliqueGraph[49].size()',cliqueGraph[49].size(), cliqueGraph[49].number_of_nodes() 
    print 'nx.graph_clique_number(windowGraph[49]):',nx.graph_clique_number(windowGraph[49])
    print 'number_of_cliques(windowGraph):',nx.graph_number_of_cliques(windowGraph[49])
    
    
    
def clique_v0(filename, threshold):
    #自底向上 
    windowGraph = {}
    cliqueGraph = {}
    windowGraph[49] = nx.Graph()
    cliqueGraph[49] = nx.DiGraph()
    data = pd.read_csv(filename, index_col = 0, sep = '\t' )
    #df = data[data.columns[0]].sort_values(ascending = False)      # 默认是最小在前 若要降序 ascending = False  
    df = data[data.columns[0]].sort_values(ascending = False)
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
            for cliques in sorted(list(nx.find_cliques(windowGraph[49])), key=lambda x:len(x)) : #cliques sorted by size of each clique
                print 'cliques size:' , len(cliques)
                gene_set = set()
                term_set = set()
                if sorted(cliques) not in dic_term.values() : #[1,2,3]=[1,2,3]  [2,1,3]!=[1,2,3]
                    #this clique is new then a term is generated for this clique                     
                    cliqueGraph[49].add_node(term, annotation = cliques, windowsize = [49])# generate a term                    
                    # find child
                    
                    for key,value in sorted(dic_term.items(), key=lambda d:d[0], reverse = True): # sorted by term id 854,853,852,851...
                        if set(value).issubset(cliques) :  
                            old_size = len(gene_set)    #old size
                            gene_set |= set(value)  #add term genes
                            if len(gene_set) > old_size :    #new size > old size
                                term_set.add(key)   # add useful term
                            if len(set(cliques).intersection(gene_set)) ==len(cliques) :   #gene_set == cliques
                                print term, 'all link to terms',gene_set.difference(cliques)
                                for child in term_set :
                                    cliqueGraph[49].add_edge(term, child)
                                    print term, child
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
    '''
    #output files
    fw = open ('49ontology_edges_v2.txt', 'w')
    fw2 = open('49ontology_term_annotation_v2.txt', 'w')
    fw3 = open('49ontology_term_annotation_v2a.txt', 'w')
    fw.write('parent' +'\t' +'child' +'\n')
    fw2.write('TermId' +'\t' + ' GeneSize'+ '\t' +'Gene'+'\n')
    fw3.write('TermId' +'\t' + ' GeneSize'+ '\t' +'Gene'+'\n')
    for i in sorted(cliqueGraph[49].edges(), key=lambda d:d[0]) :
        fw.write(str(i[0]) + '\t' + str(i[1]) +'\n')
    for i in dic_term :
        fw2.write(str(i) +'\t'+ str(len(dic_term[i])) + '\t' + ','.join(dic_term[i])+ '\n')
    # input the annotation of each term by nodes attribute---annotation
    for i in sorted(cliqueGraph[49].nodes()):
        try :
            fw3.write(str(i) + '\t'+ str(len(cliqueGraph[49].node[i]['annotation'])) + '\t'+ ' '.join(cliqueGraph[49].node[i]['annotation']) +'\n')
        except :
            continue
    fw.close()
    fw2.close()
    fw3.close()
    '''
    print 'windowGraph[49].size()',windowGraph[49].size(), windowGraph[49].number_of_nodes()
    print 'cliqueGraph[49].size()',cliqueGraph[49].size(), cliqueGraph[49].number_of_nodes() 
    print 'nx.graph_clique_number(windowGraph[49]):',nx.graph_clique_number(windowGraph[49])
    print 'number_of_cliques(windowGraph):',nx.graph_number_of_cliques(windowGraph[49])  
                        
def clique_v2(filename, threshold):
    #自底向上       
    windowGraph = {}
    cliqueGraph = nx.DiGraph()    
    data = pd.read_csv(filename, index_col = 0, sep = '\t' )
    #df = data[data.columns[0]].sort_values(ascending = False)      # 默认是最小在前 若要降序 ascending = False  
    dic_term = {}
    window_term = {}
    term = 183
    for window in range(0, data.shape[1]) :
        windowGraph[window] = nx.Graph()
    weight_value = list(np.arange(0, 1.01, 0.01))
    #weight_value = list(np.arange(0, 1, 0.05)) # [0,0.05,0.1,...0.9,0.95]
    weight_value.sort(reverse=True)
    print min(weight_value), max(weight_value)
    for t in weight_value :        
        for window in range(0, data.shape[1]) :
            window_term[window] = {}
            df = data[data[data.columns[window]]==t][data.columns[window]] #get each column(window)            
            #df = data[(data[data.columns[window]]>=t) & (data[data.columns[window]]<=t+0.05) ][data.columns[window]]
    
            for edge in range(0, df.shape[0]) : #get each row(gene)
                node_1, node_2 = df.index[edge].split('_')
                windowGraph[window].add_edge(node_1, node_2) # generate WindowGraph
            # 当前 window 在 threshold=t 构建完成
            num_clq = 0
            #print t,window,'找到clique数目：',len(list(nx.find_cliques(windowGraph[window])))
            for clique in sorted(list(nx.find_cliques(windowGraph[window])), key=lambda x:len(x)) : #cliques sorted by size of each clique
            # 判断每一个clique   []
                num_clq +=1
                gene_set = set()
                term_set = set()                                        
                # 当前阈值下 生成的cliques
                if sorted(clique) not in window_term[window].values():
                    window_term[window][num_clq] = sorted(clique) 
                        
                if sorted(clique) not in dic_term.values() : #[1,2,3]=[1,2,3]  [2,1,3]!=[1,2,3]
                    #this clique is new then a term is generated for this clique                     
                    
                    cliqueGraph.add_node(term, annotation = clique, windowsize = [window])# generate a term                                        
                    # find child
                    #print 'cliques size:' , len(cliques)
                    for key,value in sorted(dic_term.items(), key=lambda d:d[0], reverse = True): # sorted by term id 854,853,852,851...
                        if set(value).issubset(clique) :  
                            old_size = len(gene_set)    #old size
                            gene_set |= set(value)  #add term genes
                            if len(gene_set) > old_size :    #new size > old size
                                term_set.add(key)   # add useful term
                                if len(set(clique).intersection(gene_set)) ==len(clique) :   #gene_set == cliques
                                #print term, 'all link to terms',gene_set.difference(cliques)
                                    for child in term_set :
                                        cliqueGraph.add_edge(term, child)
                                    #print term, child
                                    break
                        else:
                            continue
                    if gene_set.issubset(clique) and len(gene_set)<len(clique):
                        #print len(gene_set), len(cliques)
                        #link to term
                        for child_term in term_set :
                            #print 'some', term, child_term
                            cliqueGraph.add_edge(term, child_term)
                        # link to gene     
                        #print term,'some link to genes'
                        for child_gene in set(clique)-gene_set:
                            #print term, child_gene
                            cliqueGraph.add_edge(term, child_gene)
                    dic_term[term] = sorted(clique)
                    term = term +1         
                else :
                    if sorted(clique) in dic_term.values():
                       #print clique,window
                       te =list(dic_term.keys())[list(dic_term.values()).index(sorted(clique))]
                       #print '--term id--:', te 
                       
                       
                       if window==max(cliqueGraph.node[te]['windowsize'])+1:
                           print '原来的windowsize:',cliqueGraph.node[te]['windowsize'],'max windowsize',max(cliqueGraph.node[te]['windowsize'])
                           print '合并windowsize',[window]
                           windowsize = cliqueGraph.node[te]['windowsize'] + [window]  
                           #print windowsize
                           cliqueGraph.node[te]['windowsize']= windowsize
                    # 判断是否需要合并windowsize 与上一个window_term
                    '''
                    if sorted(clique) in window_term[window-1].values():
                       print clique,window
                       te =list(dic_term.keys())[list(dic_term.values()).index(sorted(clique))]
                       print '--term id--:', te 
                       print '原来的windowsize:',cliqueGraph.node[te]['windowsize']
                       print '新增的windowsize',[window]
                       windowsize = cliqueGraph.node[te]['windowsize'] + [window]  
                       print windowsize
                       cliqueGraph.node[te]['windowsize']= windowsize
                    '''
            
            #print 'dic_term',len(dic_term)
            #print 'windowGraph[49].size()',windowGraph.size(), windowGraph.number_of_nodes()            
            #print 'cliqueGraph[49].size()',cliqueGraph.size(), cliqueGraph.number_of_nodes()   
            print t,window,'term生成了：',len(dic_term)
            
        print '----------------'
    print 'term:',term
    print 'Before '
    #print 'windowGraph[49].size()',windowGraph.size(), windowGraph.number_of_nodes()
    print 'cliqueGraph[49].size()',cliqueGraph.size(), cliqueGraph.number_of_nodes() 
    #find terms not root but do not have a parent
    print 'find terms which is not root but do not have a parent:'
    degree_0_terms = []
    for i in cliqueGraph.nodes():
        if cliqueGraph.in_degree(i)==0:
            degree_0_terms.append(i)
    print 'degree 0 terms: ',degree_0_terms
    # link these terms to a sutible parent
    for term_id in sorted(degree_0_terms) : #从小到大 从当前向后查找
        for fp  in range(term_id+1, term ):
            if set(dic_term[term_id]).issubset(dic_term[fp]):
                cliqueGraph.add_edge(fp, term_id)
                break         
    '''
    #output files
    fw = open ('49ontology_edges_v3.txt', 'w')
    fw2 = open('49ontology_term_annotation_v3.txt', 'w')
    
    fw.write('parent' +'\t' +'child' +'\n')
    fw2.write('TermId' +'\t' + ' GeneSize'+ '\t' +'Gene'+'\n')
    fw3.write('TermId' +'\t' + ' GeneSize'+ '\t' +'Gene'+'\n')
    for i in sorted(cliqueGraph[49].edges(), key=lambda d:d[0]) :
        fw.write(str(i[0]) + '\t' + str(i[1]) +'\n')
    for i in dic_term :
        fw2.write(str(i) +'\t'+ str(len(dic_term[i])) + '\t' + ','.join(dic_term[i])+ '\n')
    fw.close()
    fw2.close()
    '''
    '''
    fw3 = open('term_annotation_v1.txt', 'w')
    # input the annotation of each term by nodes attribute---annotation
    for i in sorted(cliqueGraph.nodes()):
        try :
            fw3.write(str(i) + '\t'+cliqueGraph.node[i]['windowsize']+'\t'+ str(len(cliqueGraph.node[i]['annotation'])) + '\t'+ ' '.join(cliqueGraph.node[i]['annotation']) +'\n')
        except :
            continue
    fw3.close()
    '''
    print 'After:'
   # print 'windowGraph[49].size()',windowGraph[49].size(), windowGraph[49].number_of_nodes()
    print 'cliqueGraph[49].size()',cliqueGraph.size(), cliqueGraph.number_of_nodes() 
    #print 'nx.graph_clique_number(windowGraph[49]):',nx.graph_clique_number(windowGraph[49])
    #print 'number_of_cliques(windowGraph):',nx.graph_number_of_cliques(windowGraph[49])
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
    clique_v2('G:\project2\\NPM201507\\code\\result_c5_s10_v2_weight.txt',1)
    end = time.clock()
    print 'The function run time is : %.03f seconds' % (end-start)
    
    #f_result = 'G:\project2\NPM201507\data\\result_c5_s10_2.txt'    #The results of cluster
    #data = pd.read_table(f_result, header=None,index_col = 0) 
    #print data.head()
